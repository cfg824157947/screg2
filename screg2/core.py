from sklearn.decomposition import NMF
import torch
import torch.sparse as sparse
import numpy as np
from scipy.sparse import issparse

def sparse_mean(sparse_matrix):
    """Calculate the mean of non-zero elements of a sparse matrix."""
    # Using .data to access non-zero elements directly
    # This only considers non-zero elements in the mean calculation
    return sparse_matrix.data.mean() if sparse_matrix.nnz else 0



def defaultpar_CoupledNMF_default(PeakO, W1, X, W2, Reg = None, W3 = None, beta=1, arfa=1, withoutRE=True):
    lambda1 = beta * (PeakO @ W1).mean() / (X @ W2).mean()
    if withoutRE:
        lambda2 = 0
    else:
        if Reg is None or W3 is None:
            raise ValueError("Reg and W3 must be provided if withoutRE is False")
        lambda2 = arfa * (PeakO @ W1).mean() / (W3 @ Reg.T).mean()
    return lambda1, lambda2

def perform_nmf(data, n_components, max_iter=200):
    """
    Perform Non-negative Matrix Factorization (NMF) on the given data.

    Parameters:
    - data: 2D array-like data matrix (e.g., numpy array or pandas DataFrame)
    - n_components: Number of components to extract
    - max_iter: Maximum number of iterations (default is 200)

    Returns:
    - W: Matrix of basis vectors
    - H: Matrix of coefficients
    """
    print(data.shape)
    nmf_model = NMF(n_components=n_components, max_iter=max_iter)
    W = nmf_model.fit_transform(data)
    H = nmf_model.components_
    return W, H





def nmf_cluster_joint_cross_domain(PeakO, X, lambda1, W10, W20, H0, Reg=None, K=100, maxiter=100, lambda2=None, W30=None, c1=None, c2=None, Reg_w=None, core=8, withoutRE=True):
    np.set_printoptions(precision=6)
    np.set_printoptions(suppress=True)
    np.set_printoptions(threshold=6)
    np.set_printoptions(linewidth=200)
    np.set_printoptions(formatter={'float': '{: 0.6f}'.format})

    tolfun = 1e-6
    dnorm, dnorm0 = 0, 0

    if (withoutRE == True):
        dnorm0 = np.square(PeakO - (W10 @ H0)).sum() + lambda1 * np.square(X - (W20 @ H0)).sum()
        for iter in range(maxiter):
        #W10T, W20T, W30T, H0T = W10.T, W20.T, W30.T, H0.T
            numer = (W10.T @ PeakO) + (lambda1 * (W20.T @ X)) 

            H0 = H0 * (numer /((((W10.T @ W10) + (lambda1 * W20.T @ W20)  + np.finfo(np.float64).eps)) @ H0))
            HH = H0 @ H0.T
            numer = PeakO @ H0.T
            W10 = (W10 * (numer / (W10 @ HH + np.finfo(np.float64).eps)))
            numer = X @ H0.T
            W20 = (W20 * (numer / (W20 @ HH + np.finfo(np.float64).eps)))
            if iter % 20 == 0:
                dnorm = np.square(PeakO - (W10 @ H0)).sum() + lambda1 * np.square(X - (W20 @ H0)).sum()
                if iter < maxiter:
                    if (dnorm0 - dnorm) <= tolfun or (dnorm0 - dnorm) <= 0:
                        print("dnorm0-dnorm", dnorm0 - dnorm, "is small")
                        break
                    elif iter == maxiter:
                        break
                    dnorm0 = dnorm
                print("iteration", iter)

        return {'W1': W10, 'W2': W20, 'H': H0}

    else:
        print("I have not maken withRE part yet")
        return False
        if False: 
            dnorm0 = np.square(PeakO - (W10 @ H0.T)).sum() + lambda1 * np.square(X - (W20 @ H0)).sum() + lambda2 * np.square(Reg - (W30 @ H0)).sum()
            for iter in range(maxiter):
                #W10T, W20T, W30T, H0T = W10.T, W20.T, W30.T, H0.T
                numerO = (W10.T @ PeakO) + (lambda1 * (W20.T @ X)) + (lambda2 * (W30.T @ Reg))

                H = ifzeroMCPP((H0 * (numerO / (((W10.T @ W10) + (lambda1 * W20.T @ W20) + (lambda2 * W30.T @ W30)) @ H0 + epsMCpp(numerO))))).copy()

               # HT = H.T
                HH = H @ HT

                #numerO = PeakO @ HT
                numerO = PeakO @ H.T
                W1 = ifzeroMCPP((W10 * (numerO / (W10 @ HH + epsMCpp(numerO))))).copy()

                #numerX = X @ HT
                numerX = X @ H.T
                W2 = ifzeroMCPP((W20 * (numerX / (W20 @ HH + epsMCpp(numerX))))).copy()

                Tmp1 = chooesVinMCPP(W1, c2, 0)
                Tmp2 = chooesVinMCPP(W2, c1, 0)
                numerR = Tmp1 + Tmp2
                numerR = CppoperationMA_demo(numerR, Reg_w, 2)
                Tmp1 = Reg @ H.T
                numerR += Tmp1

                W3 = ifzeroMCPP((W30 * (numerR / (W30 @ HH + W30 + epsMCpp(numerR))))).copy()

                if iter > 300:
                    if (dnorm0 - dnorm) <= tolfun or (dnorm0 - dnorm) <= dnorm0:
                        print("dnorm0-dnorm", dnorm0 - dnorm, "is small")
                        break
                    elif iter == maxiter:
                        break

                W10 = W1.copy()
                H0 = H.copy()
                W20 = W2.copy()
                W30 = W3.copy()

                if iter % 20 == 0:
                    Tmp1 = -(W1 @ H)
                    Tmp1 += PeakO
                    dnorm = np.square(Tmp1).sum()
                    Tmp1 = -(W2 @ H)
                    Tmp1 += X
                    dnorm += lambda1 * np.square(Tmp1).sum()
                    Tmp1 = -(W3 @ H)
                    Tmp1 += Reg
                    dnorm += lambda2 * np.square(Tmp1).sum()

                    print("iteration", iter)

            return {'W1': W1, 'W2': W2, 'W3': W3, 'H': H}


# PeakO and X is csr sparse matrix
def nmf_cluster_joint_cross_domain_torch_sparse_mb(PeakO, X, lambda1, W10, W20, H0, Reg=None, K=100, maxiter=40, lambda2=None, W30=None, c1=None, c2=None, Reg_w=None, core=8, withoutRE=True):
    """
    Perform Non-negative Matrix Factorization (NMF) on joint cross-domain data using PyTorch.

    Parameters:
    - PeakO (csr_matrix): Sparse matrix representing data from domain 1.
    - X (csr_matrix): Sparse matrix representing data from domain 2.
    - lambda1 (float): Regularization parameter.
    - W10 (numpy.ndarray): Initial factor matrix for domain 1.
    - W20 (numpy.ndarray): Initial factor matrix for domain 2.
    - H0 (numpy.ndarray): Initial coefficient matrix.
    - Reg (callable): Regularization function (optional).
    - K (int): Number of components for NMF.
    - maxiter (int): Maximum number of iterations.
    - lambda2 (float): Optional regularization parameter for additional regularization.
    - W30 (numpy.ndarray): Optional additional factor matrix.
    - c1 (float): Optional constant.
    - c2 (float): Optional constant.
    - Reg_w (callable): Optional regularization function for factor matrices.
    - core (int): Number of CPU cores to use.
    - withoutRE (bool): Flag indicating whether to perform with or without regularization.

    Returns:
    - dict: Dictionary containing factor matrices W1, W2, and coefficient matrix H.
    """
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    tolfun = 1e-6
    dnorm, dnorm0 = 0, 0

    indices = torch.randperm(PeakO.shape[1])

    if withoutRE:
        W10 = torch.tensor(W10, device=device)
        W20 = torch.tensor(W20, device=device)
        H0 = torch.tensor(H0, device=device)
        mb_size = mb_size = int(2 ** np.log10(PeakO.shape[1]) * 128)# Adjust this according to your preference
        n_batch = (PeakO.shape[1] + mb_size - 1) // mb_size

       # dnorm0 = ((PeakO - torch.matmul(W10, H0)).pow(2).sum() + lambda1 * (X - torch.matmul(W20, H0)).pow(2).sum()).item()

        for iter in range(maxiter):
            for i in range(n_batch):
                print("n_batch: ", n_batch)
                start_idx = i * mb_size
                end_idx = min((i + 1) * mb_size, PeakO.shape[1])
                mb_indices = indices[start_idx:end_idx]
                
                PeakOp = torch.tensor(PeakO[:, mb_indices].toarray())
                Xp = torch.tensor(X[:, mb_indices].toarray())

                print("after trans", n_batch)


                numer = torch.matmul(W10.t(), PeakOp) + lambda1 * torch.matmul(W20.t(), Xp)

                H0[:, mb_indices] *= (numer / (torch.matmul(torch.matmul(W10.t(), W10) + lambda1 * torch.matmul(W20.t(), W20), H0[:, mb_indices]) + torch.finfo(torch.float64).eps))
                HH = torch.matmul(H0[:, mb_indices], H0[:, mb_indices].t())
                numer = torch.matmul(PeakOp, H0[:, mb_indices].t())
                W10 *= (numer / (torch.matmul(W10, HH) + torch.finfo(torch.float64).eps))
                numer = torch.matmul(Xp, H0[:, mb_indices].t())
                W20 *= (numer / (torch.matmul(W20, HH) + torch.finfo(torch.float64).eps))

            #if iter % 20 == 0:
#                dnorm = ((PeakO - torch.matmul(W10, H0)).pow(2).sum() + lambda1 * (X - torch.matmul(W20, H0)).pow(2).sum()).item()
             #   if iter < maxiter:
              #      if (dnorm0 - dnorm) <= tolfun or (dnorm0 - dnorm) <= 0:
               #         print("dnorm0-dnorm", dnorm0 - dnorm, "is small")
                #        break
                 #   elif iter == maxiter:
                  #      break
                   # dnorm0 = dnorm
               # print("iteration", iter)

        return {'W1': W10.cpu().numpy(), 'W2': W20.cpu().numpy(), 'H': H0.cpu().numpy()}

