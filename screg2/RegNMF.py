import numpy as np
import pandas as pd
import scanpy as sc
from . import core
from scipy.sparse import csc_matrix, csr_matrix, find
import psutil
import os 
from anndata import AnnData
from scipy import sparse

def sparse_row_median(csr_mat):
    """
    Compute the median of each row in a sparse CSR matrix.

    Parameters:
    - csr_mat (csr_matrix): Input sparse CSR matrix.

    Returns:
    - numpy.ndarray: Array containing medians of each row.
    """
    indices = np.random.permutation(csr_mat.shape[1])
    medians = np.zeros(csr_mat.shape[0])
    mb_size = mb_size = int(2 ** np.log10(csr_mat.shape[1]) * 128)# Adjust this according to your preference
    n_batch = (csr_mat.shape[1] + mb_size - 1) // mb_size
    for i in range(n_batch):
        print("n_batch: ", n_batch)
        start_idx = i * mb_size
        end_idx = min((i + 1) * mb_size, csr_mat.shape[1])
        mb_indices = indices[start_idx:end_idx]
        row_data = np.asarray(csr_mat[:,mb_indices].todence())
        medians[mb_indices] = np.nanmedian(row_data, axis=0) + np.finfo(np.float64).eps 
    return medians

def RegNMF(rna_data, atac_data, batch_type, Meta_data=None, K=100, feature_cutperc=0.01, key_added="scReg_reduction", maxiter=40, copy='rna'):
    """
    Perform Coupled Non-negative Matrix Factorization (NMF) on RNA and ATAC data.

    Parameters:
    -----------
    rna_data : AnnData
        Annotated data for RNA.
    atac_data : AnnData
        Annotated data for ATAC.
    batch_type : str
        Type of batch information in Meta_data.
    Meta_data : pd.DataFrame, optional
        Metadata containing batch information. Default is None.
    K : int, optional
        Number of components for NMF. Default is 100.
    feature_cutperc : float, optional
        Feature cutoff precision. Default is 0.01.
    key_added : str, optional
        Key to store results in rna_data.obsm. Default is "scReg_reduction".
    maxiter : int, optional
        Maximum number of iterations for NMF. Default is 40.
    copy : str, optional
        Specify 'rna' or 'atac' to choose which data to copy results into. Default is 'rna'.

    Returns:
    --------
    AnnData
        Annotated data with NMF results added.
    """
    if isinstance(rna_data,AnnData) and isinstance(atac_data,AnnData):
        if Meta_data is None:
            scReg_reduction = RegNMF_Matrix(E=rna_data.X, O=atac_data.X, Meta_data=rna_data.obs, batch_type=batch_type, K=K, feature_cutperc=feature_cutperc, maxiter=maxiter)
        else:
            scReg_reduction = RegNMF_Matrix(E=rna_data.X, O=atac_data.X, Meta_data=Meta_data, batch_type=batch_type, K=K, feature_cutperc=feature_cutperc, maxiter=maxiter)

        if copy=='rna':
            rna_data.obsm[key_added] = scReg_reduction['H']
    
        else:
            rna_data.obsm[key_added] = scReg_reduction['H']
        return rna_data


def RegNMF_Matrix(E, O, Meta_data, batch_type, K=100, feature_cutperc=0.01, maxiter=40):


    """
    Perform Coupled Non-negative Matrix Factorization (NMF) on input matrices.

    Parameters:
    -----------
    E : csr_matrix
        Sparse matrix for RNA data.
    O : csr_matrix
        Sparse matrix for ATAC data.
    Meta_data : pd.DataFrame
        Metadata containing batch information.
    batch_type : str
        Type of batch information in Meta_data.
    K : int, optional
        Number of components for NMF. Default is 100.
    feature_cutperc : float, optional
        Feature cutoff precision. Default is 0.01.
    maxiter : int, optional
        Maximum number of iterations for NMF. Default is 40.

    Returns:
    --------
    dict
        Dictionary containing NMF results.
    """

    E = csr_matrix(E)
    O = csr_matrix(O)
    print("E_matrix: ", type(E))
    print("O_matrix: ", type(O))

    numCut = feature_cutperc * E.shape[0]
    expressed_cellN = np.array((E != 0).sum(axis=0)).squeeze()
    open_cellN = np.array((O != 0).sum(axis=0)).squeeze()

    print("E shape: ", E.shape)
    print("numCut : ", numCut)
    print("geneN: ", expressed_cellN)
    print("sum: ",(expressed_cellN > numCut).sum() )

    expressed_indices = np.where(expressed_cellN > numCut)[0]
    open_indices = np.where(open_cellN > numCut)[0]
    E = E[:, expressed_indices]
    O = O[:, open_indices]

    print("E shape after filtering: ", E.shape)
    print("O shape after filtering: ", O.shape)


    W2, H2 = core.perform_nmf(E.T, K)
    W2 = W2 / np.sqrt((H2 * H2).sum(axis=1))

    W1, H1 = core.perform_nmf(O.T, K)
    W1 = W1 / np.sqrt((H1 * H1).sum(axis=1))


    lambda1, lambda2 = core.defaultpar_CoupledNMF_default(PeakO=O,
                                                    W1 = W1, X=E,
                                                    W2 = W2, beta=1, arfa=1, withoutRE=True)



    batch_list = Meta_data[batch_type].unique()
    W20 = np.random.rand(E.shape[1], K)
    W10 = np.random.rand(O.shape[1], K)
    H0 = np.random.rand(E.shape[0], K).T

    for i in batch_list:
        W20 = np.append(W20, E[Meta_data[batch_type] == i,:].mean(axis = 0).reshape(-1,1),axis = 1)
        W10 = np.append(W10, O[Meta_data[batch_type] == i,:].mean(axis = 0).reshape(-1,1),axis = 1)
        H0 = np.append(H0,(Meta_data[batch_type] == i).values.astype('int').reshape(1,-1), axis=0)


    ans = core.nmf_cluster_joint_cross_domain_torch_sparse_mb(PeakO=O.T, X=E.T, lambda1=lambda1, W10=W10, W20=W20, H0=H0, K=100, maxiter=maxiter)
    return ans


