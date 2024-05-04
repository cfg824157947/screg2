import numpy as np
import time
import pandas as pd
import scanpy as sc
from . import core
from scipy.sparse import csc_matrix, csr_matrix, find
import psutil
import os 
from anndata import AnnData
from scipy import sparse

def sparse_row_median(X):
    """
    Compute the median of each row in a sparse CSR matrix.

    Parameters:
    - csr_mat (csr_matrix): Input sparse CSR matrix.

    Returns:
    - numpy.ndarray: Array containing medians of each row.
    """
    
    X = X.tocoo()
    df = pd.DataFrame({'Value': X.data, 'Group': X.row})
    medians = df.groupby('Group')['Value'].median()
    return medians.values

def normalize_scaling_default(X):
    b = sparse_row_median(X) + np.finfo(np.float64).eps  # Compute median along each column, ignoring NaN values
    scaled_X = X.T.multiply(1/b * np.mean(b)).T  # Normalize and scale by mean
    return (scaled_X)  # Convert back to sparse matrix format if needed

def tfidf(X):
    """
    Compute tfidf for matrix X (cell X gene)
    https://github.com/stuart-lab/signac/blob/HEAD/R/preprocessing.R
    """
    tf = X / np.log(1 + X.sum(axis = 1))
    idf = np.log(1 + (X.shape[0]/(1+X.sum(axis = 0))))
    #idf[idf == np.inf] = 0
    return tf.multiply(idf.A1)
    

def RegNMF_h5(h5_file, barcodes=None):
    """
    Perform Coupled Non-negative Matrix Factorization (NMF) on RNA and ATAC data from h5 file. 

    Parameters:
    -----------
    h5_file : str
        Single cell multiome h5 file
    barcodes : str
        The barcodes of cells that you want use
    Returns:
    --------
    AnnData
        Annotated data with NMF results added and normalized layer.
    """
    adata = sc.read_10x_h5(h5_file, gex_only=False)
    if type(barcodes) != type(None):
        adata = adata[adata.obs_names.isin(barcodes),:]
    sample = adata.obs_names.str.split('-').str[1]
    adata.obs = pd.DataFrame(sample, index=adata.obs_names, columns=["sample"])
    rna_adata = adata[:,adata.var['feature_types']=='Gene Expression']
    atac_adata = adata[:,adata.var['fecture_types']=='Peaks']
    rna_adata, atac_adata = RegNMF(rna_data=rna_adata,atac_data=atac_adata, batch_type=sample)
    adata.layers['norm'] = adata.X
    adata.layers['norm'][:,adata.var['feature_types']=='Peaks'].data = atac_adata.X.tocsr().data
    adata.layers['norm'][:,adata.var['feature_types']=='Gene Expression'].data = rna_adata.X.tocsr().data
    adata.obsm['scReg_reduction'] = rna_adata.obsm['scReg_reduction']
    return adata
    

def RegNMF(rna_data, atac_data, batch_type, Meta_data=None, K=100, feature_cutperc=0.01, key_added="scReg_reduction", maxiter=100, copy='rna', TFIDF=False, normalize=False):
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
        if TFIDF:
            atac_data.layers['norm'] = tfidf(atac_data.X).tocsr()
            atac_data.X = atac_data.layers['norm']
        if normalize:
            rna_data.layers['norm'] = normalize_scaling_default(rna_data.X).tocsr()
            rna_data.X = atac_data.layers['norm']
            atac_data.layers['norm'] = normalize_scaling_default(atac_data.X).tocsr()
            atac_data.X = atac_data.layers['norm']

        if Meta_data is None:
            scReg_reduction = RegNMF_Matrix(E=rna_data.X, O=atac_data.X, Meta_data=rna_data.obs, batch_type=batch_type, K=K, feature_cutperc=feature_cutperc, maxiter=maxiter, TFIDF=TFIDF, normalize=normalize)
        else:
            scReg_reduction = RegNMF_Matrix(E=rna_data.X, O=atac_data.X, Meta_data=Meta_data, batch_type=batch_type, K=K, feature_cutperc=feature_cutperc, maxiter=maxiter)

        if copy=='rna':
            rna_data.obsm[key_added] = scReg_reduction['H'].T
    
        else:
            rna_data.obsm[key_added] = scReg_reduction['H'].T
        return rna_data, atac_data


def RegNMF_Matrix(E, O, Meta_data, batch_type, K=100, feature_cutperc=0.01, maxiter=40,TFIDF=True,nomalize=True):


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
    E = csr_matrix(E).astype(float)
    O = csr_matrix(O).astype(float)
    print("E_matrix: ", type(E))
    print("O_matrix: ", type(O))



    numCut = feature_cutperc * E.shape[0]
    expressed_cellN = np.array((E != 0).sum(axis=0)).squeeze()
    open_cellN = np.array((O != 0).sum(axis=0)).squeeze()

    print("E shape: ", E.shape)
    print("numCut : ", numCut)
    print("geneN: ", expressed_cellN)
    print("sum: ",(expressed_cellN > numCut).sum())

    expressed_indices = np.where(expressed_cellN > numCut)[0]
    open_indices = np.where(open_cellN > numCut)[0]
    E = E[:, expressed_indices]
    O = O[:, open_indices]

    print("E shape after filtering: ", E.shape)
    print("O shape after filtering: ", O.shape)


    #E = E.tocsc()
    #O = O.tocsc()
    #print("E_matrix: ", type(E))
    #print("O_matrix: ", type(O))

    start_t = time.time()
    W2, H2 = core.perform_nmf(E.T, K)
    W2 = W2 / np.sqrt((H2 * H2).sum(axis=1))

    W1, H1 = core.perform_nmf(O.T, K)
    W1 = W1 / np.sqrt((H1 * H1).sum(axis=1))

    end_t = time.time()
    exetime = end_t- start_t
    
    print(f"NMF:  {exetime} sec")
    start_t = time.time()
    lambda1, lambda2 = core.defaultpar_CoupledNMF_default(PeakO=O,
                                                    W1 = W1, X=E,
                                                    W2 = W2, beta=1, arfa=1, withoutRE=True)


    print("lambda1: ", lambda1)
    end_t = time.time()
    exetime = end_t- start_t
    
    print(f"defaultpar_CoupledNMF_default:  {exetime} sec")
    batch_list = Meta_data[batch_type].unique()
    W20 = np.random.rand(E.shape[1], K)
    W10 = np.random.rand(O.shape[1], K)
    H0 = np.random.rand(E.shape[0], K).T

    for i in batch_list:
        W20 = np.append(W20, E[Meta_data[batch_type] == i,:].mean(axis = 0).reshape(-1,1),axis = 1)
        W10 = np.append(W10, O[Meta_data[batch_type] == i,:].mean(axis = 0).reshape(-1,1),axis = 1)
        H0 = np.append(H0,(Meta_data[batch_type] == i).values.astype('int').reshape(1,-1), axis=0)


    start_t = time.time()
    ans = core.nmf_cluster_joint_cross_domain_torch_sparse_mb(PeakO=O.T, X=E.T, lambda1=lambda1, W10=W10, W20=W20, H0=H0, K=K, maxiter=maxiter)
    end_t = time.time()
    exetime = end_t- start_t
    
    print(f"CoupledNMF:  {exetime} sec")
    ans['H'] = ans['H'][:K,:]
    return ans


