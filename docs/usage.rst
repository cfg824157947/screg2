Usage
=====

Basic Usage with Anndata
-------------------------
Load data as anndata (h5ad)

.. code-block:: python

  import pandas as pd
  import scanpy as sc
  import screg2 

  adata = sc.read_10x_mtx(path='/data2/duren_lab/palmetto/cham/toyData/filtered_feature_bc_matrix/',cache=True,gex_only=False)
  rna_adata = adata[:,toyData.var['feature_types']=='Gene Expression'] 
  atac_adata = adata[:,toyData.var['feature_types']=='Peaks']
  sc.pp.normalize_total(rna_adata) 
  sc.pp.normalize_total(atac_adata) 



  # This part will return redaction to adata_rna
  adata_rna = screg2.RegNMF(E = rna_adata, O =atac_adata, batch_type = "sample")









