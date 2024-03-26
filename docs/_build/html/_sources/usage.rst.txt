Usage
=====

Basic Usage with Anndata
-------------------------
Load data as anndata (h5ad)

.. code-block:: python

  import pandas as pd
  import scanpy as sc
  import project as scReg

  adata = sc.read_10x_mtx(path='/data2/duren_lab/palmetto/cham/toyData/filtered_feature_bc_matrix/',cache=True,gex_only=False)

  # This part will return redaction to adata_rna
  adata_rna = scReg.RegNMF(E = adata[:,toyData.var['feature_types']=='Gene Expression'], O = adata[:,toyData.var['feature_types']=='Peaks'], batch_type = "sample")









