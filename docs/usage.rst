Usage
=====

Load data form h5ad file
____________________________
If you have 2 anndata object that already finished scaling, normalizing,... You can just simply use ``screg2.RegNMF()`` to integrate rna and atac data, with batch remove. 
``screg2.RegNMF()`` will return the rna anndata object with ``scReg_reduction`` in it ``rna_adata.obsm`` 

.. code-block:: python

  import pandas as pd
  import numpy as np
  import scanpy as sc
  import screg2 

  rna_file = '/path/to/your/rna/h5ad/file'
  atac_file = '/path/to/your/atac/h5ad/file'
  rna_adata = sc.read_h5ad(rna_file)
  atac_adata= sc.read_h5ad(atac_file)
  rna_adata = screg2.RegNMF(rna_data=rna_adata, atac_data=atac_adata, Meta_data=rna_adata.obs,batch_type='sample', maxiter=100, key_added="scReg_reduction")


``screg2.RegNMF()`` will return the rna anndata object with ``scReg_reduction`` in it ``rna_adata.obsm``

.. code-block:: python
                
  rna_adata.obsm['scReg_reduction']
.. code-block:: python
                
  array([[4.13353097e-14, 1.08500913e-10, 3.27313772e-11, ...,
          0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
        [2.09097071e-14, 1.96949657e-16, 1.67404740e-11, ...,
          0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
        [1.12355639e-10, 3.41268614e-11, 6.89114276e-11, ...,
          0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
        ...,
        [6.38938638e-14, 7.88140180e-11, 5.65415185e-13, ...,
          0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
        [1.88311115e-10, 1.70259133e-10, 6.51093299e-12, ...,
          0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
        [2.15576528e-12, 9.40752216e-11, 1.57056569e-12, ...,
          0.00000000e+00, 0.00000000e+00, 0.00000000e+00]])

.. warning::
  If you do not have the meta data of cells for their sample infomation the nex code might help you

.. code-block:: python

  rna_adata.obs['sample'] = adata.obs_names.str.split("-").str[1]

Clustering and Umap
____________________________
.. code-block:: python

  sc.pp.neighbors(rna_adata, use_rep='scReg_reduction', key_added='scReg_neighbors', n_neighbors=40)
  sc.tl.leiden(rna_adata, neighbors_key='scReg_neighbors', key_added='scReg_leiden')
  sc.tl.umap(rna_adata, neighbors_key='scReg_neighbors')
  sc.pl.umap(rna_adata, color="scReg_leiden", legend_loc="on data",legend_fontsize="small",size=6, title="scReg")
  # If you want to save the figure you can
  # sc.pl.umap(rna_adata, color="scReg_leiden", legend_loc="on data",legend_fontsize="small",size=6, save="_scReg.pdf" title="scReg")




Load data from h5 file
-------------------------

.. code-block:: python

  import scanpy as sc
  h5_file = '/path/to/your/h5/file'
  adata = sc.read_10x_h5(h5_file,gex_only=False)
  # Make sample name Meta data. Makesure your cell barcode are end with -sample_name like this "AATTAATT-34"
  adata.obs['sample'] = adata.obs_names.str.split("-").str[1]

Create RNA and ATAC Anndata object 
-------------------------------------

  #start filtering cell by rna
  rna_adata= adata[:,adata.var['feature_types']=='Gene Expression']
  atac_adata = adata[:,toyData.var['feature_types']=='Peaks']
  atac_adata = 

.. code-block:: python

  import pandas as pd
  import numpy as np
  import scanpy as sc
  import screg2 
  file = "/path/to/your/h5/file"
  adata = sc.read_10x_mtx(path='/data2/duren_lab/palmetto/cham/toyData/filtered_feature_bc_matrix/',cache=True,gex_only=False)
  rna_adata = adata[:,toyData.var['feature_types']=='Gene Expression'] 
  atac_adata = adata[:,toyData.var['feature_types']=='Peaks']
  sc.pp.normalize_total(rna_adata) 
  sc.pp.normalize_total(atac_adata) 



  # This part will return redaction to adata_rna
  adata_rna = screg2.RegNMF(E = rna_adata, O =atac_adata, batch_type = "sample")









