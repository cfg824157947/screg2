# BETA
This package still in beta, please feel free to e-mail me for any question about the package and usage.

email: fchang[at]clemson.edu

# Install
``` bash
pip install screg2
```

# Usage

## Load data form h5ad file

If you have 2 anndata object that already finished scaling,
normalizing,\... You can just simply use `screg2.RegNMF()` to integrate
rna and atac data, with batch remove. `screg2.RegNMF()` will return the
rna anndata object with `scReg_reduction` in it `rna_adata.obsm`

``` python
import pandas as pd
import numpy as np
import scanpy as sc
import screg2 

rna_file = '/path/to/your/rna/h5ad/file'
atac_file = '/path/to/your/atac/h5ad/file'
rna_adata = sc.read_h5ad(rna_file)
atac_adata= sc.read_h5ad(atac_file)
rna_adata = screg2.RegNMF(rna_data=rna_adata, atac_data=atac_adata, Meta_data=rna_adata.obs,batch_type='sample', maxiter=100, key_added="scReg_reduction")
```

`screg2.RegNMF()` will return the rna anndata object with
`scReg_reduction` in it `rna_adata.obsm`

``` python
rna_adata.obsm['scReg_reduction']
```

``` python
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
```

::: warning
::: title
Warning
:::

If you do not have the meta data of cells for their sample infomation
the nex code might help you
:::

``` python
rna_adata.obs['sample'] = adata.obs_names.str.split("-").str[1]
```

## Clustering and Umap

``` python
sc.pp.neighbors(rna_adata, use_rep='scReg_reduction', key_added='scReg_neighbors', n_neighbors=40)
sc.tl.leiden(rna_adata, neighbors_key='scReg_neighbors', key_added='scReg_leiden')
sc.tl.umap(rna_adata, neighbors_key='scReg_neighbors')
sc.pl.umap(rna_adata, color="scReg_leiden", legend_loc="on data",legend_fontsize="small",size=6, title="scReg")
# If you want to save the figure you can
# sc.pl.umap(rna_adata, color="scReg_leiden", legend_loc="on data",legend_fontsize="small",size=6, save="_scReg.pdf" title="scReg")
```
