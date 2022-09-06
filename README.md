Beaconet
===========================
This project is an python implementation of Beaconet, which is introduced in the paper "A reference-free method for integration of multiple batches of scRNA-seq data".
The raw scripts and data for reproduce all results in our paper is free available in figshare doi:10.6084/m9.figshare.20764843
	
|Author|E-mail|
|---|---|
|Xu Han|hxu10670@gmail.com|
|Gao Lin|-|

## Catalogue
1. Install
2. import the nessary package and function
3. Data preprocessing
4. Correction and integration
5. Visualization and evaluation
6. Environments
---------
## Install
Our algorithm, Beaconet, is based on deep learning. Beaconet is implemented on the PyTorch framework.
If there are GPUs on your machine, it is strongly support that to install the pytorch with cuda.
Here is some example to install pytorch.

install the pytorch with cuda 11.3 using pip:
```Bash
pip3 install torch --extra-index-url https://download.pytorch.org/whl/cu113
```
install the pytorch without cuda using pip:
```Bash
pip install torch
```
See the offical site of PyTorch for more information about the installation of PyTorch.
https://pytorch.org/get-started/locally/

Install Beaconet using pip.
```Bash
pip install Beaconet
```

## import the nessary package and function
```python
import torch as t
from Beaconet import correction,visualization,get_umap,get_lmd,visualization_lmd
import pandas as pd
from glob import glob
```

## Data preprocessing
The data preprocessing is important to the data-driven algorithm. In Beaconet, we assume the inputs is the log-scaled TPM/count, and filter the non-positive values in the correction module. The details of preprocessing for Beaconet is described below:
1. filter the low-quality genes. We removed the genes that expressed in less than 3 cells. If the expression value of gene a on cell b is zero, we consider the gene does not expressed in the cell.
2. filter the low-quality cells. We removed the cells that expressed less than 300 genes.
3. normalize to library size and exclude the highly expressed values.
4. extract the top 2000 highly variable genes.
We used the preprocessing functions in scanpy package to preprocess the raw data.
```python
    adata=sc.AnnData(X=df,obs=meta)
    sc.pp.filter_genes(data=adata,min_cells=3,inplace=True)
    sc.pp.filter_cells(data=adata,min_genes=300,inplace=True)
    sc.pp.normalize_total(adata,target_sum=10000,exclude_highly_expressed=True,inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat", batch_key="batch")
    adata=adata[:,adata.var['highly_variable']]
    
    for batch,data in df.groupby(adata.obs["batch_id"]):
    	data.to_csv(f"hvg2000_batch{batch}.csv")
```
## Correction and integration
In this section, we will show you how to apply Beaconet on an example data.
the example data is saved in the folder named 'data'.

1. read the preprocessed data and meta data. It is notable that Beaconet only use the batch label information. The cell type and other information is used for visualization and evaluation of the correction results, rathter than used during correction.
```python
    # read scRNA-seq data
    dfs=[
        pd.read_csv("data/DC_batch0.csv",index_col=0),
        pd.read_csv("data/DC_batch1.csv",index_col=0)
    ]
    
    # read meta data for evaluation.
    meta = pd.read_csv("data/DC_cell_meta.csv", index_col=0)
    meta = meta.reindex(list(dfs[0].index)+list(dfs[1].index))
```
2. correction and save the result. The dfs is a list of DataFrame. the cells come from the same batch is organized in the same DataFrame.The correction function returns the corrected data.
```python
    result=correction(dfs)
    result.to_csv("test.csv")
```
## Visualization and evaluation
1. dimension reduction
```python
    ump=get_umap(result,meta,batch_col="batch",bio_col="cell_type")
```
2. calculate the local merge divergence
```python
    positive_rate,lmd=get_lmd(ump,batch_col="batch",bio_col="cell_type")
```
3. plot
```python
    visualization(ump,batch_col="batch",bio_col="cell_type")
    visualization_lmd(ump, lmd, filename="local_merge_divergence.png")
```
