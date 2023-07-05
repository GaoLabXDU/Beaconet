Beaconet: a reference-free method for integration of multiple batches of scRNA-seq data in original molecular space
===========================
This project is a python implementation of Beaconet, which is introduced in the paper "A reference-free method for integration of multiple batches of scRNA-seq data in original molecular space".

Batch effect removal is crucial for single-cell data analysis, as the undesired batch-specific systematic variations among different technologies or experimental laboratories may distort biological signals and hinder the integration of datasets. As a reference-free batch effect removal method, Beaconet is particularly beneficial for researchers who lack prior knowledge about the selection of reference when integrate multiple batches of single-cell transcriptomic datasets. Furthermore, the preserved gene features in the corrected data by Beaconet can not only effectively characterize cell types and subpopulations, but also be directly used for downstream differential expression analysis.

All the preprocessed data, code scripts,and necessary outputs in the experiments of our paper are free available in figshare doi:10.6084/m9.figshare.20764843 to ensure reproducibility.

**Keywords:** single-cell data, batch effect, reference-free, molecular feature space;

## Author and Mantainer
This software is developed and mantained by Han Xu <hxu10670@gmail.com>

## Tutorial
This is a tutorial for beginer. We will introduce the whole batch effect removal processing step-by-step using a case task, the integration of two batches of human blood dendritic cell datasets. The batch effect of these two datasets is relative small, since they were generaetd by the same laboratory (Villani et al.) using the same sequencing technology (Smart-Seq2) in different plates. CD1C DC in batch 1 and CD141 DC in batch 2 were removed refering to Tran et al.

### Catalogue
* [Install](#Install)
* [Import the nessary package and function](#import-the-nessary-package-and-function)
* [Data preprocessing](#Data-preprocessing)
* [Correct batch effect for integration](#Correct-batch-effect-for-integration)
* [Visualization-and-evaluation](#Visualization-and-evaluation)
---------

# Install

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

## Import nessary package and function
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
## Correct batch effect for integration
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

![cell_type](https://github.com/xuxiaohan/Beaconet/blob/main/bio.png)
![batch](https://github.com/xuxiaohan/Beaconet/blob/main/batch.png)
![PMD](https://github.com/xuxiaohan/Beaconet/blob/main/local_merge_divergence.png)

## Experimental environment in our study

The requirements of Beaconet is included in the 'setup.py'. these packages will be automatically installed when install Beaconet.
For convenience of the audience to reproduce the results in our paper, we provied the specific version of the required environment and package in the experiments of our paper here.

Windows 10, HP Pavilion Gaming Desktop 690-07xx

Intel(R) Core(TM) i7-9700F CPU @ 3.00GHz 3.00GHz

16 GB memory. 64bit operating system

CUDA version:

nvcc: NVIDIA (R) Cuda compiler driver

Copyright (c) 2005-2019 NVIDIA Corporation

Built on Fri_Feb__8_19:08:26_Pacific_Standard_Time_2019

Cuda compilation tools, release 10.1, V10.1.105

Python 3.6.8

numpy==1.19.5

pandas==1.1.5

seaborn==0.9.0

matplotlib==3.3.4

tqdm==4.62.2

scikit-learn==0.24.1

torch==1.8.1+cu101

umap-learn==0.5.1

scipy==1.5.4

