# the script is to calculate the LMD metric for each method
#
import os
import numpy as np
import pandas as pd
from tqdm import trange
from sklearn.model_selection import ParameterGrid
from sklearn.neighbors import kneighbors_graph
from scipy.stats import entropy

def get_lmd(df,batch_col="batch",bio_col="cell_type"):
    LMD = compute_LME(X=df[["UMAP_1", "UMAP_2"]], meta=df[[batch_col, bio_col]], batch_label=batch_col,
                      bio_label=bio_col)
    positive_rate=LMD.notna().sum()/LMD.shape[0]
    return positive_rate,LMD

def concat(x,sep="_"):
    return sep.join([str(e) for e in x])


def compute_LME(X,meta,batch_label,bio_label,k=None):
    assert isinstance(X,pd.DataFrame)

    if(k is None):
        k=max(15,int(X.shape[0]*0.01))

    batch=meta[batch_label]
    bio=meta[bio_label]

    global_dist = pd.crosstab(batch, bio, normalize="columns")
    global_dist =global_dist.sort_index()

    knn = kneighbors_graph(X, mode='distance', n_neighbors=k,
                           metric="minkowski", p=2, include_self=False)

    bio=bio.values

    merge_entropy = []
    for i in trange(X.shape[0]):
        #if (knnpredict(knn, celltype, i, mode="distance")):
        if (knnpredict(knn, bio, i, mode="count")):
            bioinfo = bio[i]
            local = batch.iloc[knn[i].nonzero()[1]]
            local=local.value_counts(normalize=True).sort_index()
            local1=pd.Series(0.0,index=global_dist.index)
            local1[local.index]=local
            LME=divergence(local1, global_dist[bioinfo], cal_type="JS")
            merge_entropy.append(LME)
        else:
            merge_entropy.append(None)

    return pd.Series(merge_entropy,index=X.index,name="LME")


def knnpredict(knn,index,i,mode="count"):
    """
    weighted
    half

    """
    center=index[i]
    _,idx=knn[i].nonzero()
    index=index[idx]
    if(mode=="count"):
        df=pd.Series(index).value_counts(normalize=True)
        if(center in df.index and df[center]>0.5):
            return True
        else:
            return False
    else:
        raise RuntimeError("mode must be in {'count'}")


def divergence(p,q,cal_type="JS",eps=1e-10):
    p[p<eps]=eps
    q[q<eps]=eps
    if(cal_type == "qp"):
        res = entropy(q,p)#(q * np.log(q) - q * np.log(p)).sum()
    elif(cal_type == "pq"):
        res = entropy(p,q)#(p * np.log(p) - p * np.log(q)).sum()
    elif(cal_type=="JS"):
        mid=(p+q)/2
        res1=entropy(p,mid)#(p * np.log(p) - p * np.log(mid)).sum()
        res2 =entropy(q,mid) #(q * np.log(q) - q * np.log(mid)).sum()
        res=(res1+res2)/2
    else:
        raise RuntimeError(f"illegal value for cal_type, must be in {set(['pq','qp','JS'])}")
    return res


