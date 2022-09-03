# the necessary tools collections. some common functions is implemented here.
from umap import UMAP
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import pandas as pd
import seaborn as sns

def _check_ndarray(emb):
    if(isinstance(emb,pd.DataFrame)):
        return emb.values
    elif(isinstance(emb,np.ndarray)):
        return emb
    else:
        raise RuntimeError(f"unknown type {type(emb)}")

def visualization_lmd(emb, lmd, filename, s=1):
    emb = _check_ndarray(emb)
    index = pd.isna(lmd)
    plt.figure()
    plt.scatter(emb[index, 0], emb[index, 1], s=s, c="red")
    pts = plt.scatter(emb[~index, 0], emb[~index, 1], s=s, c=lmd[~index], cmap="Blues_r")
    plt.colorbar(pts)

    plt.xlabel("UMAP_1")
    plt.ylabel("UMAP_2")
    plt.title(f"positive_rate: {lmd.notna().sum() / lmd.shape[0]:.4}")
    plt.tight_layout()

    if (filename is not None):
        plt.savefig(filename)
        plt.close()
    else:
        plt.show()


def visualization(emb,batch_col="batch",bio_col="cell_type"):
    ndarray_emb = emb[["UMAP_1", "UMAP_2"]].values
    plot(ndarray_emb, groupby=emb["batch"].values, filename=f"batch.png",
         dpi=500)
    plot(ndarray_emb, groupby=emb["cell_type"].values,
         filename=f"bio.png", dpi=500)

def get_umap(df,meta,batch_col="batch",bio_col="cell_type"):
    _, emb, pc = umap(df, downsampling=None)
    emb[batch_col] = pc[batch_col] = meta[batch_col]
    emb[bio_col] = pc[bio_col] = meta[bio_col]
    return emb


def umap(df,downsampling=None,n_pc=30,direct_select=False):
    """
    umap embedding that maps a high-dimensional data into 2-dimensional space.
    For the data that does not distinguish the information abundance in each dimension, it first maps the data to
    a 30-dimensional space by PCA. This is an common strategy in practice. And then the output of PCA is embeded into 2-dimensional
    space by UMAP.
    For the data, which information abundance is not different (for example, RPCI), we directly selects the first 30 principal components
    to feed to the UMAP method.
    """
    if(downsampling is not None):
        if(isinstance(downsampling,float)):
            downsampling=int(df.shape[0]*downsampling)
        #index=np.random.randint(0,df.shape[0],downsampling)
        index=np.random.choice(np.arange(df.shape[0]), downsampling, replace=False)


    row_index=None
    if(isinstance(df,pd.DataFrame)):
        row_index = df.index
        df=df.values

    if(downsampling is not None):
        df=df[index]# select n cells from ndarray
        row_index=row_index[index]
    else:
        index=np.arange(df.shape[0])

    emb_pc=None
    if(n_pc is not None):
        if(n_pc<df.shape[1]):
            if(direct_select):
                df=df[:,:n_pc]
            else:
                pca = PCA(n_components=n_pc)
                df = pca.fit_transform(df)
        emb_pc=df.copy()

    df = UMAP(
        n_neighbors=30,
        min_dist=0.3,
        metric='cosine',
        n_components=2,
        learning_rate=1.0,
        spread=1.0,
        set_op_mix_ratio=1.0,
        local_connectivity=1,
        repulsion_strength=1,
        negative_sample_rate=5,
        angular_rp_forest=False,
        verbose=False
    ).fit_transform(df)
    #df = UMAP(n_components=2).fit_transform(df)

    return index,\
           pd.DataFrame(df,index=row_index,columns=["UMAP_1","UMAP_2"]),\
           pd.DataFrame(emb_pc,index=row_index,columns=[f"PC_{i}" for i in range(emb_pc.shape[1])])



def scatterplot(df, Label1, Label2=None, fig_path=None):
    def fun(df, hue, s=13):
        # sns.set()
        sns.scatterplot(x='UMAP_1',
                        y='UMAP_2',
                        data=df,
                        hue=hue,
                        s=s)

    flag_label2=Label2 is not None

    if (fig_path is not None):
        plt.figure(figsize=(16, 8))
        if (flag_label2):
            plt.subplot(1, 2, 1)
        fun(df, hue=Label1)
        if (flag_label2):
            plt.subplot(1, 2, 2)
            fun(df, hue=Label2)

        plt.savefig(fig_path+".jpg")
        plt.close()


def plot(emb,groupby,filename,s=1,figsize=(),dpi=500):
    plt.figure()
    names = []
    for name, group in pd.DataFrame(emb).groupby(groupby):
        plt.scatter(group.values[:, 0], group.values[:, 1], s=s)
        names.append(name)
    plt.legend(names, bbox_to_anchor=(1.05, 0), loc=3, borderaxespad=0)

    plt.tight_layout()
    if (filename is not None):
        plt.savefig(filename,dpi=dpi)
        plt.close()
    else:
        plt.show()
