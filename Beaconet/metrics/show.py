# the necessary tools collections. some common functions is implemented here.
from umap import UMAP
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans

"""
This is a collection of function tools for visualizing the results. It mainly contains four useful functions, including
*visualization*, *get_umap*, *visualization_pmd*, *get_cluster*.

*visualization* is for plotting scatter figures, the dot will be colored by the given labels (cell type or batch name).
The cells have the same label will be plotted using the same color.

*get_umap* is a function for transforming the high-dimensional feature representation of cells to two-dimensional
UMAP space, which is beneficial for visualization and clustering.

*visualization_pmd* is for visualizing the Positive Merge Divergence for cells, in which the dots with blue color are
the positive cells, and the red color for negative cells. For the positive cells, the deeper blue color indicate more
better mixture. 

*get_cluster* is for grouping the cells into clusters using k-means.
"""

def visualization_pmd(emb, pmd, filename, s=1):
    """
        This is the API to visualizing the PMD metric for each cell in scatter plot. the dots with blue color are
        the positive cells, and the red color for negative cells. For the positive cells, the deeper blue color
        indicate more better mixture.

        Parameters
        ----------
        emb  : numpy.ndarray or pandas.DataFrame. The ndarray or DataFrame is expected to contain the two-dimensional
            UMAP features of cells, the shape of the matrix is expected as (n_cells, 2).

        pmd: pd.Series, each element for the merge divergence of a positive cell, and None for negative cells. Its shape
            should be (n_cells, 1)

        filename: str or None.
            save the drawn figure as file {filename} if str is passed.
            show the drawn figure if NOne is passed.

        s: float, default 1. The dot size in the figure.

        =========
        The usage demo are included in the Beaconet/test/demo.py
        =========
    """

    emb = _check_ndarray(emb)
    index = pd.isna(pmd)
    plt.figure()
    plt.scatter(emb[index, 0], emb[index, 1], s=s, c="red")
    pts = plt.scatter(emb[~index, 0], emb[~index, 1], s=s, c=pmd[~index], cmap="Blues_r")
    plt.colorbar(pts)

    plt.xlabel("UMAP_1")
    plt.ylabel("UMAP_2")
    plt.title(f"positive_rate: {pmd.notna().sum() / pmd.shape[0]:.4}")
    plt.tight_layout()

    if (filename is not None):
        plt.savefig(filename)
        plt.close()
    else:
        plt.show()


def visualization(emb,batch_col="batch",bio_col="cell_type",filename1="batch.png",filename2="bio.png"):
    """
        This is the API to visualizing cells in scatter plot. the color of dots depends on the given labels (the batch
            label, cell type label). This function will generate two figures, one for batch and the other for cell type.

        Parameters
        ----------
        emb  : pandas.DataFrame. The DataFrame is expected to contain at least four columns, including the
            two-dimensional UMAP features, a column with batch id and a column with cell type of cells.

        batch_col: str, the column name of batch label in *emb*.

        bio_col: str, the column name of cell type in *emb*.

        filename1: str or None, default 'batch.png'. Save or show batch figure.
            save the drawn figure as file {filename1} if str is passed.
            show the drawn figure if NOne is passed.

        filename2: str or None, default 'bio.png'. Save or show cell type figure.
            save the drawn figure as file {filename2} if str is passed.
            show the drawn figure if NOne is passed.

        =========
        The usage demo are included in the Beaconet/test/demo.py
        =========
    """
    ndarray_emb = emb[["UMAP_1", "UMAP_2"]].values
    plot(ndarray_emb, groupby=emb[batch_col].values, filename=filename1,
         dpi=500)
    plot(ndarray_emb, groupby=emb[bio_col].values,
         filename=filename2, dpi=500)

def get_umap(df):
    """
        This is the API to obtain the two-dimensional UMAP embedding features of the high-dimensional represented cells.
        If the number of original features is larger than 30, PCA will be applied before UMAP embedding. If the number
        of original features is smaller or equal to 30, the UMAP embedding will be directly calculated.

        Return
        ----------
        emb: pandas.DataFrame. It is a matrix containing the UMAP embedding features of given cells.

        Parameters
        ----------
        df  : pandas.DataFrame or numpy.ndarray. A matrix contains the cells with high-dimensional features.
            The shape is (n_cells, n_features).

        =========
        The usage demo are included in the Beaconet/test/demo.py
        =========
    """
    _, emb, pc = umap(df, downsampling=None)
    return emb


def get_cluster(X, k):
    """
    This is an API for obtaining clusters from given matrix X (n_cells, n_features). The features of X will be reduced
    by *get_umap* if it is larger than 30.

    Return
    ----------
    labels : ndarray of shape (n_cells,)
        Index of the cluster each sample belongs to.

    Parameters
    ----------
    X: pandas.DataFrame or numpy.ndarray. A matrix contains the cells with features.
        The shape is (n_cells, n_features).

    k: int. The expected number of clusters, which will be passed to k-means algorithm.

    =========
    The usage demo are included in the Beaconet/test/demo.py
    =========
    """
    if (X.shape[1] > 30):
        X = get_umap(X)

    return KMeans(n_clusters=k, n_init=100).fit_predict(X)

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


def _check_ndarray(emb):
    if(isinstance(emb,pd.DataFrame)):
        return emb.values
    elif(isinstance(emb,np.ndarray)):
        return emb
    else:
        raise RuntimeError(f"unknown type {type(emb)}")

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
