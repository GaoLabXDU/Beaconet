from Beaconet import correction,visualization,get_umap,get_pmd,visualization_pmd,get_cluster
import pandas as pd
from sklearn.metrics import normalized_mutual_info_score as nmi
from sklearn.metrics import adjusted_rand_score as ari

if __name__ == '__main__':
    ##########################################
    # read scRNA-seq data
    dfs=[
        pd.read_csv("data/DC_batch0.csv",index_col=0),
        pd.read_csv("data/DC_batch1.csv",index_col=0)
    ]

    # read meta data for evaluation.
    # it is notable that the meta data is not used during correction.
    meta = pd.read_csv("data/DC_cell_meta.csv", index_col=0)
    meta = meta.reindex(list(dfs[0].index)+list(dfs[1].index))

    ###########################################
    # correction
    # the dfs is a list of DataFrame. the cells come from the same batch is organized in the same DataFrame.
    # the correction function returns the corrected data.
    result=correction(dfs)

    ############################################
    #save result
    result.to_csv("test/corrected.csv")
    #result=pd.read_csv("test.csv",index_col=0)

    # dimension reducetion
    ump=get_umap(result)
    ump["batch"] = meta["batch"]
    ump["cell_type"] = meta["cell_type"]

    #calculate the positive rate and the local merge divergence of positive cells
    positive_rate,pmd=get_pmd(ump,batch_col="batch",bio_col="cell_type")

    ###########################################
    # visualization
    visualization(ump,batch_col="batch",bio_col="cell_type",filename1="test/batch.png",filename2="test/bio.png")
    visualization_pmd(ump, pmd, filename="test/positive_merge_divergence.png")

    # clustering
    cluster=get_cluster(ump[["UMAP_1","UMAP_2"]],k=4)
    #calculate NMI
    print(f"NMI: {nmi(cluster,meta['cell_type'])}")
    #calculate ARI
    print(f"ARI: {ari(cluster,meta['cell_type'])}")
