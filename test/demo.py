import torch as t
from Beaconet import correction,visualization,get_umap,get_lmd,visualization_lmd
import pandas as pd
from glob import glob

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
    result.to_csv("test.csv")
    #result=pd.read_csv("test.csv",index_col=0)
    # dimension reducetion
    ump=get_umap(result,meta,batch_col="batch",bio_col="cell_type")
    #calculate the local merge divergence
    positive_rate,lmd=get_lmd(ump,batch_col="batch",bio_col="cell_type")

    ###########################################
    # visualization
    visualization(ump,batch_col="batch",bio_col="cell_type")
    visualization_lmd(ump, lmd, filename="local_merge_divergence.png")
