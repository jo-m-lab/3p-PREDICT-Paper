#! /usr/local/bin/python3

import os
import numpy as np
import pandas as pd
import sklearn as sk
import pickle as pkl
import scanpy as sc

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap

from scipy.cluster import hierarchy as hclst

# Args
FGID13FILE = os.getenv("FGID13FILE") # 'FGID13_T1C3_SCT.loom'
CD14FILE = os.getenv("CD14FILE") # 'CD14_T1C3_SCT.loom'
FGID13MODEL = os.getenv("FGID13MODEL") # 'T1C3_FG13_SCT_model.pkl'
CD14MODEL = os.getenv("CD14MODEL") # 'T1C3_CD14_SCT_model.pkl'
FGID13META = os.getenv("FGID13META") # 'map_fgid_clusters/clmeta_FG13_all.tsv'
CD14META = os.getenv("CD14META") # 'map_crohns_clusters/clmeta_CD14_all.tsv'
CLASSCOLUMN = os.getenv("CLASSCOLUMN")
OUTFOLDER = os.getenv("OUTFOLDER")

def main():
    # read datasets
    fgid = sc.read_loom(FGID13FILE)
    cd14 = sc.read_loom(CD14FILE)
    # read models
    with open(FGID13MODEL, "rb") as fl:
        fgid_rfmod = pkl.load(fl)
    with open(CD14MODEL, "rb") as fl:
        cd14_rfmod = pkl.load(fl)
    # load meta
    fgid_meta = pd.read_csv(FGID13META, sep="\t")
    cd14_meta = pd.read_csv(CD14META, sep="\t")
    # get mtxs
    cd14_mtx = make_prediction_mtx(cd14, genenames=fgid.var.index)
    fgid_mtx = make_prediction_mtx(fgid, genenames=cd14.var.index)
    # get preds
    print("predicting...")
    cpreds = fgid_rfmod.predict(cd14_mtx)
    cprobs = fgid_rfmod.predict_proba(cd14_mtx)
    fpreds = cd14_rfmod.predict(fgid_mtx)
    fprobs = cd14_rfmod.predict_proba(fgid_mtx)

    print("storing results...")
    cd14_preds = pd.DataFrame({
        **{"CellID": cd14.obs.index,
        "manualannotation": cd14_meta\
                                .set_index("CellID")\
                                .loc[cd14.obs.index, CLASSCOLUMN],
        "cellsubsetspreds": cpreds},
        **{k: v for k, v in zip(fgid_rfmod.classes_, cprobs.T)}
    })



    print("storing results...")
    fgid_preds = pd.DataFrame({
        **{"CellID": fgid.obs.index,
        "manualannotation": fgid_meta\
                                .set_index("CellID")\
                                .loc[fgid.obs.index, CLASSCOLUMN],
        "cellsubsetspreds": fpreds},
        **{k: v for k, v in zip(cd14_rfmod.classes_, fprobs.T)}
    })

    cd14_preds.to_csv(
        f"{OUTFOLDER}/preds_CD14_on_FGID13.tsv", index=False, sep="\t")
    fgid_preds.to_csv(
        f"{OUTFOLDER}/preds_FGID13_on_CD14.tsv", index=False, sep="\t")

    # get probs as aggregated matrix
    probsFGID = fgid_preds\
        .set_index("CellID")\
        .loc[fgid_meta.set_index("CellID")\
             .loc[fgid_preds.CellID, "is_doublet"] == False]\
        .drop("cellsubsetspreds", axis=1)\
        .groupby("manualannotation")\
        .apply(np.mean, axis=0)
    probsCD14 = cd14_preds\
        .set_index("CellID")\
        .loc[cd14_meta.set_index("CellID")\
             .loc[cd14_preds.CellID, "is_doublet"] == False]\
        .drop("cellsubsetspreds", axis=1)\
        .groupby("manualannotation")\
        .apply(np.mean, axis=0)

    # %%
    ## plot confusion maps
    ##################################

    figsize = (20,20)
    sns.set(rc={'figure.figsize':figsize})
    p = sns.heatmap(probsCD14,
                xticklabels=probsCD14.columns,
                yticklabels=probsCD14.index);
    p.tick_params(axis="x", rotation=270)
    plt.tight_layout()
    plt.savefig(f"{OUTFOLDER}/confusion_heatmap_CD14.pdf", figsize=figsize)

    figsize = (20, 20)
    sns.set(rc={'figure.figsize':figsize})
    p = sns.heatmap(probsFGID,
                    xticklabels=probsFGID.columns,
                    yticklabels=probsFGID.index);
    p.tick_params(axis="x", rotation=270)
    plt.tight_layout()
    plt.savefig(f"{OUTFOLDER}/confusion_heatmap_FGID.pdf", figsize=figsize)


    # %%
    ## plot dendrograms and make linkage matrix of cellsubsets
    #############################################################

    # Has potential issues where there is only one group

    sns.set_style("ticks")
    fgid_links = hclst.linkage(probsFGID,
                               method="complete",
                               metric='correlation',
                               optimal_ordering=True)
    figsize = (5,10)
    fig, axs = plt.subplots(figsize=figsize);
    hclst.dendrogram(fgid_links,
                     show_leaf_counts=False,
                     leaf_label_func=lambda i: probsFGID.index[i], ax = axs);
    axs.set_xticklabels(labels=axs.get_xticklabels(), rotation=270, ha="right");
    plt.tight_layout()
    plt.savefig(
        f"{OUTFOLDER}/dendrogram_complete_correlationOfPredProbs_FGID13.pdf",
        figsize=figsize)


    cd14_links = hclst.linkage(probsCD14,
                               method="complete",
                               metric='correlation',
                               optimal_ordering=True)
    figsize = (10,10)
    fig, axs = plt.subplots(figsize=figsize);
    hclst.dendrogram(cd14_links,
                     show_leaf_counts=False,
                     leaf_label_func=lambda i: probsCD14.index[i], ax = axs);
    axs.set_xticklabels(labels=axs.get_xticklabels(), rotation=270, ha="right");
    plt.tight_layout()
    plt.savefig(
        f"{OUTFOLDER}/dendrogram_complete_correlationOfPredProbs_CD14.pdf",
        figsize=figsize)

    # Get ordered labels
    orderedFGIDlabels = list(probsFGID.index[hclst.leaves_list(fgid_links)])
    orderedCD14labels = list(probsCD14.index[hclst.leaves_list(cd14_links)])

    # order probs
    probsFGIDord = probsFGID.loc[orderedFGIDlabels, orderedCD14labels]
    probsCD14ord = probsCD14.loc[orderedCD14labels, orderedFGIDlabels]


    # %%
    ## Make melted frame for dotplot
    dot = (probsFGIDord - probsCD14ord.T)\
        .reset_index()\
        .melt(id_vars="manualannotation")
    dot["magnitude"] = (probsFGIDord + probsCD14ord.T)\
        .reset_index()\
        .melt(id_vars="manualannotation")\
        .value
    dot.columns = ["fgidLabel", "cd14Label", "bias", "magnitude"]
    # save full and filter based on magnitude
    dotall = dot.copy()
    th = dotall.magnitude.quantile(0.9)
    dot = dotall.loc[dotall.magnitude>=th, :]
    # get plot ready tick labels
    dot['FGID_Label'] = "FG." + dot.fgidLabel.str.ljust(max(dot.fgidLabel.str.len()), ".").str[3:]
    dot['Crohns_Label'] = "CD." + dot.loc[:, "cd14Label"].str[3:]

    #%%
    figsize = (10,10)
    sns.set(rc={'figure.figsize':figsize})
    sns.set_style("ticks")
    p = sns.jointplot(x=dotall.bias, y=dotall.magnitude)
    p.ax_joint.axhline(th, ls='--', color='red')
    p.ax_joint.text(0.1, 0, f"th: {th:2.3f}")
    plt.tight_layout()
    plt.savefig(
        f"{OUTFOLDER}/paired_confidence_distplot_th{th:2.3f}.pdf",
        figsize=figsize)


    # %%
    ## plot dotplot
    ##############################################
    mymap = ListedColormap(
        sns.diverging_palette(10, 255,
                              l=40,
                              s=99,
                              n=100,
                              center="light").as_hex())
    figsize = (20,12)
    fig = plt.figure(figsize=figsize)
    sns.set_style("white")
    sns.set_style({'font.family':'monospace'})
    p = plt.scatter(x=dot.Crohns_Label,
                    y=dot.FGID_Label,
                    c=dot.bias,
                    s=dot.magnitude*200,
                    edgecolors="black",
                    cmap=mymap,
                    norm=colors.TwoSlopeNorm(
                            vmin=min(min(dot.bias), -0.001),
                            vcenter=0,
                            vmax=max(max(dot.bias), 0.001)),
                    label="Confidence")
    plt.xticks(rotation=270, horizontalalignment='center', verticalalignment="top")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(f"{OUTFOLDER}/paired_confidence_dotplot.pdf", figsize=figsize)


# %%
def make_prediction_mtx(adata, cellnames=None, genenames=None, layername="norm_data"):
    """ From adata select layer with gene names, if gene name doesn't exist provides column of zeros"""
    from scipy import sparse
    if cellnames is None:
        cells = np.arange(adata.shape[0])
    else:
        cells = np.where(np.in1d(adata.obs.index, cellnames))[0]
        cells = cells[np.argsort(adata.obs.index[cells])]
    if genenames is None:
        mask = np.ones(adata.shape[1], dtype=bool)
        gene_map_idx = np.arange(adata.shape[1])
    else:
        mask = np.in1d(genenames, adata.var.index)
        gene_map_idx = np.array(
            [np.where(adata.var.index == g)[0][0] for g in genenames[mask]])
    # make sparse matrix with shape of selected cells and chosen genes
    mtx = sparse.lil_matrix(np.zeros((cells.shape[0], mask.shape[0])))
    # map loomfile layer mtx to specified order of cells and genes
    mtx[:, mask] = adata.layers[layername][cells, :][:, gene_map_idx]
    return sparse.csr_matrix(mtx)


if __name__ == '__main__':
    main()