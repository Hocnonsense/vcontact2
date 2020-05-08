""" Functions to build the associations between the objects from the matrices. """
import numpy as np
import pandas as pd


def cluster_taxonomy(clusters, taxonomy, level, P, R):
    """Add the associations to the cluster and taxonomy dataframes.
    Using the maximum of precision taxonomic class of each cluster
    and the maximum of recall cluster of each taxonomic class.
    If this maximum is null, the position is set to NaN.

    Args:
        clusters (dataframe): to modify
        taxonomy (dataframe): to modify
        level (str): taxonomic level of P & R
        P (sparse_matrix): Precision matrix (cluster X taxonomic class)
        R (sparse_matrix): Recall matrix  (cluster X taxonomic class)

    Returns:
        (dataframe,dataframe): dataframe clusters with added columns "pos_level"
            and "precision_level" and dataframe taxonomy with added columns
            "pos_cluster" and "recall".
    """

    df = dict()
    df["pos"] = range(P.shape[0])
    df["pos_"+level] = np.squeeze(np.asarray(np.argmax(P, 1)))
    df["precision_"+level] = np.squeeze(np.asarray(np.max(P, 1)))
    df = pd.DataFrame(df)

    # If the max of precision is null, do not associate.
    df.loc[df["precision_"+level] == 0, "pos_"+level] = np.nan

    clusters = pd.merge(clusters, df, how='left')

    df = dict()
    df["pos"] = range(R.shape[1])
    df["recall"] = np.squeeze(np.asarray(np.max(R, 0)))
    df["pos_cluster"] = np.squeeze(np.asarray(np.argmax(R, 0)))
    df = pd.DataFrame(df)

    # If the max of recall is null, do not associate.
    df.loc[df["recall"] == 0, "pos_cluster"] = np.nan

    taxonomy = pd.merge(taxonomy, df, how='left')

    return clusters, taxonomy


def contig_cluster(contigs, B):
    """
    Associate each contig with its maximal-membership cluster.
    If the maximal-membership is null, the mbship cluster position
    is set to NaN.

    Args:
        contigs (dataframe): the contigs, with:
            pos: the position of the contig in the matrix.
        B (sparse_matrix): Membership matrix

    Returns:
        dataframe: contigs with pos_cluster_mbship column added.
    """

    cm = {}
    cm["pos"] = range(B.shape[0])
    cm["membership"] = np.squeeze(np.asarray(np.max(B, 1)))
    cm["pos_cluster_mbship"] = np.squeeze(np.asarray(np.argmax(B, 1)))
    cm = pd.DataFrame(cm)

    # If the max membership is null, do not associate.
    cm.loc[cm["membership"] == 0, "pos_cluster_mbship"] = np.nan

    contigs = pd.merge(contigs, cm)
    return contigs


def contig_taxonomy(contigs, taxonomy, clusters, level,
                    cluster_choice="pos_cluster_mbship"):
    """
    Associate each contig with its predicted taxonomic class

    Args:
        contigs (dataframe): the contigs, with a "cluster_choice" column.
        taxonomy (dataframe): names of the taxonomic classes at this level.
        level (str): taxonomic level
        clustering_choice (str): column in which the cluster of a contig
            can be found.

    Returns:
        dataframe: contigs with "predicted_level" column added.
    """
    cname = "predicted_{}".format(level)
    taxonomy = taxonomy.set_index("pos")
    clusters = clusters.set_index("pos")
    for n, clust_pos in enumerate(contigs[cluster_choice]):
        if not pd.isnull(clust_pos):
            tax_pos = clusters.loc[clust_pos, "pos_{}".format(level)]
            if not np.isnan(tax_pos):
                contigs.loc[n, cname] = taxonomy.name[tax_pos]
    contigs.fillna({cname: "Non affiliated"}, inplace=True)
    return contigs
