""" Functions to export final node tables that includes VC and individual contigs """
import logging
import os
import pandas as pd
import numpy as np

# import math
# from functools import reduce
from scipy.spatial import distance
from scipy.stats import mannwhitneyu
from scipy.cluster.hierarchy import linkage, cophenet

import vcontact2.cluster_refinements

# np.warnings.filterwarnings('ignore')

pd.set_option("display.max_rows", 100)
pd.set_option("display.max_columns", 10)
pd.set_option("display.width", 10000)

logger = logging.getLogger(__name__)


def find_excluded(merged, ntw, c1_df):
    """
    Temporary placement of singleton, overlap and outlier detection. These three components should be refactored
    "closer" to the locations where they'd be initially generated. Overlap genomes *are already identified* and saved
    in a dataframe, but that's only for overlapping genome removal, the dataframe/overlaps are never saved to a file.
    The outliers/singletons require a multi-step process, as ClusterONE does not state which genomes are outliers, so
    they must be back-calculated.

    singletons = merged - ntw
    outliers = ntw - c1
    overlap = c1 overlaps

    merged: db + user sequences, represents all genomes that went into the analysis
    ntw: network, represents genomes in network, otherwise they were excluded by threshold


    :return:
    """

    # Get list of all genomes
    merged_df = pd.read_csv(merged, header=0, index_col=0)
    merged_df.rename(
        columns={
            "kingdom": "Kingdom",
            "phylum": "Phylum",
            "class": "Class",
            "order": "Order",
            "family": "Family",
            "genus": "Genus",
            "contig_id": "Genome",
        },
        inplace=True,
    )
    contigs = set(merged_df["Genome"].tolist())

    merged_df["VC Status"] = np.nan

    # Get list of genomes that made it to the network, those that didn't did not pass the thresholds and => singletons
    network_df = pd.read_csv(
        ntw,
        header=None,
        names=["source", "target", "weight"],
        index_col=None,
        delimiter=" ",
    )
    nodes = set(network_df["source"].tolist() + network_df["target"].tolist())

    # Set up final destination for singletons
    # Pseudomonas~virus~Pf1 NOT in clusters, NOT in network  -> True Singleton, but why is it VC_303?
    # Mycobacterium~phage~Patt IN clusters, IN network  -> Overlap
    singletons = contigs - nodes
    merged_df.loc[merged_df["Genome"].isin(singletons), "VC Status"] = "Singleton"

    # Overlaps
    overlap_df = c1_df.loc[pd.notnull(c1_df["pos_clusters"])].copy()

    def update_VCs(string: str):
        pieces = string.split(";")
        updates = [int(piece) + 1 for piece in pieces]  # NO?
        str_fmt = "/".join(["VC_{}".format(update) for update in updates])

        return str_fmt

    overlap_df["pos_clusters"] = overlap_df["pos_clusters"].apply(update_VCs)
    overlap_strs = [
        "Overlap ({})".format(v) for (k, v) in overlap_df["pos_clusters"].items()
    ]
    merged_df.loc[
        merged_df["Genome"].isin(overlap_df.index.tolist()), "VC Status"
    ] = overlap_strs

    # Outliers = ntw - c1, shockingly, there's 0 'overlap' between outliers and overlaps
    # Nodes that made it to the network stage (i.e. not singletons) that WEREN'T in the c1 clusters file
    c1_members = nodes - set(c1_df[pd.notnull(c1_df["pos_cluster"])].index.tolist())

    # Nodes = all genomes in network
    merged_df.loc[merged_df["Genome"].isin(c1_members), "VC Status"] = "Outlier"

    # Filter to remove non-essential data
    # Is it present?
    taxonomies = [
        taxon
        for taxon in ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
        if taxon in merged_df.columns.tolist()
    ]
    merged_df = merged_df[["Genome", "VC Status"] + taxonomies]
    merged_df = merged_df[pd.notnull(merged_df["VC Status"])]

    return merged_df


def final_summary(
    folder,
    contigs: pd.DataFrame,
    network,
    profiles_fp,
    viral_clusters: vcontact2.cluster_refinements.ViralClusters,
    excluded: pd.DataFrame,
):
    node_table = contigs.copy()

    columns = [
        "VC",
        "Size",
        "Internal Weight",
        "External Weight",
        "Quality",
        "P-value",
        "Min Dist",
        "Max Dist",
        "Total Dist",
        "Below Thres",
        "Taxon Prediction Score",
        "Avg Dist",
        "Members",
    ]

    taxonomy_ranks = [
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "subfamily",
        "genus",
    ]

    incl_taxonomy = [
        incl_taxon
        for incl_taxon in taxonomy_ranks
        if incl_taxon in node_table.columns.tolist()
    ]

    columns_wTaxon = columns + incl_taxonomy

    summary_df = pd.DataFrame(columns=columns_wTaxon)

    num_contigs = len(node_table)

    logger.info(f"Reading edges for {num_contigs} contigs")
    edges_df = pd.read_csv(
        network,
        header=None,
        index_col=None,
        delimiter=" ",
        names=["source", "target", "weight"],
    )

    # Reinforce contigs with taxonomy, and sort of overlapping clusters
    # levels = [column for column in contigs.columns if column in taxonomy_ranks]

    if incl_taxonomy:
        node_table[incl_taxonomy] = node_table[incl_taxonomy].fillna("Unassigned")

    # Build PC array
    logger.info("Building PC array")
    profiles_df = pd.read_csv(profiles_fp, header=0)

    # Get number of comparisons
    logger.info("Calculating comparisons for back-calculations")
    pseudo_matrix_counts: pd.Series = (
        profiles_df["pc_id"].value_counts().pipe(lambda s: s[s >= 3])
    )
    pseudo_matrix = profiles_df[profiles_df["pc_id"].isin(pseudo_matrix_counts.index)]
    # Length of pseudo matrix is the total number of PCs
    # pcs = len(pseudo_matrix)
    # T = 0.5 * pcs * (pcs - 1)
    # logT = np.log10(T)

    # Keep track of genomes that get clustered, but get "excluded" when their dist is greater than threshold
    clustered_singletons = {}
    for contig_cluster, contig_cluster_group in node_table.groupby(
        by="rev_pos_cluster"
    ):
        logger.debug("Processing viral cluster {}".format(contig_cluster))

        size = len(contig_cluster_group)

        cluster_contigs = (
            contig_cluster_group["contig_id"].unique().tolist()
        )  # WARN network is ~
        # selected_edges = edges_df[
        #     (edges_df['source'].isin(cluster_contigs)) & (edges_df['target'].isin(cluster_contigs))]

        # Taxonomies
        taxonomies: dict[str, list] = dict.fromkeys(incl_taxonomy, [])

        for level in incl_taxonomy:
            # paired with above, but could replace/substitute as they're the same
            if level in contig_cluster_group:
                taxonomies[level] = [
                    item
                    for item in contig_cluster_group[level].unique()
                    if not pd.isnull(item)
                ]

        vc_pc_df = profiles_df.loc[
            profiles_df["contig_id"].isin(
                [contig.replace("~", " ") for contig in cluster_contigs]
            )
        ].copy()

        crosstab = pd.crosstab(vc_pc_df["contig_id"], vc_pc_df["pc_id"])

        try:
            dist = distance.pdist(crosstab.values, metric="euclidean")
            row_linkage = linkage(dist, method="average")

            # Keep all "distance" logic here
            c, coph_dists = cophenet(row_linkage, dist)
            logger.debug("Cophenet distance: {}".format(c))

            dist_size = dist.size
            average_dist = dist.mean()
            min_dist = dist.min()  # np.min(dist[np.nonzero(dist)])  # Still want zeroes
            max_dist = dist.max()
            thres_counts = dist[dist < viral_clusters.dist].size  # np.where(dist < 8)
            frac = float(thres_counts) / dist.size
            clustered_singletons.update(
                {contig: "Clustered" for contig in cluster_contigs}
            )

        # Empty distance matrix - occurs when overlapping members leaves a cluster w/ 1 member
        except ValueError:
            dist_size = 1
            min_dist = 0
            max_dist = 0
            average_dist = 0
            thres_counts = 1  # These are newly established "singetons"
            frac = 1
            clustered_singletons[cluster_contigs[0]] = "Clustered/Singleton"

        # Get internals
        intracluster_selected_edges = edges_df[
            (edges_df["source"].isin(cluster_contigs))
            & (edges_df["target"].isin(cluster_contigs))
        ].copy()

        # ClusterONE documentation states that the sum of the internal weights and sum of external weights, but it
        # appears that it's ONLY for non-duplicate edges, so need to remove duplicates
        # https://stackoverflow.com/q/30689236
        mask = (
            intracluster_selected_edges["source"]
            < intracluster_selected_edges["target"]
        )
        intracluster_selected_edges["first"] = intracluster_selected_edges[
            "source"
        ].where(mask, intracluster_selected_edges["target"])
        intracluster_selected_edges["second"] = intracluster_selected_edges[
            "target"
        ].where(mask, intracluster_selected_edges["source"])

        intra_df = intracluster_selected_edges.drop_duplicates(
            subset=["weight", "first", "second"]
        ).drop(columns=["first", "second"])

        # Once duplicates are removed, get simple list of weights and determine stats
        internal_scores = sorted(intra_df["weight"].tolist())
        internal_weights = sum(internal_scores)

        extracluster_selected_edges = edges_df[
            (edges_df["source"].isin(cluster_contigs))
            | (edges_df["target"].isin(cluster_contigs))
        ]
        # I seriously cannot filter this with double ~, wth
        extracluster_selected_edges = extracluster_selected_edges[
            ~(
                extracluster_selected_edges["source"].isin(cluster_contigs)
                & extracluster_selected_edges["target"].isin(cluster_contigs)
            )
        ]

        mask2 = (
            extracluster_selected_edges["source"]
            < extracluster_selected_edges["target"]
        )
        extracluster_selected_edges["first"] = extracluster_selected_edges[
            "source"
        ].where(mask2, extracluster_selected_edges["target"])
        extracluster_selected_edges["second"] = extracluster_selected_edges[
            "target"
        ].where(mask2, extracluster_selected_edges["source"])

        extra_df = extracluster_selected_edges.drop_duplicates(
            subset=["weight", "first", "second"]
        ).drop(columns=["first", "second"])

        external_scores = sorted(extra_df["weight"].tolist())
        external_weights = sum(external_scores)

        try:
            quality = (internal_weights + 0.0) / (internal_weights + external_weights)
        except ZeroDivisionError:
            quality = 1

        try:  # Less will give p-value 1.0 when it's 0
            stat, pval = mannwhitneyu(internal_scores, external_scores)
        except ValueError:  # All numbers are identical
            pval = 1

        try:
            pos = len(summary_df)
            summary_df.loc[pos, columns] = (
                f"VC_{contig_cluster}",
                size,
                internal_weights,
                external_weights,
                quality,
                pval,
                min_dist,
                max_dist,
                dist_size,
                thres_counts,
                frac,
                average_dist,
                ",".join(cluster_contigs),
            )
            for level, items in taxonomies.items():
                summary_df.loc[pos, level] = items
                # It is nice knowing what their components are...

        except Exception as e:
            logger.error(e)
            exit(1)

    logger.info("Writing viral cluster overview file...")
    summary_df.to_csv(os.path.join(folder, "viral_cluster_overview.csv"))

    node_table["Viral Cluster"] = node_table["rev_pos_cluster"].astype(str)

    summary_df["items"] = summary_df["VC"].apply(lambda x: len(x.split("_")))
    summary_df = summary_df[summary_df["items"] > 2]

    translator = {
        "kingdom": "VC Kingdoms",
        "phylum": "VC Phyla",
        "class": "VC Classes",
        "order": "VC Orders",
        "family": "VC Families",
        "genus": "VC Genera",
    }

    node_columns = [
        "Genome",
        *(str(incl_taxon).capitalize() for incl_taxon in incl_taxonomy),
        "preVC",
        "VC Status",
        "VC",
        "VC Size",
        "Quality",
        "Adjusted P-value",
        "VC Avg Distance",
        "Topology Confidence Score",
        "Genus Confidence Score",
        *(map(translator.get, incl_taxonomy)),
    ]
    # Discontinue 'preVC Size' as it's no longer relevant

    node_summary_df = pd.DataFrame(columns=node_columns)

    # Iterate through final table to summary everything
    logger.info(
        f"Examining each viral cluster and breaking it down into individual genomes..."
    )
    for index, vc_df in node_table.iterrows():
        node_summary_pos = len(node_summary_df)  # Keep everything on same line

        # Minimum information - even if matches to nothing
        genome = vc_df["contig_id"]
        node_summary_df.loc[node_summary_pos, "Genome"] = genome

        for incl_taxon in incl_taxonomy:
            node_summary_df.loc[node_summary_pos, incl_taxon.capitalize()] = vc_df[
                incl_taxon
            ]  # Shouldn't need to test

        vc: str = vc_df["Viral Cluster"]

        try:
            if pd.isna(float(vc)):
                # node_summary_df.loc[node_summary_pos, 'VC Status'] = 'Unassigned (No VC)'
                continue

        except AttributeError:
            pass

        # Get VC information
        # May produce byproduct (>1, below)
        genome_df = summary_df.loc[
            summary_df["Members"].str.contains(genome, regex=False)
        ]

        if len(genome_df) == 0:
            # When VC is 'nan'  # Should this even occur if viral cluster is nan?
            continue

        if len(genome_df) > 1:
            # When genome name is subset of another... is this actually worth separating for here?
            for sub_index, sub_df in genome_df.iterrows():
                if genome in sub_df["Members"].split(","):
                    genome_df = summary_df.loc[
                        summary_df["Members"] == sub_df["Members"]
                    ]
        # If len still >1, likely not the genome we want...
        if len(genome_df) > 1:
            # Genome name is a string SUBSET of other members, but is not identical
            # Bacillus~virus~G vs Bacillus~virus~Glittering and Bacillus~virus~GA1
            logger.warning(
                f"Could not identify genome substrings for {genome}. Consider adjusting input "
                f"genomes naming so genome names aren't potential substrings of each other."
            )
            continue

        genome_s = genome_df.iloc[0]

        preVC = f'preVC_{vc.split("_")[0]}'  # Should non

        # Everything is assigned to a subcluster, meaning summary table will have dup cols
        subcluster = genome_s["VC"]
        subcluster_size = len(
            summary_df.loc[summary_df["VC"] == subcluster]["Members"].item().split(",")
        )
        quality = genome_s["Quality"]
        adj_pval = 1.0 - float(genome_s["P-value"])
        # cohesiveness = min(400, np.nan_to_num(-np.log10(genome_s['Cohesiveness'])))
        vc_avg_dist = genome_s["Avg Dist"]
        vc_conf = genome_s["Taxon Prediction Score"]
        vc_overall_conf = quality * adj_pval

        general_stats_cols = [
            "preVC",
            "VC",
            "VC Size",
            "Quality",
            "Adjusted P-value",
            "VC Avg Distance",
            "Topology Confidence Score",
        ]  # 'Genus Confidence Score'
        node_summary_df.loc[node_summary_pos, general_stats_cols] = (
            preVC,
            subcluster,
            subcluster_size,
            quality,
            adj_pval,
            vc_avg_dist,
            vc_overall_conf,
        )

        # Add in taxonomy
        for incl_taxon in incl_taxonomy:
            counts = [taxon for taxon in genome_s[incl_taxon] if taxon != "Unassigned"]
            node_summary_df.loc[node_summary_pos, translator.get(incl_taxon)] = len(
                set(counts)
            )

        # Just missing genome confidence
        genus_conf = vc_conf
        node_summary_df.loc[node_summary_pos, "Genus Confidence Score"] = genus_conf
        # or genus_conf

        try:
            status = clustered_singletons[genome]

        except KeyError:
            if genome in excluded["Genome"].tolist():
                continue  # It'll be handled by append later
            else:
                status = "Unavailable"

            logger.warning(
                f"There was an error during the handling of: {genome}\n"
                f"This is almost certainly due to being unable to identify {genome}'s VC status."
            )

        try:
            node_summary_df.loc[node_summary_pos, "VC Status"] = status

        except Exception as e:
            logger.error(
                f"There was another error during the handling of {genome}: {e}"
            )
            exit(1)

    node_summary_df["Quality"] = node_summary_df["Quality"].apply(lambda x: round(x, 4))
    node_summary_df["Adjusted P-value"] = node_summary_df["Adjusted P-value"].apply(
        lambda x: round(x, 8)
    )
    node_summary_df["Genus Confidence Score"] = node_summary_df[
        "Genus Confidence Score"
    ].apply(lambda x: round(x, 4))
    node_summary_df["Topology Confidence Score"] = node_summary_df[
        "Topology Confidence Score"
    ].apply(lambda x: round(x, 4))

    node_summary_df.set_index("Genome", inplace=True, drop=True)
    excluded.set_index("Genome", inplace=True, drop=True)

    node_summary_df.update(
        excluded, overwrite=False
    )  # Only want it to update NaN, i.e. No status

    node_summary_df.sort_index(ascending=True, inplace=True)

    logger.info(f"Writing the genome-by-genome overview file...")
    node_summary_df.to_csv(os.path.join(folder, "genome_by_genome_overview.csv"))
