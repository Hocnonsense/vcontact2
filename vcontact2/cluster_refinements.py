"""Cluster Refinements : Refining contig clusters for optimal assignments"""

import logging

import numpy as np
import pandas as pd
from scipy.spatial import distance
import scipy.cluster as sclust
from scipy.cluster.hierarchy import linkage
import vcontact2.evaluations

logger = logging.getLogger(__name__)

default_columns = [
    "contig_id",
    "proteins",
    "size",
    "origin",
    "pos",
    "pos_cluster",
    "membership",
    "pos_clusters",
    "index",
    "pos_cluster_mbship",
]

summary_headers = [
    "Distance",
    "Clustering-wise sensitivity",
    "Clustering-wise PPV",
    "Accuracy",
]


class ViralClusters(object):
    """
    Collects series of functions related to analyzing the network.

    Attributes:

    """

    def __init__(self, contigs: pd.DataFrame, profiles, optimize=False):
        """
        :param contigs: (dataframe)
        """
        self.name = "ViralClusters"

        profiles_df = pd.read_csv(profiles, header=0)
        # Build PC array
        self.metrics = pd.DataFrame(columns=summary_headers)
        self.results = {}

        self.levels = frozenset(contigs.columns) - frozenset(default_columns)
        logger.debug(
            "{} taxonomic levels detected: {}".format(
                len(self.levels), ", ".join(self.levels)
            )
        )

        if len(self.levels) == 0:
            logger.info(
                "Unable to identify any taxonomic levels. "
                "Was a reference database included? If not, then no worries"
            )

        if optimize and len(self.levels) != 0:
            # dists = np.linspace(1, 20.0, 39, endpoint=True).tolist()  # 1, 20.0, 191 = 0.1,
            dists = [
                round(dist, 2)
                for dist in np.linspace(1, 20.0, 20, endpoint=True).tolist()
            ]
        else:
            # By defining dists here, don't need to repeat performance metrics after the loop
            dists = [9]
        for dist in dists:
            logger.info("Optimizing on distance: {}".format(dist))

            adj_contigs = contigs.copy()

            for contig_cluster, contig_cluster_group in adj_contigs.groupby(
                by="pos_cluster"
            ):
                vc_pc_df = profiles_df.loc[
                    profiles_df["contig_id"].isin(
                        [
                            contig.replace("~", " ")
                            for contig in contig_cluster_group["contig_id"]
                        ]
                    )
                ].copy()

                crosstab = pd.crosstab(vc_pc_df["contig_id"], vc_pc_df["pc_id"])

                # Need to set up appropriate clusters
                for n, label in enumerate(crosstab.index.tolist()):
                    vc_pc_df.loc[vc_pc_df["contig_id"] == label, "unique_id"] = str(n)

                try:
                    row_linkage = linkage(
                        distance.pdist(crosstab.values), method="average"
                    )
                except ValueError:
                    # These are VCs whose OTHER MEMBERS are overlapping, meaning they're the ONLY remaining
                    # and you can't calculate a pdist with only 1 member
                    adj_contigs.loc[
                        contigs["pos_cluster"] == contig_cluster, "rev_pos_cluster"
                    ] = "{}_0".format(contig_cluster)
                    continue

                # Get clusters
                for n, fcluster in enumerate(
                    sclust.hierarchy.fcluster(row_linkage, dist, criterion="distance")
                ):
                    vc_pc_df.loc[vc_pc_df["unique_id"] == str(n), "fcluster"] = str(
                        fcluster
                    )

                # Previously got number of subclusters prior, then assigned them,
                # but if starting from 0 everytime, there's no point
                for n, (fcluster, fcluster_df) in enumerate(
                    vc_pc_df.groupby(by="fcluster")
                ):
                    vc_pc_df.loc[
                        fcluster_df.index, "rev_pos_cluster"
                    ] = f"{contig_cluster}_{n}"

                    contigs_items = [
                        item.replace(" ", "~")
                        for item in fcluster_df["contig_id"].unique().tolist()
                    ]
                    adj_contigs.loc[
                        adj_contigs["contig_id"].isin(contigs_items), "rev_pos_cluster"
                    ] = f"{contig_cluster}_{n}"

                    # No need to get subdistance, as it's already < threshold
                    # subtab = pd.crosstab(fcluster_df['contig_id'], fcluster_df['pc_id'])
                    # subdist = distance.pdist(subtab.values, metric='euclidean')
                    # subdist_avg = subdist.mean()

            self.results[dist] = adj_contigs

            # Performance metrics
            if "genus" in self.levels:  # If there's actual taxonomy to optimize
                evaluations = vcontact2.evaluations.Evaluations(
                    adj_contigs, levels=["genus"], focus="rev_pos_cluster"
                )

                self.metrics.loc[len(self.metrics), summary_headers] = (
                    dist,
                    evaluations.tax_metrics["genus"]["Sensitivity"],
                    evaluations.tax_metrics["genus"]["PPV"],
                    evaluations.tax_metrics["genus"]["Accuracy"],
                )
            else:
                self.metrics.loc[len(self.metrics), summary_headers] = dist, 0, 0, 0

        # Find best composite score
        self.metrics["Composite Score"] = self.metrics[
            ["Clustering-wise sensitivity", "Clustering-wise PPV", "Accuracy"]
        ].sum(axis="columns")
        self.best_score = self.metrics["Composite Score"].max()
        self.best_df = self.metrics.loc[
            self.metrics["Composite Score"] == self.best_score
        ]

        if len(self.best_df) == 1:
            self.dist = self.best_df.iloc[0]["Distance"]
            logger.info(
                "Identified a single best composite score {} for distance {}".format(
                    self.best_score, self.dist
                )
            )
        elif len(self.best_df) == 2:
            best_index = self.best_df["Distance"].tolist()
            logger.info(
                "Identified the best composite scores among two distances, "
                "selecting the larger distance: {}".format(",".join(best_index))
            )
            self.dist = self.best_df.iloc[best_index[-1]]["Distance"]
        elif len(self.best_df) > 2:
            best_index = self.best_df["Distance"].tolist()
            logger.warning(
                "Identified best composite scores among multiple distances! Optimal distance may not be "
                "calculated correctly due to a small sample size, heterogeneity in the data, or some "
                "issue with the dataset. Selecting the largest distance: {}".format(
                    ",".join([str(i) for i in best_index])
                )
            )
            self.dist = self.best_df.iloc[best_index[-1]]["Distance"]

        logger.info("Merging optimal distance determined from performance evaluations.")

        self.contigs = self.results[self.dist]

        self.performance = pd.DataFrame.from_dict(
            vcontact2.evaluations.Evaluations(
                self.contigs, focus="rev_pos_cluster"
            ).tax_metrics,
            orient="index",
        )

        logger.info(self.performance)
