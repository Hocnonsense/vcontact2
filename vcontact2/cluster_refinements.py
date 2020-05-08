"""Cluster Refinements : Refining contig clusters for optimal assignments"""

import logging
from pprint import pprint

import numpy as np
import pandas as pd
from scipy.spatial import distance
import scipy.cluster as sclust
from scipy.cluster.hierarchy import linkage, cophenet
import vcontact2.evaluations

logger = logging.getLogger(__name__)


class ViralClusters(object):
    """
    Collects series of functions related to analyzing the network.

    Attributes:

    """

    def __init__(self, contigs, profiles, optimize=False):
        """

        :param contigs: (dataframe)
        """

        self.name = "ViralClusters"
        self.contigs = contigs

        # Build PC array
        profiles_df = pd.read_csv(profiles, header=0)

        summary_headers = ['Distance', 'Clustering-wise sensitivity', 'Clustering-wise PPV', 'Accuracy']
        self.metrics = pd.DataFrame(columns=summary_headers)
        self.results = {}
        self.dist = 9

        default_columns = ["contig_id", "proteins", "size", "origin", "pos", "pos_cluster", "membership",
                           "pos_clusters", "index", "pos_cluster_mbship"]
        self.levels = frozenset(self.contigs.columns) - frozenset(default_columns)
        logger.debug("{} taxonomic levels detected: {}".format(len(self.levels), ", ".join(self.levels)))

        if len(self.levels) == 0:  #
            logger.info('Unable to identify any taxonomic levels. Was a reference database included? If not, then no'
                        'worries')

        if optimize and len(self.levels) != 0:

            # dists = np.linspace(1, 20.0, 39, endpoint=True).tolist()  # 1, 20.0, 191 = 0.1,
            dists = np.linspace(1, 20.0, 20, endpoint=True).tolist()
            dists = [round(dist, 2) for dist in dists]
        else:
            dists = [self.dist]  # By defining dists here, don't need to repeat performance metrics after the loop

        for dist in dists:

            logger.info('Optimizing on distance: {}'.format(dist))

            adj_contigs = self.contigs.copy()

            for contig_cluster, contig_cluster_group in adj_contigs.groupby(by='pos_cluster'):

                vc_pc_df = profiles_df.loc[profiles_df['contig_id'].isin([contig.replace('~', ' ') for contig in contig_cluster_group['contig_id']])].copy()

                crosstab = pd.crosstab(vc_pc_df['contig_id'], vc_pc_df['pc_id'])

                # Need to set up appropriate clusters
                for n, label in enumerate(crosstab.index.tolist()):
                    vc_pc_df.loc[vc_pc_df['contig_id'] == label, 'unique_id'] = str(n)

                try:
                    row_linkage = linkage(distance.pdist(crosstab.values), method='average')
                except ValueError:  # These are VCs whose OTHER MEMBERS are overlapping, meaning they're the ONLY remaining
                    # and you can't calculate a pdist with only 1 member
                    adj_contigs.loc[contigs['pos_cluster'] == contig_cluster, 'rev_pos_cluster'] = '{}_0'.format(contig_cluster)
                    continue

                pdist = distance.pdist(crosstab.values, metric='euclidean')
                average_dist = pdist.mean()

                # Get clusters
                fclusters = sclust.hierarchy.fcluster(row_linkage, dist, criterion='distance')
                for n, fcluster in enumerate(fclusters):
                    vc_pc_df.loc[vc_pc_df['unique_id'] == str(n), 'fcluster'] = str(fcluster)

                # Previously got number of subclusters prior, then assigned them, but if starting from 0 everytime, there's no point
                for n, (fcluster, fcluster_df) in enumerate(vc_pc_df.groupby(by='fcluster')):
                    vc_pc_df.loc[fcluster_df.index, 'rev_pos_cluster'] = '{}_{}'.format(contig_cluster, n)

                    contigs_items = [item.replace(' ', '~') for item in fcluster_df['contig_id'].unique().tolist()]
                    adj_contigs.loc[adj_contigs['contig_id'].isin(contigs_items), 'rev_pos_cluster'] = '{}_{}'.format(
                        contig_cluster, n)

                    # No need to get subdistance, as it's already < threshold
                    # subtab = pd.crosstab(fcluster_df['contig_id'], fcluster_df['pc_id'])
                    # subdist = distance.pdist(subtab.values, metric='euclidean')
                    # subdist_avg = subdist.mean()

            self.results[dist] = adj_contigs

            # Performance metrics
            if 'genus' in self.levels:  # If there's actual taxonomy to optimize
                evaluations = vcontact2.evaluations.Evaluations(adj_contigs, levels=['genus'], focus='rev_pos_cluster')
                ppv = evaluations.tax_metrics['genus']['PPV']
                sensitivity = evaluations.tax_metrics['genus']['Sensitivity']
                accuracy = evaluations.tax_metrics['genus']['Accuracy']

                self.metrics.loc[len(self.metrics), summary_headers] = dist, sensitivity, ppv, accuracy
            else:
                self.metrics.loc[len(self.metrics), summary_headers] = dist, 0, 0, 0

        # Find best composite score
        self.metrics['Composite Score'] = self.metrics[['Clustering-wise sensitivity', 'Clustering-wise PPV',
                                                        'Accuracy']].sum(axis='columns')
        self.best_score = self.metrics['Composite Score'].max()
        self.best_df = self.metrics.loc[self.metrics['Composite Score'] == self.best_score]

        if len(self.best_df) == 1:
            self.dist = self.best_df.iloc[0]['Distance']
            logger.info('Identified a single best composite score {} for distance {}'.format(
                self.best_score, self.dist))

        if len(self.best_df) == 2:
            best_index = self.best_df['Distance'].tolist()
            logger.info('Identified the best composite scores among two distances, '
                        'selecting the larger distance: {}'.format(','.join(best_index)))
            self.dist = self.best_df.iloc[best_index[-1]]['Distance']

        if len(self.best_df) > 2:
            best_index = self.best_df['Distance'].tolist()
            logger.warning('Identified best composite scores among multiple distances! Optimal distance may not be '
                           'calculated correctly due to a small sample size, heterogeneity in the data, or some '
                           'issue with the dataset. Selecting the largest distance: {}'.format(','.join([str(i) for i in best_index])))
            self.dist = self.best_df.iloc[best_index[-1]]['Distance']

        logger.info('Merging optimal distance determined from performance evaluations.')

        self.contigs = self.results[self.dist]

        self.performance = pd.DataFrame.from_dict(vcontact2.evaluations.Evaluations(
            self.contigs, focus='rev_pos_cluster').tax_metrics, orient='index')

        logger.info(self.performance)

