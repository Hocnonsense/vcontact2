"""Evaluations : A series analytics to evaluate performance of the network and resulting taxonomy"""

import logging
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class Evaluations(object):
    """
    Collects series of functions related to analyzing the network.

    Attributes:

    """

    def __init__(self, contigs, levels=None, focus='pos_cluster'):
        """
        :param contigs: (dataframe)
        """

        self.name = "Evaluations"
        self.focus = focus

        # Check if multiple taxonomic levels exist, and if so, process each level
        self.tax_metrics = {}
        if levels:
            self.levels = levels
        else:
            predefined = ['order', 'family', 'subfamily', 'genus']
            self.levels = [column for column in contigs.columns if column in predefined]
            logger.debug("{} taxonomic levels detected: {}".format(len(self.levels), ", ".join(self.levels)))

        for level in self.levels:
            logger.info("Performance evaluations at the {} level...".format(level))

            tmp_contigs_df = contigs.drop_duplicates(['contig_id'], keep='first')

            contingency_tbl = pd.crosstab(tmp_contigs_df[level], tmp_contigs_df[self.focus])

            clustering_wise_ppv, clustering_wise_sensitivity, accuracy = self.performance_metrics(contingency_tbl)

            self.tax_metrics[level] = {
                'PPV': clustering_wise_ppv,
                'Sensitivity': clustering_wise_sensitivity,
                'Accuracy': accuracy
            }

    def performance_metrics(self, contingency):
        # Calculations will assume 0's need to be counted if they exist. They technically do, but shouldn't be included
        # in the calculation
        contingency_table = contingency.replace({0: np.nan})
        # contingency_complex_weights = contingency_table.sum(axis=1)
        # contingency_cluster_weights = contingency_table.sum(axis=0)

        sensitivity_tbl = self.calc_sensitivity(contingency_table)

        ppv_tbl = self.calc_ppv(contingency_table)
        ppv_tbl = ppv_tbl.loc[ppv_tbl.index != 'unassigned']

        # Maximum "complex-wise" PPV for each cluster
        ppv_tbl.loc['cluster-wise PPV'] = ppv_tbl.max(axis=0)

        # Weighted average of cluster-wise PPV
        # https://stackoverflow.com/questions/21113384/python-numpy-average-with-nans
        ppv_indices = ~np.isnan(ppv_tbl.loc['cluster-wise PPV'].values)
        clustering_wise_ppv = np.average(
            ppv_tbl.loc['cluster-wise PPV'].values[ppv_indices])  # , weights=contingency_cluster_weights)

        # Maximum "cluster-wise" sensitivity for each complex (new column)
        sensitivity_tbl['complex-wise sensitivity'] = sensitivity_tbl.max(axis=1)

        sensitivity_tbl = sensitivity_tbl[sensitivity_tbl.index != 'unassigned']

        # Weighted average of complex-wise sensitivity
        sens_indices = ~np.isnan(sensitivity_tbl['complex-wise sensitivity'].values)
        clustering_wise_sensitivity = np.average(sensitivity_tbl['complex-wise sensitivity'].values[sens_indices])

        accuracy = self.geo_mean([clustering_wise_ppv, clustering_wise_sensitivity])

        return clustering_wise_ppv, clustering_wise_sensitivity, accuracy

    def calc_sensitivity(self, contingency_table):

        """
        Fraction of proteins in complex i found in cluster j
        Complex-wise sensitivity isn't calculated here, but is simply max value for each complex
        :param contingency_table: taxon level X VC membership matrix
        :return: sensitivity_tbl
        """

        # Get sum of rows (adds new column)
        contingency_counts = contingency_table.copy()
        contingency_counts['sum'] = contingency_counts.sum(axis=1)

        # Each cluster/complex divided by the total sum of each complex
        sensitivity = contingency_table.div(contingency_counts['sum'], axis=0)

        return sensitivity

    def calc_ppv(self, contingency_table):
        """
        Proportion of members of cluster j which belong to complex i, relative to total number of members
        :param contingency_tbl: taxon level X VC membership matrix
        :return:
        """

        # Get sum of columns (new row)
        counts = contingency_table.copy()
        counts.loc['sum'] = counts.sum(axis=0)  # adds a new ROW with sum of column

        ppv = contingency_table.div(counts.loc['sum'], axis=1)

        return ppv

    def calc_accuracy(self, sensitivity_table, ppv_tbl):
        """

        :param sensitivity_tbl:
        :param ppv_tbl:
        :return:
        """
        def geo_mean_overflow(a, b):
            res = np.log([a, b])
            return np.exp(res.sum() / len(res))

        vec_geo_mean = np.vectorize(geo_mean_overflow)

        accuracy = pd.DataFrame(vec_geo_mean(sensitivity_table, ppv_tbl))
        accuracy.columns = sensitivity_table.columns
        accuracy.index = sensitivity_table.index

        return accuracy

    def geo_mean(self, iterable):
        a = np.array(iterable)
        return a.prod() ** (1.0 / len(a))

