""" Functions to export final node tables that includes VC and individual contigs """
import logging
import os
import pandas as pd
import numpy as np
import math
from functools import reduce
from scipy.spatial import distance
from scipy.stats import mannwhitneyu
from scipy.cluster.hierarchy import linkage, cophenet
np.warnings.filterwarnings('ignore')

from pprint import pprint

logger = logging.getLogger(__name__)


def find_excluded(merged, ntw, c1):
    """
    Temporary placement of singleton, overlap and outlier detection. These three components should be refactored
    "closer" to the locations where they'd be initially generated. Overlap genomes *are already identified* and saved
    in a dataframe, but that's only for overlapping genome removal, the dataframe/overlaps are never saved to a file.
    The outliers/singletons require a multi-step process, as ClusterONE does not state which genomes are outliers, so
    they must be back-calculated.

    singletons = merged - ntw
    outliers = ntw - c1
    overlap = c1 overlaps

    :return:
    """

    # Get list of all genomes
    merged_df = pd.read_csv(merged, header=0, index_col=0)
    merged_df.rename(columns={'order': 'Order', 'family': 'Family', 'genus': 'Genus', 'contig_id': 'Genome'},
                     inplace=True)
    contigs = set(merged_df['Genome'].tolist())

    merged_df['VC Status'] = np.nan

    # Get list of genomes that made it to the network, those that didn't did not pass the thresholds and => singletons
    network_df = pd.read_csv(ntw, header=None, names=['source', 'target', 'weight'], index_col=None, delimiter=' ')
    nodes = set(network_df['source'].tolist() + network_df['target'].tolist())

    # Set up final destination for singletons
    singletons = contigs - nodes
    merged_df.loc[merged_df['Genome'].isin(singletons), 'VC Status'] = 'Singleton'

    # Overlaps requite re-reading the clusters file, because of course, I don't save it earlier
    c1_df = pd.read_csv(c1, header=0)

    # Need to adjust cluster num to match pos so it works, but doesn't need to be done here
    c1_df['Cluster'] = c1_df['Cluster'].astype(str)

    # Don't need to go through and re-name and format because we already have all members (from contigs)
    labels = {}

    for cluster, cluster_df in c1_df.groupby(by='Cluster'):  # Either this or 2 lists + enumerate
        formatter = "VC_{}".format(cluster)
        members = cluster_df['Members'].iloc[0].rstrip().split()
        for member in members:
            if member not in labels:
                labels[member] = formatter
            else:
                labels[member] = ';'.join([labels[member], formatter])

    # Overlaps. Can just use key (genome name) and set the Destination to Overlap, but it's just a little more work to
    # actually give the clusters that overlapped
    overlaps = [(k, v) for k, v in labels.items() if ';' in v]
    overlap_members = [k for (k, v) in overlaps]
    overlap_clusters = ['Overlap ({})'.format(str(v).replace(';', '/')) for (k, v) in overlaps]

    # merged_df.loc[merged_df['Genome'].isin(overlap_members), 'VC Status'] = overlap_clusters

    # Outliers = ntw - c1, shockingly, there's 0 'overlap' between outliers and overlaps
    outliers = nodes - set(labels.keys())  # Labels = all genomes in clusterONE output, regardless of clustering status
    # Nodes = all genomes in network
    merged_df.loc[merged_df['Genome'].isin(outliers), 'VC Status'] = 'Outlier'

    # Filter to remove non-essential data
    # Is it present?
    taxonomies = [taxon for taxon in ['Order', 'Family', 'Genus'] if taxon in merged_df.columns.tolist()]
    merged_df = merged_df[['Genome', 'VC Status'] + taxonomies]
    merged_df = merged_df[pd.notnull(merged_df['VC Status'])]

    return merged_df


def final_summary(folder, contigs, network, profiles, viral_clusters, excluded):

    node_table = contigs.copy()

    columns = ['VC', 'Size', 'Internal Weight', 'External Weight', 'Quality', 'P-value', 'Min Dist',
               'Max Dist', 'Total Dist', 'Below Thres', 'Taxon Prediction Score', 'Avg Dist', 'Genera', 'Families',
               'Orders', 'Members']
    summary_df = pd.DataFrame(columns=columns)

    taxon_columns = ['Level', 'Taxon', 'Min Dist', 'Max Dist', 'Total Dist', 'Below Thres', 'Taxon Prediction Score',
                     'Avg Dist']
    taxonomy_df = pd.DataFrame(columns=taxon_columns)

    num_contigs = len(node_table)

    logger.info('Reading edges for {} contigs'.format(num_contigs))
    edges_df = pd.read_csv(network, header=None, index_col=None, delimiter=' ', names=['source', 'target', 'weight'])

    # Reinforce contigs with taxonomy, and sort of overlapping clusters
    # print('Merging nodes with reference taxonomy')
    # node_table_df = node_table_df.merge(contigs_tax_df, how='left', left_on='Name', right_on='contig_id')
    predefined = ['order', 'family', 'subfamily', 'genus']
    levels = [column for column in contigs.columns if column in predefined]

    if levels:
        node_table[levels] = node_table[levels].fillna('Unassigned')

    # Build PC array
    logger.info('Building PC array')
    profiles_df = pd.read_csv(profiles, header=0)

    # Get number of comparisons
    logger.info('Calculating comparisons for back-calculations')
    pseudo_matrix_counts = profiles_df['pc_id'].value_counts()
    pseudo_matrix_counts = pseudo_matrix_counts[pseudo_matrix_counts >= 3]
    pseudo_matrix = profiles_df[profiles_df['pc_id'].isin(pseudo_matrix_counts.index)]
    # Length of pseudo matrix is the total number of PCs
    pcs = len(pseudo_matrix)
    T = 0.5 * pcs * (pcs - 1)
    logT = np.log10(T)

    for contig_cluster, contig_cluster_group in node_table.groupby(by='rev_pos_cluster'):
        logger.debug('Processing viral cluster {}'.format(contig_cluster))

        size = len(contig_cluster_group)

        cluster_contigs = contig_cluster_group['contig_id'].unique()  # WARN network is ~
        selected_edges = edges_df[
            (edges_df['source'].isin(cluster_contigs)) & (edges_df['target'].isin(cluster_contigs))]

        # Oh math, logs additive if to same base
        # sig_scores = selected_edges['weight'].tolist()
        # raw_pvals = []
        # for sig_score in sig_scores:
        #     h_pval = sig_score + logT
        #     z = math.pow(10, -h_pval)  # 3.348269070478552e+66
        #
        #     raw_pvals.append(z)
        # # https://stackoverflow.com/questions/13840379/how-can-i-multiply-all-items-in-a-list-together-with-python
        # try:
        #     cohesiveness = reduce(lambda x, y: x * y, raw_pvals)
        # except Exception:
        #     cohesiveness = 0

        # Taxonomies
        taxonomies = {
            'order': [],
            'family': [],
            'genus': []
        }
        if levels:
            for level in levels:
                values = [item for item in contig_cluster_group[level].unique() if not pd.isnull(item)]
                taxonomies[level] = values

        vc_pc_df = profiles_df.loc[profiles_df['contig_id'].isin([contig.replace('~', ' ') for contig in cluster_contigs])].copy()

        crosstab = pd.crosstab(vc_pc_df['contig_id'], vc_pc_df['pc_id'])

        try:
            dist = distance.pdist(crosstab.values, metric='euclidean')
            row_linkage = linkage(dist, method='average')
        except ValueError:
            # print('VC {}'.format(contig_cluster))
            continue  # Empty distance matrix - occurs when overlapping members leaves a cluster w/ 1 member

        c, coph_dists = cophenet(row_linkage, dist)
        logger.debug('Cophenet distance: {}'.format(c))

        average_dist = dist.mean()
        min_dist = dist.min()  # np.min(dist[np.nonzero(dist)])  # Still want zeroes
        max_dist = dist.max()
        thres_counts = dist[dist < viral_clusters.dist].size  # np.where(dist < 8)
        frac = float(thres_counts) / dist.size

        # Get internals
        intracluster_selected_edges = edges_df[
            (edges_df['source'].isin(cluster_contigs)) & (edges_df['target'].isin(cluster_contigs))].copy()

        # ClusterONE documentation states that the sum of the internal weights and sum of external weights, but it
        # appears that it's ONLY for non-duplicate edges, so need to remove duplicates
        # https://stackoverflow.com/q/30689236
        mask = intracluster_selected_edges['source'] < intracluster_selected_edges['target']
        intracluster_selected_edges['first'] = intracluster_selected_edges['source'].where(mask,
                                                                                           intracluster_selected_edges[
                                                                                               'target'])
        intracluster_selected_edges['second'] = intracluster_selected_edges['target'].where(mask,
                                                                                            intracluster_selected_edges[
                                                                                                'source'])

        intra_df = intracluster_selected_edges.drop_duplicates(subset=['weight', 'first', 'second'])
        intra_df.drop(columns=['first', 'second'], inplace=True)

        # Once duplicates are removed, get simple list of weights and determine stats
        internal_scores = intra_df['weight'].tolist()
        internal_weights = sum(internal_scores)

        extracluster_selected_edges = edges_df[
            (edges_df['source'].isin(cluster_contigs)) | (edges_df['target'].isin(cluster_contigs))]
        # I seriously cannot filter this with double ~, wth
        extracluster_selected_edges = extracluster_selected_edges[~(
        extracluster_selected_edges['source'].isin(cluster_contigs) & extracluster_selected_edges['target'].isin(cluster_contigs))]

        mask2 = extracluster_selected_edges['source'] < extracluster_selected_edges['target']
        extracluster_selected_edges['first'] = extracluster_selected_edges['source'].where(mask2,
                                                                                           extracluster_selected_edges[
                                                                                               'target'])
        extracluster_selected_edges['second'] = extracluster_selected_edges['target'].where(mask2,
                                                                                            extracluster_selected_edges[
                                                                                                'source'])

        extra_df = extracluster_selected_edges.drop_duplicates(subset=['weight', 'first', 'second'])
        extra_df.drop(columns=['first', 'second'], inplace=True)

        external_scores = extra_df['weight'].tolist()
        external_weights = sum(external_scores)

        internal_scores = sorted(internal_scores)
        external_scores = sorted(external_scores)

        try:
            quality = (internal_weights + .0) / (internal_weights + external_weights)
        except ZeroDivisionError:
            quality = 1

        try:  # Less will give p-value 1.0 when it's 0
            stat, pval = mannwhitneyu(internal_scores, external_scores)
        except ValueError:  # All numbers are identical
            pval = 1

        try:
            summary_df.loc[len(summary_df), columns] = 'VC_{}'.format(contig_cluster), size, internal_weights, \
                                                       external_weights, quality, pval, min_dist, max_dist, \
                                                       dist.size, thres_counts, frac, average_dist, \
                                                       len(taxonomies['genus']), len(taxonomies['family']), \
                                                       len(taxonomies['order']), ','.join(cluster_contigs)
        except Exception as e:
            logger.error(e)
            exit(1)

    summary_df.to_csv(os.path.join(folder, 'viral_cluster_overview.csv'))

    node_table['Viral Cluster'] = node_table['rev_pos_cluster'].astype(str)

    summary_df['items'] = summary_df['VC'].apply(lambda x: len(x.split('_')))
    summary_df = summary_df[summary_df['items'] > 2]

    node_columns = ['Genome', 'Order', 'Family', 'Genus', 'VC', 'VC Status', 'Size', 'VC Subcluster',
                    'VC Subcluster Size', 'Quality', 'Adj P-value', 'Topology Confidence Score',
                    'Genera in VC', 'Families in VC', 'Orders in VC', 'Genus Confidence Score']
    node_summary_df = pd.DataFrame(columns=node_columns)

    # Iterate through final table to summary everything
    for index, vc_df in node_table.iterrows():

        genome = vc_df['contig_id']
        host = vc_df['Host'] if 'Host' in vc_df else 'Unknown'
        order = vc_df['order'] if 'order' in vc_df else 'Unassigned'
        family = vc_df['family'] if 'family' in vc_df else 'Unassigned'
        genus = vc_df['genus'] if 'genus' in vc_df else 'Unassigned'
        vc = vc_df['Viral Cluster']

        if pd.isnull(vc):
            node_summary_df.loc[len(node_summary_df), node_columns[:6]] = genome, order, family, genus, vc, 'Unassigned'
            continue

        # Get VC information
        genome_df = summary_df.loc[summary_df['Members'].str.contains(genome)]  # May produce byproduct (>1, below)

        if len(genome_df) == 0:  # When VC is 'nan'
            continue

        if len(genome_df) > 1:  # When genome name is subset of another
            for sub_index, sub_df in genome_df.iterrows():
                if genome in sub_df['Members'].split(','):
                    genome_df = summary_df.loc[summary_df['Members'] == sub_df['Members']]

        # If len still >1, likely not the genome we want...
        if len(genome_df) > 1:  # Genome name is a string SUBSET of other members, but is not identical
            # Bacillus~virus~G vs Bacillus~virus~Glittering and Bacillus~virus~GA1
            continue

        genome_s = genome_df.iloc[0]

        size = genome_s['Size']  # Everything is assigned to a subcluster, meaning summary table will have dup cols
        subcluster = genome_s['VC']
        subcluster_size = len(summary_df.loc[summary_df['VC'] == subcluster]['Members'].item().split(','))
        quality = genome_s['Quality']
        adj_pval = 1.0 - float(genome_s['P-value'])
        # cohesiveness = min(400, np.nan_to_num(-np.log10(genome_s['Cohesiveness'])))
        vc_avg_dist = genome_s['Avg Dist']
        vc_conf = genome_s['Taxon Prediction Score']
        vc_overall_conf = quality * adj_pval
        vc_genera = genome_s['Genera']
        vc_families = genome_s['Families']
        vc_orders = genome_s['Orders']

        genus_df = taxonomy_df.loc[(taxonomy_df['Level'] == 'genus') & (taxonomy_df['Taxon'] == genus)]
        genus_conf = False
        genus_dist = False
        if genus_df.empty:
            genus_conf = vc_conf
        else:
            genus_s = genus_df.iloc[0]
            genus_conf = genus_s['Taxon Prediction Score']

        try:
            node_summary_df.loc[len(node_summary_df), node_columns] = genome, order, family, genus, vc, 'Clustered',\
                                                                      size, subcluster, subcluster_size, quality, \
                                                                      adj_pval, vc_overall_conf, \
                                                                      vc_genera, vc_families, vc_orders, genus_conf
        except Exception as e:
            logger.error(e)
            exit(1)

    node_summary_df['Quality'] = node_summary_df['Quality'].apply(lambda x: round(x, 4))
    node_summary_df['Adj P-value'] = node_summary_df['Adj P-value'].apply(lambda x: round(x, 8))
    node_summary_df['Genus Confidence Score'] = node_summary_df['Genus Confidence Score'].apply(lambda x: round(x, 4))
    node_summary_df['Topology Confidence Score'] = node_summary_df['Topology Confidence Score'].apply(lambda x: round(x, 4))

    node_summary_df = node_summary_df.append(excluded)[node_summary_df.columns.tolist()]
    node_summary_df[['Order', 'Family', 'Genus']] = node_summary_df[['Order', 'Family', 'Genus']].fillna('Unassigned')

    node_summary_df.to_csv(os.path.join(folder, 'genome_by_genome_overview.csv'))
