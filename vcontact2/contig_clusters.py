"""Contig clusters : An object to work on the similarity network and
make the affiliations"""

import os
import subprocess
import _pickle as pickle
import logging
from io import StringIO

from pprint import pprint

import numpy as np
import scipy.sparse as sparse
import pandas as pd

from . import pcprofiles
from . import matrices
from . import ml_functions
from . import associations

logger = logging.getLogger(__name__)


class ContigCluster(object):
    """ Deal with the clusters of the contig similarity network

    Attributes:
            pcs: (pandas df) Protein clusters with pos, pc_id, size, annotated, keys, nb_proteins
            contigs: (pandas df) Contigs & Reference genomes with contig_id, index, pos, proteins, origin, order,
            family, genus, pos_cluster, membership, pos_cluster_mbship
            network: (sparse matrix) Contig similarity network
            taxonomy: (pandas df) Taxonomic class
            clusters: (pandas df) Contig clusters
            mcl_results: (list of list) mcl_result[cluster][prot]
    """
    
    def __init__(self, profiles, folder, cluster_one, one_args, inflation=2, threshold=None, name=None,
                 membership_simple=False, mode='MCL'):
        """
        Init the object with a pc-profile object and perform the clustering

        Args:
            profiles: PCprofiles object or a tuple
                 (pcs (df), contigs (df), network (sparse matrix)).
            inflation (float): inflation for mcl.
            folder (str): path where to save files
            threshold (float): minimal significativity. 
            name (str): A name to identify the object.
            membership_simple (bool): if false use non boolean membership. 
        """
        self.mode = mode
        if self.mode not in ['ClusterONE', 'MCL']:
            logger.error('Unable to identify which clustering mode is enabled.')

        if mode == 'MCL':
            self.name = "cc_sig{}_mcl{}".format(threshold, inflation) if name is None else name
        if mode == 'ClusterONE':
            self.name = "c1" if name is None else name

        self.inflation = inflation
        self.thres = threshold
        self.folder = folder
        self.one_opts = one_args
        self.cluster_one = cluster_one
        self.df = None  # Keep copy of contig dataframe to identify singletons, outliers and overlaps later

        if isinstance(profiles, pcprofiles.PCProfiles):
            self.pcs = profiles.pcs.copy()
            self.contigs = profiles.contigs.copy()
            self.network = profiles.ntw
        else:
            logging.debug("Reading input from tuple")
            self.pcs, self.contigs, self.network = profiles[0].copy(), profiles[1].copy(), profiles.network

        if threshold is not None:
            before = self.network.getnnz()
            self.network = self.network.multiply(self.network >= self.thres)
            delta = before-self.network.getnnz()
            if delta:
                logger.debug("Filtered {} edges according to the sig. threshold {}.".format(delta, self.thres))

        # Columns != contig_id, proteins, size, origin and pos
        self.levels = frozenset(self.contigs.columns) - frozenset(["contig_id", "proteins", "size", "origin", "pos"])
        logger.debug("{} taxonomic levels detected: {}".format(len(self.levels),", ".join(self.levels)))

        # Set
        self.taxonomy = self.extract_taxonomy(self.contigs, levels=self.levels)

        # Export to MCL, run MCL, return clusters
        if mode == 'ClusterONE':
            self.clusters, self.cluster_results = self.one_cluster(os.path.join(self.folder, self.name),
                                                                   self.cluster_one, self.one_opts)

        if mode == 'MCL':
            self.clusters, self.cluster_results = self.mcl_cluster(os.path.join(self.folder, self.name),
                                                                   self.inflation)

        self.matrix = {}
        logger.info("Computing membership matrix...")
        self.permissive = not membership_simple
        if membership_simple:
            self.matrix["B"] = matrices.bool_membership(self.contigs)
        else:
            self.matrix["B"] = matrices.membership(self.cluster_results, self.network, self.contigs)

        self.contigs = associations.contig_cluster(self.contigs, self.matrix["B"])

    def __repr__(self):
        return "GenomeCluster object {}, {} contigs and {} clusters".format(
            self.name, len(self.contigs), len(self.clusters))

    #--------------------------------------------------------------------------#
    # PART 1 : IMPORT, EXTRACTION, CLUSTERING
    #--------------------------------------------------------------------------#
    def extract_taxonomy(self, contigs=None, levels=None):
        """ Build the taxonomy dataframe.

        Args:
            contigs (pandas.DataFrame): with column "level" and "pos".
            levels (list): column name in contigs

        Returns:
            dict: A dictionary of pandas.DataFrame, one key by taxonomic level.
        """
        contigs = self.contigs if contigs is None else contigs
        levels = self.levels if levels is None else levels
        tax = {}

        for t in levels:  # order, family, genus, ...
            # tax[t] = pd.DataFrame(contigs.groupby(t).pos.size(), columns=["references"])
            tax[t] = contigs.groupby(t).pos.count().to_frame()
            tax[t].columns = ['references']
            tax[t].index.name = "name"
            tax[t].reset_index(inplace=True)

            tax[t]["pos"] = tax[t].index

        return tax

    def mcl_cluster(self, basename, I, force=False):
        """Export the matrix, Run MCL and load the results

        Args:
            basename: (str) Path for the exported files
            I: (float) inflation for mcl
            force: (bool) overwrite existing file

        Returns:
            See self.load_clusters.

        Side-Effects:
           Save basename.ntw the network file
           Save basename.clusters the clustering results
           self.contig: add "cluster" column
        """

        fi_ntw = basename + ".ntw"
        fi_clusters = basename + ".clusters"

        # Export for MCL
        logger.info("Exporting for MCL")
        if not os.path.exists(fi_ntw) or force:
            self.to_clusterer(self.network, fi_ntw)
        else:
            logger.debug("Network file already exist.")

        # MCL
        logger.info("Clustering the pc similarity-network")

        if not os.path.exists(fi_clusters or force):
            subprocess.call("mcl {0} -o {1} --abc -I {2}".format(fi_ntw, fi_clusters, I), shell=True)
            logger.debug("MCL({}) results are saved in {}.".format(I, fi_clusters))
        else:
            logger.debug("MCL({}) file already exist.".format(I, fi_clusters))

        # Load clusters
        return self.load_mcl_clusters(fi_clusters)

    def one_cluster(self, basename, cluster_one, options, force=False):
        """Export the matrix, Run ClusterONE and load the results

        Args:
            basename: (str) Path for the exported files
            options: (float) inflation for mcl
            force: (bool) overwrite existing file

        Returns:
            See self.load_clusters.

        Side-Effects:
           Save basename.ntw the network file
           Save basename.clusters the clustering results
           self.contig: add "cluster" column
        """

        fi_ntw = basename + ".ntw"
        fi_clusters = basename + ".clusters"

        # Export for ClusterONE
        logger.info("Exporting for ClusterONE")
        if not os.path.exists(fi_ntw) or force:
            self.to_clusterer(self.network, fi_ntw)
        else:
            logger.debug("Network file already exist.")

        # ClusterONE
        logger.info("Clustering the PC Similarity-Network using ClusterONE")

        if not os.path.exists(fi_clusters) or force:

            # Disable --fluff as it's not in published algorithm or used in published manuscript
            if 'jar' in cluster_one:
                cluster_one_cmd = 'java -jar {} {} --input-format edge_list --output-format csv'.format(cluster_one,
                                                                                                        fi_ntw)
            else:
                cluster_one_cmd = '{} {} --input-format edge_list --output-format csv'.format(cluster_one, fi_ntw)

            for opt, val in options.items():
                cluster_one_cmd += ' {} {}'.format(opt, val)

            cluster_one_cmd += ' > {}'.format(fi_clusters)

            logger.info('Running clusterONE: {}'.format(cluster_one_cmd))
            subprocess.call(cluster_one_cmd, shell=True)
            logger.debug("ClusterONE results are being saved to {}.".format(fi_clusters))

        else:
            logger.debug("ClusterONE file {} already exist.".format(fi_clusters))

        # Load clusters
        return self.load_one_clusters(fi_clusters)

    def to_clusterer(self, matrix, fi, names=None):
        """Save a network in a file ready for MCL and/or ClusterONE

        Args:
            matrix (scipy.sparse_matrix): network.
            fi (str): filename .
            names (pandas.dataframe): with the columns
                "pos":  (int) is the position in the matrix.
                "id": (str) column contain the id of the node.
                If None, self.contigs is used.

        Returns:
            str: filename
        """

        names = self.contigs if names is None else names
        names = names.set_index("pos").contig_id
        with open(fi, "wt") as f:
            matrix = sparse.dok_matrix(matrix)
            for r, c in zip(*matrix.nonzero()):
                f.write(" ".join([str(x) for x in (names[r], names[c], matrix[r, c])]))
                f.write("\n")

        logger.debug("Saving network in file {0} ({1} lines).".format(fi, matrix.getnnz()))
        return fi

    def load_mcl_clusters(self, mcl_fi):
        """ Load clusters from the mcl results

        Args:
            mcl_fi (str): path to the MCL result file.

        Returns:
            df (pandas.DataFrame): give for each contig cluster
                its name, size and position in the matrix.

        Side-Effect:
            Modify self.contig to add the column "pos_cluster"
            giving the pos of the cluster it belongs to.

        The file fi was probably generated using :
        "mcl <file>.ntw --abc -I 2 -o <file>.clusters".
        """

        # Read the files
        with open(mcl_fi) as f:
            c = [line.rstrip("\n").split("\t") for line in f]
        c = [x for x in c if len(c) > 1]
        nb_clusters = len(c)
        formatter = "CC_{{:>0{}}}".format(int(round(np.log10(nb_clusters)) + 1))
        name = [formatter.format(str(i)) for i in range(nb_clusters)]
        size = [len(i) for i in c]
        pos = range(nb_clusters)

        logger.info("{} clusters loaded (singletons and non-connected nodes are dropped).".format(len(c)))

        # Update self.contigs (To refactor)
        self.contigs.reset_index(inplace=True)
        self.contigs.set_index("contig_id", inplace=True)
        self.contigs["pos_cluster"] = np.nan

        for i, cluster in enumerate(c):
            for n in cluster:
                self.contigs.loc[n, "pos_cluster"] = i

        self.contigs.reset_index(inplace=True)  # self.contigs = id, index, pos, proteins, pos_cluster

        return pd.DataFrame({"id": name, "size": size, "pos": pos}), c

    def load_one_clusters(self, one_fn):
        """ Load clusters from ClusterONE results

        Args:
            one_fn (buffer): File buffer to the ClusterONE result file.

        Returns:
            df (pandas.DataFrame): give for each contig cluster
                its name, size and position in the matrix.

        Side-Effect:
            Modify self.contig to add the column "pos_cluster"
            giving the pos of the cluster it belongs to.

        The file fi was probably generated using :
        "mcl <file>.ntw --abc -I 2 -o <file>.clusters".
        """

        # Read file buffer
        # fh_clusters = StringIO(one_fh.decode('ascii', errors='ignore'))
        clusters_df = pd.read_csv(one_fn, header=0)
        clusters_df['Cluster'] = clusters_df['Cluster'].astype(int) - 1
        clusters_df['Cluster'] = clusters_df['Cluster'].astype(str)

        c = [line.rstrip().split() for line in clusters_df['Members']]
        c = [x for x in c if len(c) > 1]  # ClusterONE won't export singletons, they'll be outliers
        nb_clusters = len(c)
        formatter = "CC_{{:>0{}}}".format(int(round(np.log10(nb_clusters)) + 1))
        name = [formatter.format(str(i)) for i in range(nb_clusters)]
        size = [len(i) for i in c]
        pos = range(nb_clusters)

        logger.info("{} clusters loaded (singletons and non-connected nodes are dropped).".format(len(c)))

        # Update self.contigs (To refactor)
        self.contigs.reset_index(inplace=True)
        self.contigs.set_index("contig_id", inplace=True)
        self.contigs["pos_cluster"] = np.nan
        self.contigs["pos_clusters"] = np.nan

        for i, cluster in enumerate(c):
            for n in cluster:
                if pd.isnull(self.contigs.loc[n, "pos_cluster"]):  # If never seen before
                    self.contigs.loc[n, "pos_cluster"] = str(i)
                elif pd.isnull(self.contigs.loc[n, "pos_clusters"]):  # If has been seen, check if others have been seen
                    self.contigs.loc[n, "pos_clusters"] = ';'.join([str(self.contigs.loc[n, "pos_cluster"]), str(i)])
                else:  # If has been seen, along with others (NOTE the (s) at suffix)
                    self.contigs.loc[n, "pos_clusters"] = ';'.join([str(self.contigs.loc[n, "pos_clusters"]), str(i)])

        self.df = self.contigs.copy()

        # Remove overlapping contigs
        self.contigs.loc[pd.notnull(self.contigs["pos_clusters"]), "pos_cluster"] = np.nan

        # Side-effect of casting multiple pos as string (above), still need to process non-overlapping members
        self.contigs.loc[self.contigs['pos_cluster'].notnull(), 'pos_cluster'] = self.contigs.loc[self.contigs['pos_cluster'].notnull(), 'pos_cluster'].astype(int)

        self.contigs.reset_index(inplace=True)

        return pd.DataFrame({"id": name, "size": size, "pos": pos}), c

    #--------------------------------------------------------------------------#
    # PART 2: Affiliations
    #--------------------------------------------------------------------------#

    def total_affiliation(self, levels=None):
        """Routine of the analysis using all the dataset.

        Args:
            levels (tuple): Taxonomic levels to consider.

        Returns:
            dataframe: Classification metrics. One line by taxonomic level.

        Warning:
            This function directly modifies the object attributes.
        """

        results = []
        levels = self.levels if levels is None else levels

        for level in levels:
            logger.info("Affiliation at the {} level...".format(level))
            self.matrix[level] = {}
            # Taxonomic reference matrix
            self.matrix[level]["K"] = matrices.reference_membership(level, self.contigs, self.taxonomy[level])
            # (2, 0) True

            # recall, precision and F-measure matrix
            (self.matrix[level]["Q"],
             self.matrix[level]["R"],
             self.matrix[level]["P"],
             self.matrix[level]["F"]) = matrices.correspondence(self.matrix[level]["K"],  # sparse matrix, (#,#) bool
                                                                self.matrix["B"])  # 0, 1 matrix

            # Associate clusters and taxonomic classes.
            self.clusters, self.taxonomy[level] = associations.cluster_taxonomy(self.clusters,
                                                                                self.taxonomy[level],
                                                                                level,
                                                                                self.matrix[level]["P"],
                                                                                self.matrix[level]["R"])
            # Associate the contigs with the classes
            # Now contains contigs + predicted_*
            self.contigs = associations.contig_taxonomy(self.contigs,
                                                        self.taxonomy[level],
                                                        self.clusters,
                                                        level)

            # Compute classification metrics.
            logger.info("Computing the classification metrics")
            results.append(ml_functions.classification_metrics(self.contigs.loc[:, [level, "predicted_"+level]],
                                                               ref_col=level,
                                                               pred_col="predicted_"+level))
        return pd.DataFrame(results, levels)

    def cross_validation_affiliation(self, level="family", folds=10):
        """Cross validation affiliation.

        Cut the dataset (where a reference taxonomy exist) into <fold> equal parts
        and use alternativly <fold-2> of them to do the affiliation. Compute the
        classification metrics on the learning set and on the two remaining sets:
        cross validation (used to select the model) and test (used to have an 
        unbiased estimate of the error of the classification.)

        Note: 
            The splitting is stratified to keep the relative taxonomic classes in 
            the same proportions in each fold. 

        Args:
            level (str): Taxonomic level to consider.
            folds (int): number of folds for the cross-validation.

        Returns:
            dict: Dict of dataframes, one entry by set (learning,cv,test). 
                In the dataframes are the classification metrics, one row by
                selected cv-set. 
        """
        results = {"train_set": [],
                   "cv_set": [],
                   "test_set": []}
        conditions = {}
        contigs = self.contigs
        taxonomy = self.taxonomy[level]
        clusters = self.clusters
        
        # Stratified split according to the taxonomic level. 
        contigs = ml_functions.split_dataset(contigs,level,folds)
        logger.info("{} folds cross-validation".format(folds))
        all_sets = frozenset(range(folds))

        for cv_set, test_set in zip(range(folds), range(1, folds) + [0]):
            logger.info("Cross-validation fold {:2} ({:.0%})".format(cv_set, cv_set/float(folds)))

            # Filtering conditions: 
            train_set = list(all_sets - frozenset([cv_set,test_set]))
            conditions["cv_set"] = "cvset_{}=={}".format(level, cv_set)
            conditions["test_set"] = "cvset_{}=={}".format(level, test_set)
            conditions["train_set"] = "cvset_{} in {}".format(level, train_set)

            # Affiliation: 
            K = matrices.reference_membership(level, contigs, taxonomy, conditions["train_set"])
            _,R,P,_ = matrices.correspondence(K, self.matrix["B"])
            clusters, taxonomy = associations.cluster_taxonomy(clusters, taxonomy,
                                                               level, P, R)
            contigs = associations.contig_taxonomy(contigs, taxonomy, clusters, level)

            # Computing metrics:
            for set_ in conditions.keys():
                df = contigs.query(conditions[set_]).loc[:, [level, "predicted_"+level]]
                results[set_].append(ml_functions.classification_metrics(df,
                                                                         ref_col=level,
                                                                         pred_col="predicted_"+level))

            # Cleaning the dataset for next step:
            contigs = contigs.drop("predicted_"+level,1)
            taxonomy = taxonomy.drop(["pos_cluster","recall"],1)
            clusters = clusters.drop(["pos_"+level,"precision_"+level],1)
        for set_ in results.keys():
            results[set_] = pd.DataFrame(results[set_],[range(folds)])
        return results

    def learning_curve_affiliation(self, level="family", folds=10):
        """Learning curve affiliation.

        Cut the dataset (where a reference taxonomy exist) into <fold> equal parts
        and increasingly one to <fold-1> of them to do the affiliation. Compute the
        classification metrics on the learning set and on a randomly picked
        cross validation set in the non used subset. 

        Note: 
            The splitting is stratified to keep the relative taxonomic classes in 
            the same proportions in each fold. 

        Args:
            levels (str): Taxonomic level to consider.
            folds (int): number of folds for the cross-validation.

        Returns:
            dict: Dict of dict of dataframes, one entrie by level then one entry
                by set (learning,cv). In the dataframes are the classification
                metrics), one row by learning set size . 
        """
    
        results = {"train_set": [],
                   "cv_set": []}
        conditions = {}
        contigs = self.contigs
        taxonomy = self.taxonomy[level]
        clusters = self.clusters
        
        # Stratified split according to the taxonomic level. 
        contigs = ml_functions.split_dataset(contigs, level, folds)
        logger.info("{} folds cross-validation".format(folds))

        for cv_set in range(folds):
            logger.info("Cross-validation set {:2} ({:.0%})".format(cv_set, cv_set/float(folds)))

            remaining_sets = range(folds)
            remaining_sets.pop(cv_set)
            train_sets = [[remaining_sets[y] for y in range(x)] for x in range(1,folds)]
            for train_set in train_sets:
                logger.info("Training set of size {:2}".format(len(train_set)))

                # Filtering conditions: 
                conditions["cv_set"] = "cvset_{}=={}".format(level, cv_set)
                conditions["train_set"] = "cvset_{} in {}".format(level, train_set)

                # Affiliation: 
                K = matrices.reference_membership(level, contigs,
                                                  taxonomy, conditions["train_set"])
                _,R,P,_ =  matrices.correspondence(K, self.matrix["B"])
                clusters, taxonomy = associations.cluster_taxonomy(clusters, taxonomy,
                                                                   level, P, R)
                contigs =  associations.contig_taxonomy(contigs, taxonomy, clusters,
                                                        level)

                # Computing metrics:
                for set_ in results.keys():
                    df = contigs.query(conditions[set_]).loc[:,[level,"predicted_"+level]]
                    metrics = ml_functions.classification_metrics(df,
                                                                  ref_col=level,
                                                                  pred_col="predicted_"+level)
                    metrics["train_size"] = len(train_set)
                    metrics["cv_set"] = cv_set
                    results[set_].append(metrics)

                # Cleaning the dataset for next step:
                contigs = contigs.drop("predicted_"+level, 1)
                taxonomy = taxonomy.drop(["pos_cluster","recall"], 1)
                clusters = clusters.drop(["pos_"+level,"precision_"+level], 1)
        for set_ in results.keys():
            results[set_] = pd.DataFrame(results[set_])
        return results

    #--------------------------------------------------------------------------#
    # PART 4: Pickle-save
    #--------------------------------------------------------------------------#

    def to_pickle(self,path=None):
        """ Pickle (serialize) object to file path."""
        path = self.name + ".pkle" if path is None else path
        with open(path, 'wb') as f:
            pickle.dump(self, f)


def read_pickle(path):
    """Read pickled object in file path."""
    with open(path, 'rb') as fh:
        return pickle.load(fh)
