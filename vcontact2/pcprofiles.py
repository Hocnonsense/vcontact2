""" Object containing the pc-profiles of the contigs and able to
    compute the similarity network on the contigs and the pcs. """
import logging
import sys
import pandas as pd
import numpy as np
import scipy.stats as stats
import scipy.sparse as sparse
import networkx
import _pickle as pickle

import multiprocessing as mp
from itertools import zip_longest, chain
import math

logger = logging.getLogger(__name__)

# Divide by zero for singletons
np.seterr(divide='ignore')


class PCProfiles(object):
    """
    Protein Cluster presence/absence matrix, or PCprofiles.
    Used to build the similarity networks.

    Attributes:
        contigs (pandas.DataFrame):
        pcs (pandas.DataFrame):
        name (str): 
        profiles (sparse.matrix):
        singletons (sparse.matrix): 
        contig_ntw (sparse.matrix):
        modules_ntw (sparse.matrix):
    """

    def __init__(self, contigs, pcs, profiles, threads, name=None, sig=1.0, max_sig=300, sig_mod=1.0, mod_shared_min=3):
        """
        Args:
            contigs (dataframe): Contig info, required fields are pos, id.
            pcs (dataframe): Protein clusters info, required fields are pos, id.
            profiles (dict): Required field are matrix (the contig X PC profiles
                matrix) and singletons (the contigs X1 matrix giving the number
                of singletons by contig.)`
            threads (int): Number of CPUs
            sig (float): Sig. threshold in the contig similarity network.
            max_sig (int): Maximum sig ceiling
            sig_mod (float): Sig. threshold in the pc similarity network.
            mod_shared_min (float): Minimal number of contigs a pc must appear into
                to be taken into account in the modules computing.
            name (str): name the object (useful in interactive mode)
        """

        if name is None:
            self.name = "PCprofiles"

        self.name = "PCprofiles"
        self.threads = threads

        # Get the data
        self.contigs = contigs   # pos, id, proteins
        self.pcs = pcs   # pos, id, pc_id, nb_contigs
        self.matrix, self.singletons = profiles["matrix"], profiles["singletons"]

        # Store parameters
        self.sig = sig
        self.sig_mod = sig_mod
        self.mod_shared_min = mod_shared_min

        # Compute networks
        self.ntw = self.network(self.matrix, self.singletons, thres=sig, max_sig=max_sig, threads=self.threads)
        self.ntw_modules = self.network_modules(self.matrix, thres=sig_mod, mod_shared_min=mod_shared_min,
                                                threads=self.threads)

    def __repr__(self):
        return ("PC-profiles object {0}\n{1:10} contigs by {2} shared protein clusters and {3} singletons.\n"
                "Contig sig threshold: {4:10}\nModule sig threshold: {5:10}\nModule shared min {6:10}".format(
            self.name, self.matrix.shape[0], self.matrix.shape[1], self.singletons.sum(), self.sig, self.sig_mod,
            self.mod_shared_min))

    def network(self, matrix, singletons, thres=1, max_sig=1000, threads=1):
        """
        Compute the hypergeometric-similarity contig network.

        Args:
            matrix (scipy.sparse)x: contigs x protein clusters :
                M(c,p) == True <-> PC p is in Contig c.
            thres (float): Minimal significativity to store an edge value.
            max_sig (int): Maximum significance score
        
        Return
            scipy.sparse: S symmetric lil matrix, contigs x contigs.
          S(c,c) = sig(link)
        """

        # There are 
        contigs, pcs = matrix.shape
        pcs += singletons.sum()
        
        # Number of comparisons
        T = 0.5 * contigs * (contigs - 1)
        logT = np.log10(T)

        # Number of protein clusters in each contig
        # = # shared pcs + #singletons
        number_of_pc = matrix.sum(1) + singletons
        number_of_pc = number_of_pc.A1  # Transform into a flat array

        # Number of common protein clusters between two contigs, tuple + commons
        commons_pc = matrix.dot(sparse.csr_matrix(matrix.transpose(), dtype=int))

        S = sparse.lil_matrix((contigs, contigs))

        total_c = float(commons_pc.getnnz())

        if threads == 1:
            i = 0  # Display
            for A, B in zip(*commons_pc.nonzero()):  # For A & B sharing contigs
                if A != B:
                    # choose(a, k) * choose(C - a, b - k) / choose(C, b)
                    # sf(k) = survival function = 1 -cdf(k) = 1 - P(x<k) = P(x>k)
                    # sf(k-1)= P(x>k-1) = P(x>=k)
                    # It is symmetric but I put the smallest before to avoid numerical bias.
                    a, b = sorted([number_of_pc[A], number_of_pc[B]])
                    pval = stats.hypergeom.sf(commons_pc[A, B] - 1, pcs, a, b)
                    sig = min(max_sig, np.nan_to_num(-np.log10(pval) - logT))
                    if sig > thres:
                        S[min(A, B), max(A, B)] = sig
                    # Display
                    i += 1
                    if i % 1000 == 0:
                        sys.stdout.write(".")
                    if i % 10000 == 0:
                        sys.stdout.write("{:6.2%} {}/{}\n".format(i / total_c, i, total_c))
        else:
            pool = mp.Pool(processes=threads)

            results = []
            # https://stackoverflow.com/questions/26104512/how-to-obtain-the-results-from-a-pool-of-threads-in-python
            commons_contigs = [pair for pair in zip(*commons_pc.nonzero())]

            for partition in self.grouper(commons_contigs, math.ceil(len(commons_contigs) / threads)):
                non_null_part = [item for item in partition if item]  # Remove None from grouper
                results.append(
                    pool.apply_async(self.hypergeom, args=(number_of_pc, commons_pc, pcs, logT, thres, max_sig, non_null_part)))

            pool.close()
            pool.join()

            final_results = list(chain.from_iterable([r.get() for r in results]))

            for (A, B, sig) in final_results:
                if sig > thres:
                    S[min(A, B), max(A, B)] = sig

        S += S.T  # Symmetry
        S = S.tocsr()
        if len(S.data) != 0:
            logger.debug("Hypergeometric contig-similarity network:\n {0:10} contigs,\n {1:10} edges (min:{2:.2}"
                         "max: {3:.2}, threshold was {4})".format(contigs, S.getnnz(), S.data.min(), S.data.max(), thres))
        else:
            raise ValueError("No edge in the similarity network !") 
        return S

    def network_modules(self, matrix=None, thres=1, mod_shared_min=3, threads=1):
        """
        Compute the hypergeometric-similarity pc network.

        Warning:
            Use only the PCs that are present in 3 contigs or more.

        Args: 
            matrix (scipy.sparse): contigs x protein clusters :
                M(c,p) == True <-> PC p is in Contig c.
            thres (float): Minimal significativity to store an edge value.
            mod_shared_min (float): Minimal number of contigs a pc must appear into
                to be taken into account in the modules computing.

        Returns:
            scipy.sparse: Symmetric lil_matrix, PCs x PCs
                S(c,c) = sig(link)
        """

        matrix = self.matrix if matrix is None else matrix
        contigs = matrix.shape[0]

        # Number of contig in which a given PC is found
        number_of_contigs = np.squeeze(np.asarray(np.transpose(matrix.sum(0))))

        # We only keep the pcs that are presents in more than "mod_shared_min" contigs
        pos_pcs_in_modules = [i for i, x in enumerate(number_of_contigs) if x >= mod_shared_min]
        pcs_in_modules = len(pos_pcs_in_modules)
        logger.debug("{} PCs present in strictly more than {} contigs".format(pcs_in_modules, mod_shared_min))

        # Filtering the matrix to only include pcs with >= 3 contigs
        matrix = matrix[:, pos_pcs_in_modules]

        # Number of comparisons
        T = 0.5 * pcs_in_modules * (pcs_in_modules - 1)
        logT = np.log10(T)

        # Number of common contigs between two pcs
        commons_contigs = sparse.csr_matrix(matrix, dtype=int).transpose().dot(matrix)

        # Initiate empty matrix
        S = sparse.lil_matrix((pcs_in_modules, pcs_in_modules))

        total_c = float(commons_contigs.getnnz())

        if threads == 1:

            i = 0
            for A, B in zip(*commons_contigs.nonzero()):
                if A != B:
                    # choose(a, k) * choose(C - a, b - k) / choose(C, b)
                    # sf(k) = survival function = 1 -cdf(k) = 1 - P(x<k) = P(x>k)
                    # sf(k-1)= P(x>k-1) = P(x>=k)
                    # It is symmetric but I put the smallest before to avoid numerical biases.

                    a, b = sorted([number_of_contigs[A], number_of_contigs[B]])
                    pval = stats.hypergeom.sf(commons_contigs[A, B]-1, contigs, a, b)
                    sig = min(300, np.nan_to_num(-np.log10(pval)-logT))

                    if sig > thres:
                        S[min(A, B), max(A, B)] = sig

                    i += 1
                    if i % 1000 == 0:
                        sys.stdout.write(".")
                    if i % 10000 == 0:
                        sys.stdout.write("{:6.2%} {}/{}\n".format(i/total_c, i, total_c))
        else:
            pool = mp.Pool(processes=threads)

            results = []

            commons_modules = [pair for pair in zip(*commons_contigs.nonzero())]

            for partition in self.grouper(commons_modules, math.ceil(commons_contigs.shape[0] / threads)):
                non_null_part = [item for item in partition if item]
                results.append(pool.apply_async(self.hypergeom,args=(number_of_contigs, commons_contigs, contigs, logT, thres, 300, non_null_part)))

            pool.close()
            pool.join()

            final_results = list(chain.from_iterable([r.get() for r in results]))

            for (A, B, sig) in final_results:
                if sig > thres:
                    S[min(A, B), max(A, B)] = sig

        logger.debug("Hypergeometric PCs-similarity network : {0} pcs, {1} edges".format(pcs_in_modules, S.getnnz()))
        S += S.T  # Symmetry

        return S

    def nodes_properties(self, matrix):
        """ Compute several node specific statistics.
        
        Args:
            matrix: (sparse.matrix) contig-similarity network
        
        Side-effect:
            self.contigs: Add the following columns:
                degree: Number of edge to the node. 
                clustering_coefficient: Proportion of existing edges
                    over possible edges in the neighborhood of the node
                betweeness_centrality: sum of the fraction of all-pairs 
                    shortest path that pass trhough the node.
        """

        D = networkx.from_scipy_sparse_matrix(matrix)
        
        bc = pd.Series(networkx.betweenness_centrality(D), name="betweeness_centrality")
        degr = pd.Series(networkx.degree(D), name="degree")
        clcoef = pd.Series(networkx.clustering(D), name="clustering_coef")
                
        df = pd.concat([bc, degr, clcoef], axis=1)
        self.contigs = pd.merge(self.contigs, df, left_on="pos", right_index=True)

    def hypergeom(self, num_pcs, common_pcs, total_pcs, total_comparisons, thres, max_thres, shard):
        """
        Multithreader helper function

        Args:
            num_pcs (numpy.ndarray):
            common_pcs (scipy.sparse.csr.csr_matrix): Number of common protein clusters between two contigs
            total_pcs (int): Total PCs
            total_comparisons (float): Log transform of total comparisons between contigs x contigs
            thres (float): Minimum significance score to be retained
            max_thres (float): Maximum significance score allowed
            shard (list): List of contig x contig tuples of contig positions to be evaluated

        Returns:
            list: contig contig significance scores passing threshold
        """

        results = []

        for (A, B) in shard:

            if A != B:
                a, b = sorted([num_pcs[A], num_pcs[B]])
                pval = stats.hypergeom.sf(common_pcs[A, B] - 1, total_pcs, a, b)
                sig = min(max_thres, np.nan_to_num(-np.log10(pval) - total_comparisons))

                if sig > thres:
                    results.append([A, B, sig])

        return results

    def grouper(self, iterable, n, fillvalue=None):
        "Collect data into fixed-length chunks or blocks"
        # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
        args = [iter(iterable)] * n
        return zip_longest(*args, fillvalue=fillvalue)

    def to_pickle(self, path=None):
        """ Pickle (serialize) object to file path."""
        path = self.name + ".pkle" if path is None else path
        with open(path, 'wb') as f:
            pickle.dump(self, f)


def build_pc_matrices(profiles, contigs, pcs):
    """
    Build the pc profiles matrices (shared & singletons) from dataframes.

    Args:
        profiles (dataframe): required fields are contig_id and pc_id. # pos, contig_id, pc_id
        contigs (dataframe): contigs info, required field are proteins, pos and id. # pos, contig_id, proteins
        pcs (dataframe): pcs info, required field are pos and id.  # pos, id, size, annotated

    Returns:
        (tuple of sparse matrix): Shared PCs and singletons matrix.
    """

    pc_by_cont = profiles.groupby("contig_id").count().pc_id
    pc_by_cont = pc_by_cont.reset_index()
    pc_by_cont = pd.merge(contigs.sort_values("pos").loc[:, ["pos", "contig_id", "proteins"]], pc_by_cont, how="left",
                          left_on="contig_id", right_on="contig_id").fillna(0)
    singletons = (pc_by_cont.proteins - pc_by_cont.pc_id).values
    singletons = sparse.lil_matrix(singletons).transpose()

    # Matrix
    profiles.index.name = "pos"
    profiles.reset_index(inplace=True)

    # pc_id or contig?
    profiles = pd.merge(profiles, pcs.loc[:, ["pc_id", "pos"]], left_on="pc_id", right_on="pc_id", how="inner",
                            suffixes=["", "_pc"])  # pos, contig_id, pc_id, id (pc), pos_pc

    profiles = pd.merge(profiles, contigs.loc[:, ["contig_id", "pos"]], left_on="contig_id", right_on="contig_id", how="inner",
                            suffixes=["", "_contig"])

    profiles = profiles.loc[:, ["pos_contig", "pos_pc"]]

    matrix = sparse.coo_matrix(([1]*len(profiles), (zip(*profiles.values))), shape=(len(contigs), len(pcs)),
                               dtype="bool")

    return matrix.tocsr(), singletons.tocsr()


def read_pickle(path):
    """Read pickled object in file path."""
    with open(path, 'rb') as fh:
        return pickle.load(fh)