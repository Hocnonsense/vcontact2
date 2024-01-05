"""
Modules are groups of protein families
"""
import logging
import subprocess
import os

import pandas as pd
import numpy as np
import scipy.sparse as sparse
import scipy.stats as stats

from .pcprofiles import PCProfiles

# import pcprofiles

logger = logging.getLogger(__name__)


class Modules(object):
    def __init__(
        self,
        pcp: PCProfiles,
        folder,
        inflation=5,
        threshold=10,
        shared_min=3,
        name=None,
        threads=0,
    ):
        """
        Args:
            profiles: PCProfiles object or a tuple:
                (
                    features (df),
                    contigs (df),
                    pc similarity-network (sp.matrix),
                    pc-profiles (sp.matrix)
                ).
            gc (df): (gc (df) and contigs (df))
            threshold: (int) Minimal sig value to take into account an edge.
            name: (str) A name to identify the object.
            folder (str): path where to save files to.
        """
        self.thres = threshold
        self.name = (
            "mod_sig{}_i{}".format(threshold, inflation) if name is None else name
        )
        self.folder = folder
        self.inflation = inflation
        self.shared_min = shared_min

        # profiles.contigs = pos, id (contigs), proteins
        # profiles.pcs = pos, id (pc), pc_id (same!), nb_proteins

        if isinstance(pcp, PCProfiles):
            self.pcs = pcp.pcs.copy()  # same
            self.contigs = pcp.contigs.copy()  # same
            network = pcp.ntw_modules
            self.matrix: sparse.spmatrix = pcp.matrix
        else:
            logging.debug("Reading input from tuple")
            self.pcs, self.contigs, network, self.matrix = (
                pcp[0].copy(),
                pcp[1].copy(),
                pcp[2],
                pcp[3],
            )

        # self.pcs = pos, id (pc), pc_id (same), nb_contigs
        # self.contigs = pos, id (contig), proteins (#)
        # self.matrix = (contig, PC) True/False

        # Filter the network according to the threshold:
        before = network.getnnz()
        self.network: sparse.csr_matrix = network.multiply(network >= self.thres)
        logger.debug(
            "Filtered {} edges according to the sig. threshold {}.".format(
                before - self.network.getnnz(), self.thres
            )
        )

        # Define the modules, runs MCL on PCs -> modules,
        # df w/ annotated_proteins, id, pos, position, size
        self.modules = self.define_modules(
            self.network, self.folder, self.inflation, threads
        )
        # Calls self.matrix, returns matrix w/ proportion of module's PCs in contig
        self.matrix_module = self.module_in_contigs()

    def __repr__(self):
        return (
            "Modules object {}, {} modules, (contigs, pc) : {}," "sig. threshold {} "
        ).format(self.name, len(self.modules), self.matrix.shape, self.thres)

    def define_modules(self, matrix, folder, I, threads=0):
        """
        Save the pc network in a file ready for MCL
        Run MCL
        Load clusters from the MCL results

        Args:
            matrix: (scipy.sparse matrix) network.
            folder: (str) folder.
            I: (float) Inflation for mcl.

        Returns:
            A dataframe containing:
                id: the id of the module.
                size: the number of protein clusters in the module.
                pos: the position of the module in the matrix.
                proteins: the number of proteins in the module.
                annotated_proteins: the number of annotated proteins in the module.

        Side-Effects:
            self.pcs: Add the column "module".

        Saved Files:
            name.ntwk: The pc similarity network.
            name_mcl_I.clusters: mcl results.
            name_mcl_I_modules.pandas: the module dataframe.
            name_mcl_I_pcs.pandas: the pc dataframe.
        """
        basename = "modules"
        fi_in = os.path.join(folder, "{}.ntwk".format(basename))
        fi_out = os.path.join(folder, "{}_mcl_{}.clusters".format(basename, I))
        fi_dataframe = os.path.join(
            folder, "{}_mcl_{}_modules.pandas".format(basename, I)
        )
        fi_feat = os.path.join(folder, "{}_mcl_{}_pcs.pandas".format(basename, I))

        # Save for MCL
        logger.info("Exporting the PC-network for MCL")
        if not os.path.exists(fi_in):
            names = self.pcs.set_index("pos")["pc_id"]
            with open(fi_in, "wt") as f:
                matrix = sparse.dok_matrix(matrix)
                for r, c in zip(*matrix.nonzero()):
                    f.write(
                        " ".join([str(x) for x in (names[r], names[c], matrix[r, c])])
                    )
                    f.write("\n")

            logger.debug(
                "Saving network in file {0} ({1} lines)".format(fi_in, matrix.getnnz())
            )
        else:
            logger.debug("Network file {} already exist.".format(fi_in))

        # Run MCL
        logger.info("Clustering the PC similarity-network")
        if not os.path.exists(fi_out):
            subprocess.call(
                f"mcl {fi_in} -o {fi_out} --abc -I {I} -te {threads}",
                shell=True,
            )
            logger.debug("MCL({}) results are saved in {}.".format(I, fi_out))
        else:
            logger.debug(f"MCL({I}) file already exist.")

        # Read MCL results
        logger.info("Loading the clustering results")
        if not os.path.exists(fi_dataframe) or not os.path.exists(fi_feat):
            with open(fi_out) as f:
                clusters = [
                    x
                    for x in (line.rstrip("\n").split("\t") for line in f)
                    if len(x) > 1
                ]  # Drop singletons
            nb_modules = len(clusters)
            formatter = "MD_{{:>0{}}}".format(int(round(np.log10(nb_modules)) + 1))
            module_names = [formatter.format(i) for i in range(nb_modules)]  # MD_XXXX
            module_size = [len(i) for i in clusters]
            pos = range(nb_modules)

            proteins = np.zeros(nb_modules)
            annotated_proteins = np.zeros(nb_modules)

            self.pcs = self.pcs.reset_index().set_index(
                "pc_id"
            )  # index = PC_XXXX, index, pos, pc_id, nb_contigs, module
            self.pcs["module"] = np.nan  # Setting NaN here allows for dropping later

            # print(self.contigs)  # pos, id, proteins
            # need the new gc.contigs, which = id, index, pos, proteins, pos_cluster, membership (T/F), pos_cluster_mbs

            for i, cluster in enumerate(
                clusters
            ):  # Clusters is list-of-list of PCs, with module "groupings"
                for n in cluster:  # For PC_xxxx in the PC group, i.e. in the module i
                    self.pcs.loc[n, "module"] = i
                    proteins[i] += self.pcs.loc[
                        n, "nb_proteins"
                    ]  # How many proteins in PC_XXXX?
                    annotated_proteins[i] += self.pcs.loc[n, "annotated"]

            self.pcs.dropna(subset=["module"], inplace=True)

            dataframe = pd.DataFrame(
                {
                    "id": module_names,
                    "size": module_size,
                    "pos": pos,
                    "proteins": proteins,
                    "annotated_proteins": annotated_proteins,
                }
            )

            dataframe.to_pickle(fi_dataframe)
            self.pcs.to_pickle(fi_feat)
            logger.debug(
                ("Saving {} modules containing {} " " protein clusters in {}.").format(
                    len(module_names), sum(module_size), fi_dataframe
                )
            )
        else:
            dataframe = pd.read_pickle(fi_dataframe)
            self.pcs = pd.read_pickle(fi_feat)
            logger.debug(
                "Read {} modules from {}.".format(len(dataframe), fi_dataframe)
            )
        return dataframe

    def module_in_contigs(self):
        """
        Compute the presence of modules in each contig

        Args:
            matrix (scipy.sparse):, M[contigs,pc] =  "pc in contig" (bool)
            contigs (pd.DataFrame):
            modules (pd.DataFrame): with columns : pos,

        Returns:
            matrix: S, S[contig, module] = proportions of module's pcs
                present in contig (in [0,1])
        """

        matrix: sparse.csc_matrix = self.matrix.tocsc()
        # Number of pcs of module m in each contig.
        N = sparse.lil_matrix((len(self.contigs), len(self.modules)))
        for m, data in self.pcs.reset_index().groupby("module"):
            pos = data["pos"].values
            N[:, int(m)] = matrix[:, pos].sum(1)

        # Number of pcs in each module
        # sizes = np.matrix([self.modules.size.values]*N.shape[0])
        sizes = np.matrix([self.modules["size"].values] * N.shape[0])
        # Compute
        S = sparse.csc_matrix(N.todense() / sizes)

        return S

    def link_modules_and_clusters(
        self,
        clusters: pd.DataFrame,
        contigs: pd.DataFrame,
        matrix_module: sparse.csc_matrix,
        thres=1.0,
        own_threshold=0.5,
    ):
        """
        Link the modules with the contigs clusters
        using the hypergeometric formula.

        Args:
            matrix_modules (scipy.sparse):
                bool(M[contig,module]>=own_threshold): module is present in contig
            modules (pd.DataFrame): information about the protein modules
            clusters (pd.DataFrame): information about the contig clusters
            contigs (pd.DataFrame): information about the contigs
            thres: (float) significativity threshold.
            own_threshold: minimal proportion of PCs to "own" the module.

        Returns:
            S (scipy.sparse) S[module,cluster] = sig.

        Formula:
            P(X>=c) ~ H(n,a,b)
                c: number of contigs in cluster c owning module m.
                n: number of contigs.
                a: number of contigs in the cluster c.
                b: number of contigs owning module m.
        """

        # Filtering
        non_filtered = matrix_module.getnnz()
        matrix_module = matrix_module >= own_threshold
        logger.info(
            "{} contigs-modules owning association, "
            "{} filtered (a contig must have {:.0%} of the PCs to own a module).".format(
                matrix_module.getnnz(),
                non_filtered - matrix_module.getnnz(),
                own_threshold,
            )
        )

        nb_contigs, nb_modules = matrix_module.shape
        nb_clusters = len(clusters)
        matrix_module = matrix_module.tocsc()
        logger.info(
            "Linking {} modules with {} contigs clusters...".format(
                nb_modules, nb_clusters
            )
        )

        # Phage in a given cluster
        a_values = clusters.sort_values(by="pos")["size"].values
        # Phage displaying a given module
        b_values = matrix_module.sum(0).A1
        # contig in cluster
        xy = (
            contigs.reset_index(drop=True)
            .loc[:, ["pos_cluster", "pos"]]
            .dropna(subset=["pos_cluster"])
            .sort_values(by="pos")
            .values
        )
        pa_matrix = sparse.coo_matrix(
            ([1] * len(xy), zip(*xy)), shape=(nb_clusters, nb_contigs)
        )
        # Phage in a given cluster displaying a given module
        c_values: sparse.csr_matrix = pa_matrix.dot(matrix_module)

        # Number of comparisons
        logT = np.log10(nb_clusters * nb_modules)

        S = sparse.lil_matrix((nb_clusters, nb_modules))
        for A, B in zip(*c_values.nonzero()):
            # choose(a, k) * choose(C - a, b - k) / choose(C, b)
            # sf(k) = survival function = 1 -cdf(k) = 1 - P(x<k) = P(x>k)
            # sf(k-1)= P(x>k-1) = P(x>=k)
            pval = stats.hypergeom.sf(
                c_values[A, B] - 1, nb_contigs, a_values[A], b_values[B]
            )
            sig = np.nan_to_num(-np.log10(pval) - logT)
            if sig > thres:
                S[A, B] = min(300, sig)

        logger.info(
            "Network done {0[0]} clusters, {0[1]} modules and {1} edges.".format(
                S.shape, S.getnnz()
            )
        )
        return S

    def link_modules_and_clusters_df(
        self,
        clusters: pd.DataFrame,
        contigs: pd.DataFrame,
        thres=1,
        own_threshold=0.7,
    ):
        """Returns the dataframe giving the association between clusters and modules

        Args:
            matrix_modules (scipy.sparse):
                bool(M[contig,module]>=own_threshold): module is present in contig
            modules (pd.DataFrame): information about the protein modules
            clusters (pd.DataFrame): information about the contig clusters
            contigs (pd.DataFrame): information about the contigs
            thres: (float) significativity threshold.
            matrix_module: ()
            own_threshold: minimal proportion of PCs to "own" the module.

        Returns:
            A dataframe containing the position of the module,
            of the cluster, the significativity of the association
            and merged with the lines from modules and clusters..
        """
        modules: pd.DataFrame = self.modules
        matrix_module: sparse.csc_matrix = self.matrix_module  # Module is in contig

        matrix = (
            self.link_modules_and_clusters(
                clusters, contigs, matrix_module, thres, own_threshold
            )
            .todok()
            .items()
        )
        data = (
            pd.DataFrame(
                {
                    "pos_module": [x[0][1] for x in matrix],
                    "pos_cluster": [x[0][0] for x in matrix],
                    "sig": [x[1] for x in matrix],
                }
            )
            .merge(modules, left_on="pos_module", right_on="pos", how="left")
            .drop(columns="pos")
            .merge(
                clusters,
                left_on="pos_cluster",
                right_on="pos",
                how="left",
                suffixes=["_module", "_cluster"],
            )
            .drop(columns="pos")
        )
        return data
