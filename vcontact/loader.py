""" Loading tools """
import pandas as pd
import logging
import scipy.sparse as sparse
import os
import _pickle as pickle

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def csv(contigs_fi, pcs_fi, pcprofiles_fi, folder, name, force=False):
    """
    Import csv file and build the info table and pc-profiles matrix.
    Save everything into h5 and pickle files. 

    Args:
        contigs_fi (str): path to the csv file containing the contigs
        pcs_fi (str): path to the csv file containing the pcs
        pcprofiles_fi (str): path to the csv file containing the pc_profiles
        folder (str): output folder path 
        name (str): experiment name 
        force (str): overwrite existing files
    """
    store = pd.HDFStore(os.path.join(folder, '{}.h5'.format(name)))

    if "contigs" not in store or force:
        if contigs_fi:
            contigs = pd.read_csv(contigs_fi, sep=',', low_memory=False)
            contigs['contig_id'] = contigs['contig_id'].str.replace(' ', '~')
            logger.debug("Read {} entries from {}".format(len(contigs), contigs_fi))
            contigs.index.name = "pos"
            contigs.reset_index(inplace=True)
            store.append('contigs', contigs, format='table')
        else:
            raise ValueError("Need contig info file")

    pcs = None 
    if "pcs" not in store or force:  # and not pcs_fi:
        pcs = pd.read_csv(pcs_fi, engine='python', sep=',')
        logger.debug("Read {} entries from {}".format(len(pcs), pcs_fi))
    elif "pcs" in store:
        pcs = store.pcs
        
    if not os.path.exists(os.path.join(folder, "profiles.pkle")) or force:
        if pcprofiles_fi is not None:
            profiles = pd.read_csv(pcprofiles_fi, engine='python', sep=',')
            profiles['contig_id'] = profiles['contig_id'].str.replace(' ', '~')

            if pcs is None:
                pcs = pd.DataFrame(profiles.pc_id.drop_duplicates())
                pcs.columns = ["pc_id"]

            # Filtering the PC profiles that appears only once
            before_filter = len(profiles)
            cont_by_pc = profiles.groupby("pc_id").count().contig_id.reset_index()

            # get the number of contigs for each pcs and add it to the dataframe
            cont_by_pc.columns = ["pc_id", "nb_proteins"]
            pcs = pd.merge(pcs, cont_by_pc, left_on="pc_id", right_on="pc_id", how="left")
            pcs.fillna({"nb_proteins": 0}, inplace=True)

            # Drop the pcs that <= 1 contig from the profiles.
            pcs = pcs[pcs['nb_proteins'] > 1]  # .query("nb_contigs>1")
            at_least_a_cont = cont_by_pc[cont_by_pc['nb_proteins'] > 1]  # cont_by_pc.query("nb_contigs>1")

            # profiles = profiles.query("pc_id in at_least_a_cont.pc_id")
            profiles = profiles[profiles['pc_id'].isin(at_least_a_cont.pc_id)]

            logger.info("Read {} entries (dropped {} singletons) from {}".format(
                        len(profiles), (before_filter-len(profiles)), pcprofiles_fi))
        
            pcs = pcs.reset_index(drop=True)
            pcs.index.name = "pos"
            pcs = pcs.reset_index()
            
            matrix, singletons = _matrix(profiles, store.contigs, pcs)  # Build PC profiles matrices from 3 df's

            # save
            profiles = {"matrix": matrix, "singletons": singletons}
            store.append('pcs', pcs, format='table')
            with open(os.path.join(folder, "profiles.pkle"), "wb") as f:
                pickle.dump(profiles, f)
        else:
            raise ValueError("Need profiles file")
        
    else:  # If the pickle file exist
        with open(os.path.join(folder, "profiles.pkle"), "r") as f:
            profiles = pickle.load(f)

        matrix, singletons = profiles["matrix"], profiles["singletons"]

    logger.info("{} contains : \n {:10} contigs, \n {:10} protein-clusters, \n {:10} singletons \n {:10} "
                "profile-entries.".format(folder, len(store.contigs), len(store.pcs), singletons.sum(),
                                          matrix.getnnz()))

    if (len(store.contigs) != matrix.shape[0] or (len(store.contigs) != singletons.shape[0])
        or (len(store.pcs) != matrix.shape[1])):
        logger.error("profile matrix: {}, singletons matrix {}".format(matrix.shape, singletons.shape))
        logger.debug(matrix.todense())
        logger.debug(singletons.todense())
        raise ValueError("Number of contigs and/or PC inconsistent with the profiles data.")

    return store.contigs, store.pcs, profiles


def _matrix(profiles, contigs, pcs):
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
