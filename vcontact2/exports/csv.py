""" Functions to export as csv """
import os
import pandas
import numpy as np


def complete(folder, pcm=None, gc=None, mod=None, link=None):
    """ Export all the things!
   
    Args:
        pcm (vcontact.pcprofiles.PCProfiles): Profiles
        gc (vcontact2.genome_clusters.GenomeClusters): Contigs
        mod (vcontact2.modules.Modules object): Modules
        link (pandas.DataFrame): Link between contig clusters and modules
    """

    permissive = ""
    if gc is not None:
        if gc.permissive:
            permissive = "_permissive"

        fn = "sig{}_mcl{}{}".format(gc.thres, gc.inflation, permissive)
        gc.contigs.to_csv(os.path.join(folder, "{}_contigs.csv".format(fn)), index=False)
        gc.clusters.to_csv(os.path.join(folder, "{}_clusters.csv".format(fn)), index=False)

    if mod is not None:
        fn = "sig{}_mcl{}_minshared{}".format(mod.thres, mod.inflation, mod.shared_min)
        mod.modules.to_csv(os.path.join(folder, "{}_modules.csv".format(fn)), index=False)

    if link is not None:
        fn = "sig{}_mcl{}_modsig{}_modmcl{}_minshared{}".format(gc.thres, gc.inflation, mod.thres, mod.inflation, mod.shared_min)
        link.to_csv(os.path.join(folder, "{}_link_mod_cluster.csv".format(fn)), index=False)


def summary(folder, gc=None, pcm=None, category="RefSeq-85"):
    
    output = []
    contigs = None
    permissive = ""
    if gc is not None:
        if gc.permissive:
            permissive = "_permissive"

        contigs = gc.contigs.set_index("pos")
        pos_refseq = contigs.query("origin==category").index
        pos_tara = contigs.query("origin!=category").index

        output.append(("{tot} Contigs: {nb_cat} in {cat} and {nb_ncat} not in {cat}"
                       "").format(nb_ncat=len(pos_tara), nb_cat=len(pos_refseq),
                                  tot=len(pos_tara)+len(pos_refseq),
                                  cat=category))

        clusters_refseq = frozenset(contigs.loc[pos_refseq, "pos_cluster"].drop_duplicates())
        clusters_tara = frozenset(contigs.loc[pos_tara, "pos_cluster"].drop_duplicates())
        print(("{tot} Clusters: {nb_cat} with a {cat} sequence "
               "{nb_ncat} with a non {cat} sequence"
               "").format(cat=category,
                          nb_cat=len(clusters_refseq),
                          nb_ncat=len(clusters_tara),
                          tot=len(contigs.pos_cluster.drop_duplicates())))

    if pcm is not None:
        if contigs is None:
            contigs = pcm.contigs.set_index("pos")
            pos_refseq = contigs.query("origin==category").index
            pos_tara = contigs.query("origin!=category").index

        refseq_pcs = (np.squeeze(np.asarray(pcm.matrix[pos_refseq,:].tocsc().sum(0))) > 0).sum()
        tara_pcs = (np.squeeze(np.asarray(pcm.matrix[pos_tara,:].tocsc().sum(0))) > 0).sum()
        commons_pcs = int(np.dot(np.asarray((np.asarray(pcm.matrix[pos_refseq, :].tocsc().sum(0)) > 0), dtype=int),
                                 np.asarray(np.transpose((np.asarray(pcm.matrix[pos_tara, :].tocsc().sum(0)) > 0)), dtype=int)
                             ))
        print(("{tot} shared PCs (+ {sing}, for a total of {allpcs}): "
               "{nb_cat} pcs in {cat} sequences and "
               "{nb_ncat} pcs in a non {cat} sequences "
               "and {common} in common").format(nb_cat=refseq_pcs,
                                                nb_ncat=tara_pcs,
                                                common=commons_pcs,
                                                tot=pcm.matrix.shape[1],
                                                sing=pcm.singletons.sum(),
                                                allpcs=pcm.matrix.shape[1]+pcm.singletons.sum(),
                                                cat=category))

    fp = os.path.join(folder, "sig{}_mcl{}{}".format(gc.thres, gc.inflation, permissive))
    with open(fp, "w") as f:
        f.write("\n".join(output))
