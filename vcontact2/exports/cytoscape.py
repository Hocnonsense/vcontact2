"""Export functions to use the data generated with cytoscape.
http://www.cytoscape.org/"""
import logging
from itertools import combinations
import csv

import scipy.sparse as sparse
import numpy as np
import pandas

from .. import options

logger = logging.getLogger(__name__)


def contigs(network, contigs, clusters, path):
    """ Export the network and node information in a format
    that can be used with cytoscape (tab files).

    Args:
        network (sparse-matrix): Similarity network.
        contigs (dataframe): Nodes information.
        clusters (list): Id of clusters to exports. 
        path (str): basename of the exported files.

    Returns:
        tuple: Tuple of the path to the edge files and 
            nodes information files
    """
    fi_ntw = path+".ntw"
    fi_info = path+".info"
    line = 0
    
    contigs = contigs.query("pos_cluster in clusters")

    logger.debug(("Exporting the network ({} contigs, "
                 "{} clusters)").format(len(contigs),len(clusters)))

    contigs.set_index("name").to_csv(fi_info,sep="\t")
    logger.info("Wrote {} lines in {}.".format(len(contigs),fi_info))    

    with open(fi_ntw,"w") as f:
        f.write("Node_1\t Type\t Node_2\t Edge_weight \n")
        for r,c in combinations(contigs.pos,2):
            if network[r,c]:
                n1 = contigs.ix[r,"name"] 
                n2 = contigs.ix[c,"name"] 
                f.write("\t".join([str(n1),
                                   "contig-contig",
                                   str(n2),
                                   str(network[r,c])]))
                f.write("\n")
                line += 1
    logger.info("Wrote {} lines in {}.".format(line,fi_ntw))

    return fi_ntw,fi_info


def membership(fi, B, contigs, clusters, clusters_list=None,
               criterion="predicted_family"):
    """Save the membership matrix in an input format for cytoscape (csv)

    Args:
        B (sparse_matrix): Membership matrix.
        contigs (pandas.df): Contigs informations (required columns: 
            pos, pos_cluster).
        clusters (pandas.df): Clusters information (required columns:  pos,
            criterion).
        cluster_list (list): The output is restircted to thoses clusters.
        criterion (str): Column on wich group the membership. 
        fi (str): Path to the output file.
    """
    if clusters_list is not None:
        contigs = contigs.query("pos_cluster in clusters_list")
    contigs.sort("pos",inplace = True)

    connected_clusters = [i for i,x
                          in enumerate(np.squeeze(np.asarray(B[contigs.pos,:].sum(0)))) if x]

    clusters = clusters.query("pos in connected_clusters")

    groups = clusters[criterion].drop_duplicates()

    for group, data in clusters.groupby(criterion):
        contigs[group] = B[contigs.pos,:][:,data.pos].sum(1)

    with open(fi,"wb") as f:
        f.write("Node\t")
        f.write("\t".join(groups))
        f.write("\tChart")
        f.write("\n")
        for n,l in contigs.iterrows():
            f.write("{0}\t ".format(l["name"]))
            f.write("\t".join([str(l[g]) for g in groups]))
            f.write('\t piechart: attributelist="{0}" showlabels=false colorlist="{1}"'.format(",".join(groups),",".join([options.COLORS[g] if g in options.COLORS else options.COLORS["pink"] for g in groups])))
            f.write("\n")
    logger.debug("Wrote {} lines in {}".format(n,fi))


def clusters(cluster_network, contigs, clusters, criterion, path):
    """ Export the network and node information in a format
    that can be used with cytoscape (tab files).

    Args:
        cluster_network (sparse-matrix): cluster network.
        contigs (dataframe): Nodes information.
        clusters (list): Id of clusters to exports. 
        path (str): basename of the exported files.
        criterion (str): criterion to put on the piechart.
    Returns:
        tuple: Tuple of the path to the edge files and 
            nodes information files
    """
    fi_ntw = path+".ntw"
    fi_info = path+".info"
    line = 0
        
    gbc = contigs.dropna(subset=[criterion]).loc[:,[criterion,"pos_cluster"]].groupby("pos_cluster")
    keys = []
    lines = []
    ref = contigs.query("origin=='refseq_jan14'")
    
    for clust, data in gbc:
        keys.append(clust)
        lines.append(data.groupby(criterion).count()[criterion])
    criterion_by_cluster = pandas.DataFrame(lines,keys)
    criterion_by_cluster.fillna(0, inplace=True)
    criterion_by_cluster["Graph"] = ('piechart: attributelist="{0}" showlabels=false "'
                                     'colorlist="contrasting"').format(",".join(criterion_by_cluster.columns))
    info = pandas.merge(clusters, criterion_by_cluster,
                        how="left",
                        left_on="pos",right_index=True)
    
    info.set_index("name").to_csv(fi_info,sep="\t",quote="")

    clusters = clusters.set_index("pos")
    with open(fi_ntw,"w") as f:
        f.write("Cluster_1\tCluster_2\tinter_ov_intra\n")
        for r in range(cluster_network.shape[0]):
            for c in range(r):
                if cluster_network[r,c]:
                    f.write("\t".join([clusters.ix[r,"name"],
                                       clusters.ix[c,"name"],
                                       str(cluster_network[r,c])]))
                    f.write("\n")        
                    line += 1
                    
    logger.info("Wrote {} lines in {}.".format(line,fi_ntw))
    return fi_ntw,fi_info





'''
def cluster_network_cytoscape(self,fi_ntw=None,fi_clust_info=None,network=None):
        #TODO: REFACTOR THIS 

        fi_ntw = self.name+"_contigsclusters_network.ntw"
        fi_clust_info = self.name + "_contigsclusters_network.info"
        info = self.clusters

        
        info = pandas.merge(info,self.taxonomy.query("level=='family'"),
                            how="left",
                            left_on="pos_family",right_on="pos",
                            suffixes=["","__family"]
                            )
        #print info
        info = pandas.merge(info,self.taxonomy.query("level=='genus'"),
                            how="left",
                            left_on="pos_genus",right_on="pos",
                            suffixes=["","__genus"]
                            )
        info.reset_index(inplace=True)
        info.set_index("name",inplace=True)
        #print info
        info = info.loc[:,["size","name__family","name__genus"]]
        info.to_csv(fi_clust_info,sep="\t")
       
        network = self.network if network is None else network
  

        L = self.link_clusters(network=network)
        
        cluster_idx = self.clusters.reset_index()
        cluster_idx.set_index("pos",inplace=True)
        #print cluster_idx
        with open(fi_ntw,"w") as f:
            f.write("Cluster_1\t")
            f.write("Cluster_2\t")
            f.write("inter_ov_intra")
            f.write("\n")
            for r in range(L.shape[0]):
                for c in range(r):
                    if L[r,c]:
                        n1 = cluster_idx.ix[r,"name"] if r != L.shape[0]-1 else "non-clustered" 
                        n2 = cluster_idx.ix[c,"name"] if c != L.shape[0]-1 else "non-clustered"
                        f.write("\t".join([str(n1),
                                           str(n2),
                                           str(L[r,c])]))
                        f.write("\n")        

def export_network(fi_ntw, network, dataframe,
                   pos, col_cluster,
                   other_clusters=None, membership=None):
    """ Export the network in a cytoscape file """ 
    line =1
    with open(fi_ntw,"w") as f:
        f.write("Node_1\t")
        f.write("Type\t")
        f.write("Node_2\t")
        f.write("Edge_weight")
        f.write("\n")
        for r,c in combinations(pos,2) :
                if network[r,c]:
                    n1 = dataframe.ix[r,"name"] 
                    n2 = dataframe.ix[c,"name"] 
                    f.write("\t".join([str(n1),
                                       "contig-contig",
                                       str(n2),
                                       str(network[r,c])]))
                    f.write("\n")
                    line +=1
        if other_clusters != None and membership != None:            
                for c in other_clusters:
                    pos_c = list(dataframe.query("{}=={}".format(col_cluster,c)).index)
                    if network[pos_c,:][:,pos_c].sum():
                        for n in pos:
                            s = network[n,pos_c].sum()
                            if s and membership[n,c] > 0.1:
                                f.write("\t".join([dataframe.ix[n,"name"],
                                                   "contig-cluster",
                                                   "cluster_{}".format(c),
                                                   str(s)]))
                                f.write("\n")
                                line += 1
    print("Wrote {} lines in {}.".format(line,fi_ntw))


def membership(B, contigs, fi):
    """Save the membership matrix in an input format for cytoscape (csv)
    INPUT:
    - genome_cluster (GenomesCluster) object
    - fi (str) filename (default self.name without parenthesis) 
    """

    groups = ["g{0}".format(n) for n in range(B.shape[1])]

    with open(fi,"wb") as f:
        f.write("Node\t")
        f.write("\t".join(groups))
        f.write("\tChart")
        f.write("\n")
        for r in range(B.shape[0]):
            f.write("{0}\t ".format(contigs.query("pos=={}".format(r)).name.values[0]))
            for c in range(B.shape[1]):
                f.write("{}\t".format(B[r,c]))
            f.write('piechart: attributelist="{0}" colorlist="contrasting"'.format(",".join(groups)))
            f.write("\n")


def network(network, dataframe, fi,
                      pos=None, cluster=None,
                      col_cluster="mcl_cluster", col_pos="pos", col_name="name",
                      membership=None):
    """Save the network defined by the matrix with
    the names given by dataframe into a file fi for cytoscape
    
    the position in the matrix and the node names are given by
    the pos_col and name_col resp.
    
    Save also all the informations about a node given by the dataframe
    in a info file"""

    fi_ntw = fi+"_cytoscape.ntw"
    fi_info = fi + "_cytoscape.info"
    fi_mb = fi + "_membership_cytoscape.info"
    fi_info_clust = fi + "_clusters_cytoscape.info"

    ### Get the position of the contigs in the selected clusters 
    dataframe = dataframe.reset_index().set_index(col_pos)
    other_clusters = list(dataframe.loc[:,col_cluster].drop_duplicates().dropna().values)
    if pos == None:
        if cluster != None:
                pos = [] 
                for c in cluster:
                        pos += list(dataframe.query("{}=={}".format(col_cluster,c)).index)
                        other_clusters.pop(other_clusters.index(c))
        else:
            raise ValueError("Need pos list or cluster index")
    print "{} contigs in those {} cluster(s), {} other clusters.".format(len(pos),len(cluster),len(other_clusters))

    ### Export the network 
    export_network(fi_ntw,network,dataframe,
                   pos,col_cluster,
                   other_clusters,membership) 

    ### Export the membership info 
    if membership != None :
        export_membership(fi_mb,membership,dataframe,pos)
    
    ### Export node informations
    dataframe = dataframe.loc[pos,:]
    dataframe.index.name = "pos2"
    dataframe = dataframe.reset_index().set_index(col_name)
    dataframe.to_csv(fi_info,sep="\t")
    print("Wrote {} lines in {}.".format(len(dataframe),fi_info))    

    ### Export clusters informations
    info_clust = pandas.DataFrame(dataframe.groupby("mcl_cluster").size(),columns=["size"])
    info_clust.to_csv(fi_info_clust,"\t")
    print("Wrote {} lines in {}.".format(len(info_clust),fi_info_clust))    

def export_membership(fi_mb, membership, dataframe, pos):
    """ Export the membership matrix """ 
    membership = membership[pos,:]
    groups = [n for n,i in enumerate((np.squeeze(np.asarray(membership.sum(0)))>0)) if i]
    membership = membership[:,groups]
    groups = ["g_{}".format(g) for g in groups]

    with open(fi_mb,"wb") as f:
        f.write("Node\t")
        f.write("\t".join(groups))
        f.write("\tChart")
        f.write("\n")
        line = 1
        for i,r in enumerate(pos):
                f.write("{0}\t ".format(dataframe.ix[r,"name"]))
                f.write("\t".join([str(i) for i in np.squeeze(np.asarray(membership[i,:]))]))

        f.write('\tpiechart: attributelist="{0}" showlabels=false colorlist="contrasting"'.format(",".join(groups)))
        f.write("\n")
        line += 1
    print("Wrote {} lines in {}.".format(membership.shape[0],fi_mb))    
    
    '''


