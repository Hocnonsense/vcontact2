""" Useful functions and snippets """
import networkx
import numpy as np


def summary(matrix, nodeinfo, criterion="origin"):
    """ Summary of some graph proterties
    
    Args:
        matrix: (sparse matrix) adjacency matrix
        nodeinfo: (dataframe) with:
            name: node name (str).
            pos: node position in the matrix (int)
            "criterion": a column to distinguish categories of nodes. 
        criterion: column to use to distinguish nodes. 
    Returns:
        A dict. containing informations about the graph.
    """

    info = {}
    categories = nodeinfo[criterion].drop_duplicates().values
    
    graph = networkx.from_scipy_sparse_matrix(matrix)
    #networkx.relabel_nodes(graph,
    #                       dict(zip(nodeinfo.pos,
    #                                nodeinfo.name)),
    #                       copy=False)

    info["nodes"] = graph.order() 
    info["edges"] = graph.size()

    cc = networkx.connected_components(graph)
    info["connected_components"] = len(cc)
    
    lcc = [len(x) for x in cc]
    info["non_connected_nodes"] = lcc.count(1)
    info["biggest_connected_component"] = max(lcc)

    non_connected = [x[0] for x,len_ in zip(cc,lcc) if len_ == 1]
    info_non_connected = nodeinfo.set_index("pos").loc[non_connected,:].fillna({criterion: "No category"})
                                                                               
    info["non_connected_by_{}".format(criterion)] = dict(info_non_connected.groupby(criterion).index.count())

    return info 


def filtering(matrix, threshold):
    """ Filter a sparse matrix according to a threshold
    
    Args:
        matrix: (sparse matrix)
        threshold: (float)
    
    Returns:
        The sparse matrix where every cell <= threshold are now null
    """
    return matrix.multiply(matrix>=threshold)


def summary_assoc(matrix, info, criterion=None, lines=True, thres=0, name_thres=10):
    """ """
    report = {}
    if not lines:
        matrix = matrix.T
    boolmatrix = matrix > thres
    try:
        report["entries"] = matrix.getnnz()
        report["assoc"] = boolmatrix.getnnz()
    except AttributeError:
        report["entries"] = matrix.shape[0] * matrix.shape[1]
        report["assoc"] = boolmatrix.sum()
    report["objects"] = matrix.shape[0]
    report["ref_objects"] = len(info)

    info = info.set_index("pos")
    assoc_by_element = boolmatrix.sum(1)
    
    non_connected = [n for n,x in enumerate(assoc_by_element) if x==0]
    report["associated_to_nothing"] = len(non_connected)

    if report["associated_to_nothing"]<= name_thres:
        report["associated_to_nothing_names"] = [info.name[p] for p in non_connected] 

    report["most_associated_links"] = np.max(assoc_by_element)
    pos_ma = np.argmax(assoc_by_element)
    report["most_associated_name"] = info.name[pos_ma]
    if criterion is not None:
        info_non_connected = info.loc[non_connected,:].fillna({criterion: "No category"})
        report["associated_to_nothing_by_{}".format(criterion)] = dict(info_non_connected.groupby(criterion).index.count())
        report["most_associated_{}".format(criterion)] = info.loc[np.argmax(assoc_by_element),criterion]

    return report
