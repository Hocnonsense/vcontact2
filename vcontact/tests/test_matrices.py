import numpy as np
import scipy.sparse as sparse

from .. import matrices
F= {}
def setup():
    F["network_hypergeom"] = sparse.lil_matrix(np.matrix([[0,0,1,1,0,0],
                                                          [0,0,0,0,0,1],
                                                          [1,0,0,1,0,0],
                                                          [1,0,1,0,0,0],
                                                          [0,0,0,0,0,0], 
                                                          [0,1,0,0,0,0]]))
    
    F["B"] = np.matrix([[1,0],
                        [0,1],
                        [1,0],
                        [1,0],
                        [0,0],
                        [0,1]],dtype=float)

    F["K_family"] = sparse.lil_matrix(np.matrix([[1,0,0],
                               [0,1,0],
                               [1,0,0],
                               [1,0,0],
                               [0,0,1],
                               [0,1,0]]))
    F["K_genus"] = sparse.lil_matrix(np.matrix([[1,0,0,0],
                              [0,1,0,0],
                              [0,0,0,0],
                              [1,0,0,0],
                              [0,0,1,0],
                              [0,0,0,1]]))
    F["Q_genus"] = np.matrix([[2,0,0,0],
                              [0,1,0,1]])
                                
    F["R_genus"] = np.matrix([[1,0,0,0],
                              [0,1,0,1]])
    
    F["P_genus"] = np.matrix([[1,0,0,0],
                              [0,0.5,0,0.5]])
    
    F["F_genus"] = np.matrix([[1,0,0,0],
                              [0,2./3.,0,2./3.]])

    F["Q_family"] = np.matrix([[3,0,0],
                               [0,2,0]])
                                
    F["R_family"] = np.matrix([[1,0,0],
                               [0,1,0]])

    F["P_family"] = np.matrix([[1,0,0],
                               [0,1,0]])

    F["F_family"] = np.matrix([[1,0,0],
                               [0,1,0]])


def test_correspondence():
    for level in ("family","genus"):
        obtained = matrices.correspondence(F["K_{}".format(level)], F["B"])
        wanted = [F["{}_{}".format(i,level)] for i in ("Q", "R", "P", "F")]
        print(level)
        for i,m in enumerate(("Q", "R", "P", "F")):
            print(m)
            print("Obtained: \n {}".format(obtained[i]))
            print("Wanted: \n {}".format(wanted[i]))
            np.testing.assert_array_equal(obtained[i],wanted[i])
