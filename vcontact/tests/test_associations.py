""" Unit test for the association module"""
from .. import associations
import numpy.testing 
import numpy as np
import pandas 
import pandas.util.testing as ptest

F = {} # Fixtures 
def setup():
    F["contigs"] = pandas.DataFrame({"pos":range(10),
                                     "pos_cluster_mbship":[10, 12, 14, 16, 18, 11, 13, 15, 17, 19]})
    F["clusters"] = pandas.DataFrame({"pos":[x + 10 for x in range(10)],
                                      "pos_family":[x + 20 for x in range(10)[::-1]]})
    F["taxonomy"] = pandas.DataFrame({"pos":[x + 20 for x in range(10)],
                                      "name":["class_{}".format(x + 20)for x in range(10)[::-1]]})
    F["contigs"].ix[2,"pos_cluster_mbship"] = np.nan
    
    F["clusters"].ix[5,"pos_family"] = np.nan
    F["contigs_out"] = pandas.DataFrame({"pos":F["contigs"].pos.values,
                                         "pos_cluster_mbship":F["contigs"].pos_cluster_mbship.values,
                                         "predicted_family":["class_{}".format(x + 20)for x in range(10)[::2]+range(10)[1::2]]})
    
    F["contigs_out"].ix[7,"predicted_family"] = "Non affiliated"
    F["contigs_out"].ix[2,"predicted_family"] = "Non affiliated"
    
def test_contig_taxonomy():
    print "FIXTURES:"
    print F["contigs"], "(contigs)"
    print F["clusters"], "(clusters)"
    print F["taxonomy"], "(taxonomy)"

    print "EXPECTED:"
    print F["contigs_out"]
    print F["contigs_out"].dtypes
    contigs = associations.contig_taxonomy(F["contigs"], F["taxonomy"], F["clusters"], "family",
                                           cluster_choice="pos_cluster_mbship")
    print "OUTPUT:"
    print contigs
    print contigs.dtypes
    ptest.assert_frame_equal(contigs,F["contigs_out"],check_dtype=False)
                             
                                            
    
