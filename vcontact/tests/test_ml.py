""" Unit test for the ml_functions module"""
from .. import ml_functions
import numpy.testing 
import numpy as np
import pandas 

F = {} # Fixtures 
def setup():
    F["classification"] = pandas.DataFrame({"reference":[1     ,1     ,1,1,1,np.nan,2,3,2,8,np.nan],
                                            "predicted":[np.nan,np.nan,2,1,1,3     ,2,1,2,3,np.nan]})
    F["metrics"] = {"entries": 11,
                    "common_classes":3,
                    "reference_entries":9,
                    "classified_entries":8,
                    "reference_classes":4,
                    "classified_classes":3,
                    "precision": 1.0/3,
                    "recall":7.0/20,
                    "specificity":195.0/224,
                    "fmeasure":14.0/41,
                }

def test_metrics():
    metrics = ml_functions.classification_metrics(F["classification"])
    for k in metrics.keys():
        numpy.testing.assert_almost_equal(metrics[k],F["metrics"][k],err_msg=": "+k)
