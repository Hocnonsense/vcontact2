"Machine learning functions - Some helper functions from ML"
from sklearn.model_selection import StratifiedKFold
import pandas
import numpy as np


def split_dataset(dataset, criterion, fold):
    """Split the dataset in <fold> startified parts to use as training and
    cross validation test. Stratified means that the proportions of criterion
    in the original dataset are conserved as much as possible in the subsets.

    Args:
        dataset: (dataframe) the data.
        fold: (int) number of parts in which to split the dataset.
        criterion: (str) column on which perform the stratification.

    Returns:
        dataset: with the column "cvset_criterion" (int) in [0-fold-1] added.
            gives the id of the set this contig belong to. Set to nan to the
            contigs that do not meet the criterion.
    """

    cname = "cvset_{}".format(criterion)
    cv = StratifiedKFold(dataset[criterion], fold)
    for n, (training_set, validation_set) in enumerate(cv):
        dataset.ix[validation_set, cname] = n

    dataset.loc[pandas.isnull(dataset[criterion]), cname] = np.nan

    return dataset


def classification_metrics(dataset, ref_col="reference", pred_col="predicted"):
    """ Returns a series of metric given a classification.

    The metrics are "macro-averaged", which means that we take the average
    of metric between the classes. As a consequence each class as the 
    same weight in the final result regardless of its size.
    
    Non affiliated references (NaN) are dropped.
    
    Predicted classes absent from the references or that are not defined (NaN)
    are counting as negative prediction for all classess (and thus dispatched
    between True and False negatives). 

    By convention if the divider is null the metric is null.

    Args:
        dataset (pandas.DataFrame): one row by entry.
            reference: the reference classification of the entry
            predicted: the predicted classification of the entry
    ref_col (str): column where the reference classes are.
        pred_col (str): column where the predicted classes are.

    Returns:
        dict: A dictionary including:
            entries: number of entries.
            reference_entries: number of entries with a reference class.
            reference_classes: number of different classes. (used for 
                the macro averaging).
            common_classes: number of classes found in reference and prediction. 
            classified_entries: number of entries with a predicted class.
            classified_classes: number of different predicted classes. 
            recall: TP/(TP+FN) (positive effectively predicted as positive). 
                (a.k.a, sensitivity, power or 1-type II error).
            precision: TP/(TP+FP) (predicted positive effectively positive).
                (a.k.a. positive predictive value).
            specificity: TN/(FP+FN) negative effectively predicted as negative.
                (a.k.a 1-type I error)
            fmeasure: 2PR/(P+R) harmonic mean of precision and recall.
    """

    ref = ref_col
    pred = pred_col
    out = {}

    out["entries"] = len(dataset)

    out["reference_entries"] = dataset[ref].count()
    out["classified_entries"] = dataset[pred].count()

    ref_class = dataset[ref].drop_duplicates()
    pred_class = dataset[pred].drop_duplicates()
    out["reference_classes"] = dataset[ref].drop_duplicates().count()
    out["classified_classes"] = dataset[pred].drop_duplicates().count()

    out["common_classes"] = len(frozenset(ref_class) & frozenset(pred_class))
    
    dataset = dataset.dropna(subset=[ref])
 
    P, R, S = 0.0, 0.0, 0.0
    if out["reference_classes"]:
        pred_count = dataset.groupby(pred).size()  # NO LONGER dataset.groupby(pred).count()[pred]
        for cat, data in dataset.groupby(ref):
            TP = np.sum(data[pred]==cat)
            FN = np.sum(data[pred]!=cat)
            if cat in pred_count.index:
                PP = pred_count[cat] 
            else:
                PP = 0
            FP = PP - TP
            TN = out["reference_entries"] - PP - FN

            P += np.nan_to_num(TP/float(TP+FP))
            R += np.nan_to_num(TP/float(TP+FN))
            S += np.nan_to_num(TN/float(FP+TN))

        out["precision"] = P/out["reference_classes"]
        out["recall"] = R/out["reference_classes"]
        out["specificity"] = S/out["reference_classes"]
        out["fmeasure"] = 2 * (out["precision"] * out["recall"])/(out["precision"]+out["recall"])
    else:
        out["precision"] = 0
        out["recall"] = 0
        out["specificity"] = 0
        out["fmeasure"] = 0

    return out