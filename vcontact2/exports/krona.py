"""Export functions for krona.
http://sourceforge.net/projects/krona/"""

import logging

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

logger = logging.getLogger(__name__)


def textfile(dataframe, fi, columns=("predicted_family", "predicted_genus"),
             placeholder="Non affiliated"):
    """ Write a text file from a dataframe to be used with Krona.

    Args:
        dataframe (pandas.DataFrame): data. 
        fi (str): path to the file to write.
        columns (iterable): Columns to use as wedges levels. 
        placeholder (str): String to use when the value is N/A.

    Returns:
        str: Path to the file. 

    Note:
        Use `ktImportText test_krona` 
    """
    fillna = dict([(x, placeholder) for x in columns])
    count = dataframe.fillna(fillna).groupby(columns).count()[columns[0]]
    output = StringIO()
    count.to_csv(output, sep="\t")

    # Save with the count in first position. 
    with open(fi, "w") as f:
        for j,x in enumerate(output.getvalue().split("\n")):
            x = x.split("\t")
            f.write("\t".join([x[-1]]+x[:-1]))
            f.write("\n")
    logger.debug("Saved {} lines to {}".format(j,fi))
    return fi 
