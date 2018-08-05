.. _tutorial_pcs:

Protein clustering tutorial
===========================

This tutorial explain how to generate PC-profiles using
vContact/MCL. To learn how to analyse the profile see the main :ref:`tutorial`

Prepare your data for input
----------------------------

You will need a file giving the association between each protein and
each contig. Optionally, this file can contains a function field
containing keywords (separated by a semi column) about the functional
annotation of this protein.  If it is the case the resulting protein
clusters will be annotated with the count of the occurance of those
keywords.

========== ========= =================
id [#f1]_  contig_id keywords
========== ========= =================
YC.330JH1E NC.000001 tail; fiber;
YC.6567899 NC.000001
YC.789U666 NC.000002 helicase
========== ========= =================

.. [#f1] Exact same id as in the fasta/blast file.

You will also need a fasta file or output of a all versus all BLASTp.

Run the analysis
----------------

To generate the profiles from a already computed blast, use the command::

  vcontact-pcs -p proteins.csv -b blastresults.tab   -o output_dir

To ask vcontact to run the blast for you use ::

  vcontact-pcs -p proteins.csv -f fasta.faa  -o output_dir -e 0.0005


In both case it will output ::

	output_dir/
		contigs.csv
		pcs.csv
		profiles.csv

See the :ref:`tutorial` to see what they contains.
