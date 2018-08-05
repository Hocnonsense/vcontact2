.. _tutorial:

vContact Tutorial
=================

This tutorial explain how to analyze a virome for wich you already
have a pc-profile. To use vContact/MCL to generate the profiles from
the blast results see this (not written yet) :ref:`protein
clustering tutorial <tutorial_pcs>`.

Prepare your data for input 
----------------------------

You need at least a file contiaining the *pc-profiles* (association of
contigs and clusters) and one containing information about the
*contigs* (the number of predicted proteins is the only required
field). You can also optionnaly provide a file describing the modules.

The supported fields are :

================ ============== ========= ==========================================
File              Field          Required  Description
================ ============== ========= ==========================================
**profiles**     contig_id      yes       uniq contig name.
--               pc_id          yes       uniq protein cluster name.
**contigs info** id             yes       uniq contig name.
--        	     proteins       yes       number of protein on the contig.
--               size           no        size of the contig bp.
--     	    	 origin         no        origin of the contig (defaut is filename).
--               *<others>*     no        reference taxonomy.
**pc info**      id             yes       uniq protein cluster name.
--               size           no        number of proteins.
--               *<others>*     no        additional description.
================ ============== ========= ==========================================


Example files:

`profiles.csv`

========= ==================
contig_id pc_id
========= ==================
NC_550500 PC_000001
NC_550500 PC_000002
NC_280000 PC_000002
NC_280000 PC_000004
NC_280000 PC_000006
========= ==================


`contigs.csv`

========= ======== ================== =============
id        proteins family             origin
========= ======== ================== =============
NC_550500       10  Myoviridae        refseq_jan14
NC_280000        9  Myoviridae        refseq_jan14 
NC_550000       55  Podoviridae       refseq_jan14
000000001       14                    myVirome
========= ======== ================== =============

`protein clusters.csv`

================== =============
id                 keywords
================== =============
PC_000001              reductase   
PC_000002           
PC_000002
PC_000004
PC_000006                  tail
================== =============

.. WARNING:: For a given contig, if the number of proteins in
   `contigs.csv` file is higher than the number of times it appears in
   the `profiles.csv` file, the remaining proteins are treated as
   unclustered singletons (and change the significativity in the
   similarity network). To account for multiple proteins of one
   cluster on a single contig, repeat the line in the `protein.csv`
   file. (It will only change the number of singletons).

.. NOTE:: Proteins clusters that do not appear (or only once) in
   `profiles.csv` are considered as singletons not stored.

   
Run the analysis
----------------

Simplest way
^^^^^^^^^^^^
The main command provided by vcontact is ::

  vcontact -p profiles.csv -c contigs.csv --pcs protein_clusters.csv -o output_dir 

If the output directory does not exist, it will be created. The output is structured as follow::
	
  output_dir/
	contigs.h5
	protein_clusters.h5
	clusters_i2_sig1.h5
	modules_i5_sig1.h5
	profiles.pkle
	clusters_i2_sig1.csv
	modules_i5_sig1.csv
	affiliation_i2_sig1.csv

The results are stored in the csv files and you can parse them with
your favorite tools to visualise them. vContact provides also several
exporting functions to speed up the process for common tools see :ref:`export`.
	
.. NOTE::
   The `.h5` (HDF5 binary files) and `.pkle` (Python serialized
   object) are used to speed up new computation on the same output
   directory. If they exists you do not need to precise the arguments
   `-p`, `-c` and `-pcs` anymore when working with the same `-o`.  You
   can use the option `-f` to force their overwriting if you have
   changed the data without changing the output directoy name.
	   
Change parameters
^^^^^^^^^^^^^^^^^

You can change the following parameters for the contigs clustering :

  * `-i` : Inflation in the MCL clustering of the contigs. (default 2)
  * `-sig` : Significativity threshold (default 1)

For the modules:
	
  * `--i_mod` : Inflation in the MCL clustering of the modules (default 5)
  * `--sig_mod` : Significativity threshold
  * `--shared_min` : Minimal number of contigs sharing a protein cluster to be used in the modules (default 3)

For the linking of modules and clusters
  * `--own_thres` : Minimal proportion of protein cluster of the module to consider that the module is present.
  * `--assoc-sig` : Significativity threshold

.. _export:
	
Exporting
--------

.. Warning:: Not implemented in the command line yet. You can use
  these function if your are using vcontact in an interactive python
  session (see the API)

You can ask for an export in a folder where the computation has already been done ::

  vcontact -o output_dir --export cytoscape
  vcontact -o output_dir --export cytoscape=cluster_id
  vcontact -o output_dir --export krona





