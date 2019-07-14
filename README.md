vConTACT : Viral CONtig Automatic Cluster Taxonomy 2
====================================================

vConTACT2 is a tool to perform guilt-by-contig-association automatic classification of viral contigs.

## Requirements

vConTACT requires numerous python packages to function correctly, and each must be properly installed and working for 
vConTACT to also work.

 * python ~3.6 (not python 2.7 safe)
 * networkx>=1.11
 * numpy>=1.12.1
 * scipy>=0.19.0
 * pandas>=0.19.2
 * scikit-learn>=0.18.1
 * biopython>=1.68
 * hdf5>=1.8.17
 * pytables>=3.3.0
 * psutils

vConTACT also requires several executables, depending on use.

 * MCL (always required)
 * BLASTP (only if using BLASTP for PC construction)
 * DIAMOND (only if using DIAMOND for PC construction)
 * ClusterONE (only if using for PC or VC construction)

Generally you want these tools to be in your system or user PATHs. vConTACT will use any user-provided paths before 
searching through system PATHs.

Hardware requirements *can be considerable* (exceeding 48 GB!), depending  mainly on the size and complexity of the 
dataset. (Relationship between memory requirements and sequences analyzed forthcoming)

## Installation

Installing vConTACT dependencies may seem daunting, but the instructions below should work for the vast majority of 
users. While Windows is not officially supported, users have been able to run vConTACT on these machines (usually 
through some sort of virtual machine).

#### Singularity

A singularity image is provided accompanying this documentation. Please use this file to create and bootstrap the 
vConTACT container.

```bash
sudo singularity build vConTACT2.simg vConTACT2.def
```

The build process can take a significant amount of time depending on available hardware and network speed. Builds can 
take anywhere from 5 to 30 minutes. If you see *Finalizing Singularity container* at the end of bootstrapping, you're 
probably good to go.

Once built, the container can be run via:

```bash
singularity run vConTACT2.img <args>
```

#### Conda-based installation (Mac/Linux) (recommended)

We highly recommend using a python environment when installing this software, as dependency versions can (and often do)
 conflict with each other.

For this, we'll be installing everything into a single directory that the user has access to. This is generally the 
user's home directory. In the example, we'll be installing to the "conda" directory under the user's home folder.

First, grab our favorite manager, Anaconda/Miniconda and install. If it offers the option of adding the installation to 
the user's $PATH, then do so. Otherwise, follow the instructions (at the end of the install) to ensure the install is 
activated.

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# Install into $HOME/conda
conda create --name vContact2 python=3
source activate vContact2
```

Install python dependencies.

```bash
conda install -y -c conda-forge hdf5 pytables pypandoc biopython networkx numpy pandas scipy scikit-learn psutil
conda install -y -c bioconda mcl blast diamond
```

Install ClusterONE

```bash
wget http://www.paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar
cp cluster_one-1.0.jar $HOME/conda/bin/
```

Note: DIAMOND is highly recommended over BLASTP for any large-scale analysis. It’s much faster and shows 
little/no difference in the final VCs. *This hasn't been officially benchmarked, but a sufficient number of in-house 
analyses have been performed to recommend.*

Finally, install vConTACT2 from source file.

```bash
wget https://bitbucket.org/MAVERICLab/vcontact2/get/master.tar.gz
tar xvf MAVERICLab-vcontact2-XXXXXXX.tar.gz
cd MAVERICLab-vcontact2-XXXXXXX && pip install .
```

(Some users have mentioned that their version of pip installs to a non-conda location. In this case, run "pip install 
--prefix=$HOME/conda/bin)

Alternatively, install from bitbucket.

```bash
git clone bitbucket.org/MAVERICLab/vcontact2
cd vcontact2 && pip install .
```

#### Pip-based installation (Mac/Linux)

Pip should automatically install vConTACT2 alongside all of its dependencies, provided python3.x is correctly installed.

```bash
git clone bitbucket.org/MAVERICLab/vcontact2
cd vcontact2 && pip install .
```

You might encounter an issue where pip install doesn't install ALL of the database files. In this case, you'll have to 
manually copy the database files to wherever pip is installing vContact2 to.

```bash
cp vcontact2/vcontact/data/ViralRefSeq-prokaryotes-v??.* $HOME/conda/lib/python3.7/site-packages/vcontact/data/
```

*Your installation path might be at a different location.* Usually it's some form of 
"<where-you-installed-conda>/lib/python3.X/site-packages/vcontact/data/"

## Usage

#### For the impatient

Singularity

```bash
singularity run vConTACT2.img --raw-proteins [proteins file] --rel-mode ‘Diamond’ --proteins-fp [gene-to-genome mapping file] --db 'ProkaryoticViralRefSeq85-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin [path to ClusterONE] --output-dir [target output directory]
```

"--c1-bin" *must* point to the jarfile. For the singularity image, it's been set up to be located in /usr/local/bin/

Local installs

```bash
vcontact --raw-proteins [proteins file] --rel-mode ‘Diamond’ --proteins-fp [gene-to-genome mapping file] --db 'ProkaryoticViralRefSeq85-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin [path to ClusterONE] --output-dir [target output directory]
```

"--c1-bin" *must* point to the jarfile, not to the directory it's in.

## Input files and formats

vConTACT2 tries to alleviate some of the challenges created by the complex file format of the original vConTACT and 
provide a more seamless "pipeline" to go from raw data to a finished network. It also allows the user flexibility in 
providing files at various points along the processing pipeline. Generally speaking, if a user provides an intermediary 
file, vConTACT2 will skip the steps it would of taken to get to that file. This is useful because it allows partial 
recovery and re-start of any interrupted analysis. That said, *if there is an issue and you're trying to troubleshoot, 
always re-run and specify a new output directory*.

##### If starting from raw proteins

The only required files are:

1. A FASTA-formatted amino acid file.

```
>ref|NP_039777.1| ORF B-251 [Sulfolobus spindle-shaped virus 1]
MVRNMKMKKSNEWLWLGTKIINAHKTNGFESAIIFGKQGTGKTTYALKVAKEVYQRLGHE
PDKAWELALDSLFFELKDALRIMKIFRQNDRTIPIIIFDDAGIWLQKYLWYKEEMIKFYR
IYNIIRNIVSGVIFTTPSPNDIAFYVREKGWKLIMITRNGRQPDGTPKAVAKIAVNKITI
IKGKITNKMKWRTVDDYTVKLPDWVYKEYVERRKVYEEKLLEELDEVLDSDNKTENPSNP
SLLTKIDDVTR
>ref|NP_039778.1| ORF D-335 [Sulfolobus spindle-shaped virus 1]
MTKDKTRYKYGDYILRERKGRYYVYKLEYENGEVKERYVGPLADVVESYLKMKLGVVGDT
PLQADPPGFEPGTSGSGGGKEGTERRKIALVANLRQYATDGNIKAFYDYLMNERGISEKT
AKDYINAISKPYKETRDAQKAYRLFARFLASRNIIHDEFADKILKAVKVKKANADIYIPT
```

2\. A "gene-to-genome" mapping file, in either tsv (tab)- or csv (comma)-separated format.

```
protein_id,contig_id,keywords
ref|NP_039777.1|,Sulfolobus spindle-shaped virus 1,ORF B-251
ref|NP_039778.1|,Sulfolobus spindle-shaped virus 1,ORF D-335
ref|NP_039779.1|,Sulfolobus spindle-shaped virus 1,ORF E-54
ref|NP_039780.1|,Sulfolobus spindle-shaped virus 1,ORF F-92
ref|NP_039781.1|,Sulfolobus spindle-shaped virus 1,ORF D-244
ref|NP_039782.1|,Sulfolobus spindle-shaped virus 1,ORF E-178
ref|NP_039783.1|,Sulfolobus spindle-shaped virus 1,ORF F-93
ref|NP_039784.1|,Sulfolobus spindle-shaped virus 1,ORF E-51
ref|NP_039785.1|,Sulfolobus spindle-shaped virus 1,ORF E-96
```

(Multiple keywords must be separated using ";":)

Note: These keywords have no effect on the calculations. However, vConTACT will aggregate these gene keywords and in 
the output files you can find how many keywords were found in the same VC or PC. So, for example, if there are 5 phage 
tail fiber proteins in the same PC, that could be indicative of something pretty interesting! Or perhaps (more 
likely) there are novel genomes in your analysis. 4 of them have keywords associated to tail fibers, 1 of them is from 
a novel virus - well, if they're in the same PC then there's a good chance that "unknown" ORF is a tail fiber. vConTACT 
can't replace existing tools for understanding protein homology or similar functions, but it can help guide you through 
looking at viruses from their VC point-of-view.

```
protein_id,contig_id,keywords
ref|NP_039777.1|,Sulfolobus spindle-shaped virus 1,Fuselloviridae;Alphafusellovirus;Sulfolobus spindle-shaped virus 1;ORF B-251
ref|NP_039778.1|,Sulfolobus spindle-shaped virus 1,Fuselloviridae;Alphafusellovirus;Sulfolobus spindle-shaped virus 1;ORF D-335
ref|NP_039779.1|,Sulfolobus spindle-shaped virus 1,Fuselloviridae;Alphafusellovirus;Sulfolobus spindle-shaped virus 1;ORF E-54
ref|NP_039780.1|,Sulfolobus spindle-shaped virus 1,Fuselloviridae;Alphafusellovirus;Sulfolobus spindle-shaped virus 1;ORF F-92
ref|NP_039781.1|,Sulfolobus spindle-shaped virus 1,Fuselloviridae;Alphafusellovirus;Sulfolobus spindle-shaped virus 1;ORF D-244
ref|NP_039782.1|,Sulfolobus spindle-shaped virus 1,Fuselloviridae;Alphafusellovirus;Sulfolobus spindle-shaped virus 1;ORF E-178
ref|NP_039783.1|,Sulfolobus spindle-shaped virus 1,Fuselloviridae;Alphafusellovirus;Sulfolobus spindle-shaped virus 1;ORF F-93
ref|NP_039784.1|,Sulfolobus spindle-shaped virus 1,Fuselloviridae;Alphafusellovirus;Sulfolobus spindle-shaped virus 1;ORF E-51
ref|NP_039785.1|,Sulfolobus spindle-shaped virus 1,Fuselloviridae;Alphafusellovirus;Sulfolobus spindle-shaped virus 1;ORF E-96
```

And the run command:

```
vcontact --raw-proteins [proteins file] --rel-mode ‘Diamond’ --proteins-fp [gene-to-genome mapping file] --db 'ProkaryoticViralRefSeq85-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin [path to ClusterONE] --output-dir [target output directory]
```

##### If starting with a BLASTP or Diamond results file

In addition to the gene-to-genome mapping file (above), users must provide a tab-delimited (i.e. "tabular") BLASTP 
(-outfmt 6) or Diamond file (--outfmt 0).

```
NP_039777.1	NP_039777.1	100.0	251	0	0	1	251	1	251	2.1e-144	510.8
NP_039777.1	YP_003331457.1	49.2	238	114	4	5	239	1	234	1.4e-55	215.7
NP_039777.1	YP_003331489.1	49.6	234	111	4	2	232	20	249	3.2e-55	214.5
NP_039777.1	NP_944455.1	48.7	228	111	3	8	232	3	227	9.2e-55	213.0
NP_039777.1	YP_001552190.1	48.5	227	111	3	9	232	4	227	7.8e-54	209.9
NP_039777.1	YP_077262.1	24.6	211	131	10	29	230	89	280	6.5e-08	57.4
NP_039777.1	YP_007348313.1	24.2	211	132	10	29	230	89	280	7.1e-07	53.9
NP_944455.1	NP_039777.1	48.7	228	111	3	3	227	8	232	7.2e-54	209.9
NP_963931.1	NP_039777.1	33.8	219	133	4	17	233	16	224	5.9e-30	130.6
NP_963972.1	NP_039777.1	44.5	236	116	5	5	230	12	242	8.8e-44	176.4
```

Note that the example above could come from either Diamond or BLASTP!

```
vcontact --blast-fp [BLASTP/Diamond file] --rel-mode ‘Diamond’ --proteins-fp [gene-to-genome mapping file] --db 'ProkaryoticViralRefSeq85-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin [path to ClusterONE] --output-dir [target output directory]
```

##### If starting with contig, PC and PC profile info files

Existing vConTACT users will recognize these are the output files from vConTACT-PCs, the tool that parsed BLASTP output
 files and a gene-to-genome mapping file and generated these 3 files. Although vConTACT-PCs has been fully integrated 
 (with improved options), we still want to be able to support analyses performed with the original vConTACT (and if 
 users wish to compare v1 and v2).

**pcs.csv**: File with information about each PC. The size of the PC, how many ORFs/genes were annotated, and counts 
these annotations for the "keywords" column.

```
pc_id,size,annotated,keys
PC_00000,138,138.0,"UvsW RNA-DNA and DNA-DNA helicase ATPase (3); UvsW RNA-DNA and DNA-DNA helicase (2); RNA-DNA and DNA-DNA helicase (12); hypothetical protein Aes012_171 (1); hypothetical protein Aes508_160 (1); hypothetical protein CC2_292 (1); UvsW.1 conserved hypothetical protein (9); unnamed protein product (3); hypothetical protein PHG25ORF166w (1); uvsW.1 hypothetical protein (1); hypothetical protein ST44RRORF175w (1); DNA helicase (16)"
PC_00002,115,115.0,"gp60plus39 DNA topoisomerase subunit (8); DNA topoisomerase subunit (10); topoisomerase II large subunit (20); DNA topoisomerase II large subunit (18); unnamed protein product (2); gp39plus60 DNA topoisomerase II large subunit (2); gp60plus39 (1); putative DNA topoisomerase II (1); Putative phage DNA topoisomerase (large subunit) (2); DNA topoisomerase large subunit (6); DNA topoisomerase (1)"
PC_00003,113,113.0,Phage_cluster_17 (3); HNH_3 (1); AP2 (2); hypothetical protein JJJB_0065 (1); putative endonuclease (3); HNH homing endonuclease (6); predicted homing endonuclease (1); putative homing endonuclease (4); putative HNH endonuclease (11); hypothetical protein ABY59_0200052 (1); endonuclease (4); HNH endonuclease (15); homing endonuclease (5); putative homing endonuclease RB16 3 (1); hypothetical protein (10); gp51 (1)
PC_00004,111,111.0,Phage_cluster_43 (1); hypothetical protein PBI_ABROGATE_510 (1); hypothetical protein PBI_ABROGATE_520 (1); hypothetical protein AENEAS_54 (1); hypothetical protein AENEAS_55 (1); hypothetical protein PBI_ALSFRO_57 (1); hypothetical protein PBI_ALSFRO_56 (1); hypothetical protein ALVIN_53 (1); hypothetical protein ALVIN_52 (1); hypothetical protein ANUBIS_59 (1); hypothetical protein CKC_55 (1); hypothetical protein BARRIGA_53 (1)
PC_00005,105,105.0,"NrdA aerobic NDP reductase large subunit (6); aerobic NDP reductase large subunit (15); aerobic ribonucleoside diphosphate reductase large subunit (2); NrdA ribonucleotide reductase A subunit (4); unnamed protein product (1); NrdA aerobic ribonucleoside diphosphate reductase large subunit (1); NrdA-B aerobic NDP reductase large subunit (1); NrdA-A aerobic NDP reductase large subunit (1); ribonucleotide reductase of class Ia (aerobic) alpha subunit (11)"
PC_00006,102,102.0,DexA exonuclease A (16); exonuclease A (44); DNA exonuclease A (4); unnamed protein product (2); dexA exonuclease A (1); exonuclease (11); putative exonuclease A (2); Putative exonuclease A (1); DexA (1); dexA gene product (1); hypothetical protein Ea357_064 (1); hypothetical protein ECML134_011 (1); hypothetical protein MX01_12 (1); hypothetical protein QL01_13 (1); hypothetical protein WG01_13 (1); putative DexA exonuclease A (2)"
PC_00007,101,101.0,"gp41 replication and recombination DNA helicase (7); DNA primase-helicase subunit (27); gp41 DNA primase-helicase subunit (8); replication and recombination DNA helicase (12); unnamed protein product (2); gp41 DNA helicase (1); 41 helicase (2); DNA primase/helicase (11); DNA helicase (7); putative DNA primase-helicase subunit (4); Putative DNA helicase (1); 41 gene product (1); hypothetical protein ECML134_039 (1)"
PC_00008,101,101.0,"RnlB RNA ligase 2 (15); RNA ligase 2 (41); putative RnlB RNA ligase 2 (1); putative RNA ligase 2 (3); unnamed protein product (1); RnlB-B RNA ligase 2 (1); RnlB-A RNA ligase 2 (1); RNA ligase (21); rnlB gene product (1); RnlB 2nd RNA ligase (1); hypothetical protein ECML134_173 (1); hypothetical protein HY03_0045 (1); hypothetical protein HY03_0044 (1); putative RNA ligase (2); phage-associated RNA ligase (1)"
```

Note: parentheses ("") for rows that include ","

**profiles.csv**: Each ORF gets assigned to a PC (unless it's a singleton, in which case it's empty) and that ORF 
"position" inherits the PC it was assigned.

```
contig_id,pc_id
Sulfolobus spindle-shaped virus 1,PC_06169
Sulfolobus spindle-shaped virus 1,PC_19100
Sulfolobus spindle-shaped virus 1,PC_07015
Sulfolobus spindle-shaped virus 1,PC_08048
Sulfolobus spindle-shaped virus 1,PC_06170
Sulfolobus spindle-shaped virus 1,PC_05061
Sulfolobus spindle-shaped virus 1,PC_05065
Sulfolobus spindle-shaped virus 1,PC_06171
```

**contigs.csv**: How many proteins are associated with each contig/genome.

```
contig_id,proteins
Sulfolobus spindle-shaped virus 1,31
Sulfolobus spindle-shaped virus 2,34
Sulfolobus spindle-shaped virus 4,34
Sulfolobus spindle-shaped virus 5,34
Sulfolobus spindle-shaped virus 6,33
Sulfolobus spindle-shaped virus 7,33
Sulfolobus turreted icosahedral virus 1,36
Sulfolobus turreted icosahedral virus 2,34
```

```
vcontact --contigs-fp [contig csv file] --pcs-fp [PCs csv file] --pcprofiles-fp [PC profile csv file] --vcs-mode ClusterONE --c1-bin [path to ClusterONE] --output-dir [target output directory]
```

Note: using this method will disallow the use of reference databases!

### Example files

Example files are provided in the test_data/ directory. To use vConTACT2 with them, run the following command:

```
vcontact2 --raw-proteins test_data/VIRSorter_viral_prots.faa --rel-mode ‘Diamond’ --proteins-fp test_data/proteins.csv --db 'ProkaryoticViralRefSeq85-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin [path to ClusterONE] --output-dir VirSorted_Outputs
```

You should find a large assortment of input, intermediary and final output files in the "VirSorted_Outputs" directory. 
*Most important* are **viral_cluster_overview.csv** and **genome_by_genome_overview.csv** file. They contain a list of 
VC-by-VC and genome-by-genome processed, as well as any reference databases used. 

## Output files

There are a lot of output files generated by vConTACT2, *most of these are temporary or intermediate files that are not 
useful to the general user*. The most important files are the network and annotation files.

**genome_by_genome_overview.csv**: Contains all the taxonomic information to *reference genomes*, as well as all the 
clustering information (initial VC (VC_22), refined VC (VC_22_1)), confidence metrics, and misc scores.

One important note is that the taxonomic information *is not included* for user sequences. This means that each user 
will need to find their genome(s) of interest and check to see if reference genomes are located in the same VC. If the 
user genome is within the same VC subcluster as a reference genome, then there's a very high probability that the user 
genome is part of the same genus. If the user genome is in the same VC *but not the same subcluster as a reference*, 
then it's highly likely the two genomes are related at roughly genus-subfamily level. If there are no reference genomes 
in the same VC or VC subcluster, then it's likely that they are not related at the genus level at all. That said (more 
below), it is possible they could be related at a higher taxonomic level (subfamily, family, order)

**c1.ntw**: Contains source / target / edge weight information for all genome pairs higher than the significance 
threshold as determined by the probability that those two genomes would share N genes. The lowest value in this file 
must be greater than the minimum significance threshold (default: 1). To create a network figure in Gephi or Cytoscape, 
the user will want to import *this file* into their favorite program.

Once imported, a user can add an "annotation file" - that will be the genome_by_genome file. The annotation information 
for each genome will be added to each node/genome in the network. Then, the user can color the network figure by any 
attribute in the annotation file. This will allow the user to "quickly" create a colorful image suitable for publication.

This is documented in the protocols.io protocol, available at https://dx.doi.org/10.17504/protocols.io.x5xfq7n

Note: Many times a user will notice that their genome is connected to another (possibly reference) genome in the 
network but those two genomes won't be in the same VC subcluster or even the same VC. This doesn't mean that they 
aren't related, it just means they did not share a sufficiently significant proportion of their genes to be of the same 
genus. *They could very much be related at the subfamily or family level.* However, that's for the researcher to 
decide.

## Citation

If you find vConTACT2 useful, please cite our newly published article in Nature Biotech!:

Bin Jang, H., Bolduc, B., Zablocki, O., Kuhn, J. H., Roux, S., Adriaenssens, E. M., … Sullivan, M. B. (2019). Taxonomic assignment of uncultivated prokaryotic virus genomes is enabled by gene-sharing networks. Nature Biotechnology. https://doi.org/10.1038/s41587-019-0100-8

The theory and original code this work was based on:

Bolduc, B., Jang, H. Bin, Doulcier, G., You, Z., Roux, S., & Sullivan, M. B. (2017). vConTACT: an iVirus tool to classify double-stranded DNA viruses that infect Archaea and Bacteria. PeerJ, 5, e3243. https://doi.org/10.7717/peerj.3243

## Known Issues

* The resume function is still a bit buggy. In the case of a job failure, please start *with a fresh* output directory.
