#!/usr/bin/env python

"""
The logic behind this script has changed more times than there are numbers in the universe. Most recently, the thinking
is that Viral RefSeq should be parsed, then incorporate NCBI's virus summary table information (through SQL-like inner
join), then go through and use NCBI's taxonomy dump information, resolving any issues between the 3 (up to this point).
Finally, add ICTV taxonomic information. The problem is that every stage needs significant if's and for-loops, which
slows down this entire process. And let's not mention synonyms, misspellings, or discrepancies even within NCBI! The
final iteration takes a lazier and maybe (?) brute-force approach that loads all the information from all the sources,
then processes Viral RefSeq, querying each database for information. To resolve synonyms, a new info database was
created that basically builds all possible variations for every name, so when initially searching the databases - all
the variations are searched. This is because I don't know which name will agree with which name in which database.
ICTV could have the 'correct' naming, but ICTV's summary info uses a different name for the same virus, and ANOTHER name
for its taxonomy.
"""

import os
import sys
import gzip
import re
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
from tqdm import tqdm
from pathlib import Path
from pprint import pprint

from ete3 import NCBITaxa

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 10)
pd.set_option('display.width', 10000)

# Load NCBI taxonomy
ncbi = NCBITaxa()

# ncbi.update_taxonomy_database()  # Ensure current

version = '211'
refseq_dir = Path(f'/Users/bolduc.10/Research/Databases/ViralRefSeq-v{version}/')


def build_viral_refseq(refseq_fp_list):

    # Build NCBI Viral Refseq V85
    refseq_dict = {}

    for refseq_fp in refseq_fp_list:
        print(f'Processing {refseq_fp}')

        with gzip.open(refseq_fp, 'rt') as refseq_fh:  # MUST BE rt, NOT rb, r or any other r*
            for record in SeqIO.parse(refseq_fh, 'fasta'):

                # NCBI taxdmp provides taxID and a number of synonyms, mispellings, etc
                accession = record.id
                rdesc = record.description
                rseq = str(record.seq)
                rprotein = rdesc.split(' ', 1)[-1].split('[')[0].strip()

                virus_name = re.findall(r'\[([^]]*)\]', rdesc)[-1]

                if '[' in virus_name:
                    virus_name = virus_name.split(' - ')[0].split('-[')[0].split(' [')[0]
                    if '[' in virus_name:
                        print(virus_name)

                refseq_dict[accession] = {
                    'Virus Name': virus_name,
                    'Accession': accession,
                    'Seq': rseq,
                    'Protein': rprotein,
                    'Description': rdesc,
                }

    print('Building dataframe from viral refseq')
    return pd.DataFrame.from_dict(refseq_dict, orient='index')


refseq_fps = [
    os.path.join(refseq_dir, 'viral.1.protein.faa.gz'),
    os.path.join(refseq_dir, 'viral.2.protein.faa.gz'),
    os.path.join(refseq_dir, 'viral.3.protein.faa.gz'),
    os.path.join(refseq_dir, 'viral.4.protein.faa.gz')
]
# RefSeq has all the possible genomes/genes to which we have access
refseq_df = build_viral_refseq(refseq_fps)


def get_tax_id(species):
    """to get data from ncbi taxomomy, we need to have the taxid. we can
    get that by passing the species name to esearch, which will return
    the tax id"""
    search = Entrez.esearch(term=species, db="taxonomy", retmode="xml")
    record = Entrez.read(search)

    return record['IdList'][0]


def get_tax_data(taxid):
    """once we have the taxid, we can fetch the record"""
    search = Entrez.efetch(id=taxid, db="taxonomy", retmode="xml")
    return Entrez.read(search)


Entrez.email = "bolduc.10@osu.edu"

# Don't need to "clean up" refseq dataframe, it's not really a system burden
refseq_refs_df = refseq_df.drop_duplicates(subset=['Virus Name'], keep='first').copy()

unknown = []
ranks = ['Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Subfamily', 'Genus', 'Species']
for index, series in tqdm(refseq_refs_df.iterrows(), total=len(refseq_refs_df), unit=' genomes'):
    virus_name = series['Virus Name']

    taxID = ncbi.get_name_translator([virus_name])

    if len(taxID) > 1:
        sys.stderr.write('More than 1 tax id found')
        sys.exit(1)
    elif len(taxID) == 0:

        try:
            taxID = get_tax_id(virus_name)  # Backup
        except IndexError or len(taxID) == 0:
            unknown.append(virus_name)
            taxID = '10239'

    else:
        # print('Unable to find taxID for {} initially, but found {}'.format(organism, taxID))
        taxID = list(taxID.values())[0][0]
        refseq_refs_df.loc[index, 'TaxID'] = str(taxID)

    if isinstance(taxID, np.ndarray):  # Else int
        taxID = str(taxID[0])

    try:
        lineage = ncbi.get_lineage(taxID)  # Had 1 [0] there... tax_id.values()[0][0]
    except IndexError:
        print('IndexError')
        lineage = [1, 10239]

    lineage_names = ncbi.get_taxid_translator(lineage)  # Dictionary with taxID: name ... does not have taxon LEVEL

    for lineage_num, lineage_name in lineage_names.items():
        rank_name_dict = ncbi.get_rank([lineage_num])  # Dictionary with taxID: rank name

        for lineage_n, rank_name in rank_name_dict.items():
            if rank_name in [ele.lower() for ele in ranks]:  # != 'no rank':
                refseq_refs_df.loc[index, 'NCBI-{}'.format(rank_name)] = lineage_name

for unknown_ in unknown:
    print(f"Unable to identify the tax ID for: {unknown_}")


def build_ncbi_virus_report(report_fp):
    # Load virus table - NCBI's version with taxonomic and genomic information

    # The problem with this table is that - while it has host information - it is missing viruses/prophage that are in
    # refseq

    print('Reading NCBI viruses information file {}'.format(report_fp))
    report_df = pd.read_csv(report_fp, delimiter='\t', header=0,
                            dtype={'TaxID': str, 'BioProject ID': str})  # 7450 unique organism/name, 313620 proteins
    # Organism/Name, TaxID, BioProject, Accession, BioProject ID, Group (strand), SubGroup (family?), Size (Kb), GC%,
    # Host, Segmemts, Genes, Proteins, Release Date, Modify Date, Status
    report_df.rename(columns={'#Organism/Name': 'Organism/Name'}, inplace=True)  # Fix annoying comment

    return report_df


# NCBI virus report file with summary information on a subset of viruses, especially host
ncbi_report_fp = refseq_dir / 'viruses.txt'  # Organism
ncbi_report_df = build_ncbi_virus_report(ncbi_report_fp)  # Organism/Name

# Minor MANUAL adjustments - ETE3 will warn if necessary

ncbi_report_df = ncbi_report_df.drop_duplicates(subset=['Organism/Name'], keep='first')

# Regardless of how complete the NCBI virus report is - if this script can't parse the virus name or get a TaxID from
# the amino acid sequences, then there's no point. A virus w/out proteins/amino acid sequences is discarded
merged_report_df = refseq_refs_df.merge(ncbi_report_df, how='left', on='TaxID')

# Host, Report-Name, Report-Family

# After merge between NCBI viral refseq proteins and NCBI virus report ON TaxID, either "Virus Name" or "Species" will
# need to match ICTV "Species" report. Since everything is based on refseq proteins, MIGHT AS WELL just incorporate
# virushostdb BEFORE ICTV taxonomy..

virus_host_db_fp = os.path.join(refseq_dir, 'virushostdb.tsv')
virus_host_db_df = pd.read_csv(virus_host_db_fp, delimiter='\t', header=0, dtype={'virus tax id': str,
                                                                                  'host tax id': str})
virus_host_db_df = virus_host_db_df.drop_duplicates(subset=['virus tax id'], keep='first')

merged_report_df = merged_report_df.merge(virus_host_db_df, how='left', left_on='TaxID', right_on='virus tax id')

merged_report_df = merged_report_df.drop(columns=['Seq'])  # It goes to 1 protein from many, no need to keep


def build_ictv_report(ictv_report):
    print(f'\nReading ICTV taxonomy file {ictv_report}')
    # Sort, Order, Family, Subfamily, Genus, Species, Type Species? Exemplar Accession Number, Exemplar Isolate,
    # Genome Composition, Last Change, MSL of Last Change, Proposal, Taxon History URL
    report_df = pd.read_csv(ictv_report_fp, delimiter=',', header=0, quotechar='"', dtype=object)

    # Remove all the "sub" from ranks
    report_df = report_df.drop(columns=[col for col in report_df.columns.tolist() if 'Sub' in col])
    report_df = report_df.dropna(axis='columns', how='all')

    return report_df


# The LAST part of this merger is adding in ICTV taxonomy
# Load ICTV taxonomy
ictv_report_fp = refseq_dir / 'ICTV_Master_Species_List_2021.v1.csv'
ictv_report_df = build_ictv_report(ictv_report_fp)

# Clean up
ictv_report_df = ictv_report_df.drop_duplicates(subset=['Species'], keep='first')
# Don't want to confuse the two taxonomies
ictv_report_df.columns = [f'ICTV-{col}' if col in ranks else col for col in ictv_report_df.columns]

# 'Virus Name' = NCBI viral refseq proteins
# 'Organism/Name = NCBI virus report
# 'NCBI-Species' = ETE3 of TaxID... with TaxID coming from BOTH/MERGED above 2

# Find which name is used by ICTV, then create column that uses that as a join key
merged_report_df['join name'] = np.nan
virus_names = merged_report_df['Virus Name'].tolist()
organism_names = merged_report_df['Organism/Name'].tolist()
ete_names = merged_report_df['NCBI-species'].tolist()

for idx, series in tqdm(ictv_report_df.iterrows(), unit=' ICTV names', total=len(ictv_report_df)):
    ictv_name = series['ICTV-Species']

    if ictv_name in virus_names:
        merged_report_df.loc[merged_report_df['Virus Name'] == ictv_name, 'join name'] = ictv_name
    elif ictv_name in organism_names:
        merged_report_df.loc[merged_report_df['Organism/Name'] == ictv_name, 'join name'] = ictv_name
    elif ictv_name in ete_names:
        merged_report_df.loc[merged_report_df['NCBI-species'] == ictv_name, 'join name'] = ictv_name

# Try merge on NCBI name(s), then ETE3 species name
merged_report_df = merged_report_df.merge(ictv_report_df, how='left', left_on='join name', right_on='ICTV-Species')
merged_report_df['host'] = merged_report_df['host lineage'].apply(lambda x: x.split(';')[0] if not pd.isnull(x) else np.nan)

# Filter to remove Eukaryota - this doesn't really affect anything downstream, other than having a smaller table
# Might as well rename
merged_report_df['Host'] = merged_report_df['Host'].apply(lambda x: x.capitalize() if not pd.isnull(x) else np.nan)  # vs .replace ??

merged_report_df = merged_report_df[(merged_report_df['host'].isin(['Bacteria', 'Archaea'])) |
                                    (merged_report_df['Host'].isin(['Bacteria', 'Archaea']))]

merged_report_fp = refseq_dir / 'merged_report.csv'
merged_report_df.to_csv(merged_report_fp, sep=',', quotechar='"')

for domain in ['Archaea', 'prokaryotes']:

    host_container = domain
    if host_container == 'prokaryotes':
        host_container = 'Bacteria|Archaea'

    domain_df = merged_report_df[(merged_report_df['Host'].str.contains(host_container, na=False)) |
                                 (merged_report_df['host'].str.contains(host_container, na=False))].copy()

    # WRITE OUT ['Organism/Name', 'origin', 'order', 'family', 'genus'] for TAXONOMY
    taxonomy_columns = ['origin', 'Organism/Name'] + [rank.lower() for rank in ranks]
    taxonomy_df = pd.DataFrame(index=domain_df['Virus Name'].tolist(),
                               columns=taxonomy_columns)

    taxonomy_df['origin'] = f'RefSeq-v{version}'
    taxonomy_df['Organism/Name'] = domain_df['Virus Name'].tolist()

    def better_col(ref_comp, other_comp):

        # If anything in reference is there, it supersedes anything else
        if pd.notnull(ref_comp):
            if ref_comp != 'unassigned':
                return ref_comp  # If something's there and it's not unassigned...

        if pd.notnull(other_comp):
            if other_comp != 'unassigned':
                return other_comp

        return np.nan

    ictv_ranks = [f'ICTV-{rank}' for rank in ranks]  # No need for .capitalize
    ncbi_ranks = [f'NCBI-{rank.lower()}' for rank in ranks]

    for ictv_rank, ncbi_rank in zip(ictv_ranks, ncbi_ranks):
        if (ictv_rank in domain_df.columns.tolist()) and (ncbi_rank in domain_df.columns.tolist()):
            values = domain_df.apply(lambda x: better_col(x[ictv_rank], x[ncbi_rank]), axis=1).values
            taxonomy_df[ncbi_rank.split('-')[-1].lower()] = values  # set entire df of that rank

    # Clean up empty columns
    taxonomy_df = taxonomy_df.dropna(axis='columns', how='all')

    # Don't need species...
    taxonomy_df = taxonomy_df.drop(columns=['species'])

    # It's been years since ICTV updated their naming scheme and - since then - NCBI has updated to reflect that
    taxonomy_fp = os.path.join(refseq_dir, f'ViralRefSeq-{domain.lower()}-v{version}.Merged-reference.csv')
    taxonomy_df.to_csv(taxonomy_fp, index=False)

    # Can now build the proteins and gene-to-genome files
    domain_refseq = refseq_df[refseq_df['Virus Name'].isin(taxonomy_df.index.tolist())]

    protein_records = []
    gene2genome_d = {}
    for i, series in domain_refseq.iterrows():

        seq = series['Seq']
        genome = series['Virus Name']
        accession = series['Accession']
        protein = series['Protein']
        description = f'{protein} [{genome}]'

        record = SeqRecord(Seq(seq), id=accession, name=accession, description=description)

        protein_records.append(record)

        gene2genome_d[accession] = {
            'protein_id': accession,
            'contig_id': genome,
            'keywords': protein
        }

    print(f'There are {len(gene2genome_d)} proteins in the dictionary, {len(protein_records)} to write')

    viral_refseq_faa_fp = os.path.join(refseq_dir, f'ViralRefSeq-{domain}-v{version}.faa')
    with open(viral_refseq_faa_fp, 'w') as viral_refseq_faa_fh:
        SeqIO.write(protein_records, viral_refseq_faa_fh, 'fasta')

    # WRITE OUT GENE-TO-CONTIGS FILE - for both database AND benchmarking inputs
    protein2contig_df = pd.DataFrame.from_dict(gene2genome_d, orient='index')
    protein2contig_df = protein2contig_df.sort_values(by='contig_id')
    protein2contig_fp = refseq_dir / f'ViralRefSeq-{domain}-v{version}.protein2contig.csv'
    protein2contig_df.to_csv(protein2contig_fp, index=False)

print('Program Complete')
