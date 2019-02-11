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
from Bio.Alphabet import IUPAC
from Bio import Entrez

from ete3 import NCBITaxa

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 10)
pd.set_option('display.width', 10000)

from pprint import pprint

# Load NCBI taxonomy
ncbi = NCBITaxa()

db_dir = os.path.abspath(os.path.join(os.getcwd(), '..', 'Databases'))
# refseq_dir = os.path.join(db_dir, 'ViralRefSeq/v85/')
refseq_dir = os.path.join('/Users/bolduc.10/Research/Databases/ViralRefSeq-v88/')


def build_viral_refseq(refseq_fp_list):

    # Build NCBI Viral Refseq V85
    refseq_dict = {}

    for refseq_fp in refseq_fp_list:
        print('Processing {}'.format(refseq_fp))

        with gzip.open(refseq_fp, 'rt') as refseq_fh:  # MUST BE rt, NOT rb, r or any other r*
            for record in SeqIO.parse(refseq_fh, 'fasta'):

                # NCBI taxdmp provides taxID and a number of synonyms, mispellings, etc

                accession = record.id
                rdesc = record.description
                rseq = str(record.seq)
                rprotein = rdesc.split(' ', 1)[-1].split('[')[0].strip()

                rorganism = re.findall(r'\[([^]]*)\]', rdesc)[-1]

                if '[' in rorganism:
                    rorganism = rorganism.split(' - ')[0].split('-[')[0].split(' [')[0]
                    if '[' in rorganism:
                        print(rorganism)

                refseq_dict[accession] = {
                    'Organism': rorganism,
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
    # os.path.join(refseq_dir, 'Campylobacter_virus_IBB35.faa.tar.gz'), # ORFs spread across 5 contigs
    # os.path.join(refseq_dir, 'Enterobacteria_phage_HX01.faa.tar.gz')  # Only 62 of 269 genes in refseq?!
]
# RefSeq has all the possible genomes/genes to which we have access
refseq_df = build_viral_refseq(refseq_fps)

# Only (???) phage mistaken naming in RefSeq? This phage name exists in proteins, but its genome name cant be found?!?!
# refseq_df['Organism'] = refseq_df['Organism'].str.replace('Enterobacteria phage phiK', 'Escherichia virus phiK')


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
lite_refseq_df = refseq_df.loc[:, ['Organism', 'Accession']].copy()
lite_refseq_df.drop_duplicates(subset=['Organism'], keep='first', inplace=True)

total_organisms = len(lite_refseq_df)
for n, (index, series) in enumerate(lite_refseq_df.iterrows(), start=1):
    organism = series['Organism']
    sys.stdout.write('\rProcessing {} ({} of {})'.format(organism, n, total_organisms))

    taxID = ncbi.get_name_translator([organism])

    if len(taxID) > 1:
        sys.stderr.write('More than 1 tax id found')
        sys.exit(1)
    elif len(taxID) == 0:

        try:
            taxID = get_tax_id(organism)  # Backup
        except IndexError or len(taxID) == 0:
            print('Unable to identify taxID for {}'.format(organism))
            taxID = '10239'

    else:
        # print('Unable to find taxID for {} initially, but found {}'.format(organism, taxID))
        taxID = list(taxID.values())[0][0]
        lite_refseq_df.loc[index, 'TaxID'] = str(taxID)

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
            if rank_name not in ['no rank', 'phylum', 'class', 'subspecies', 'subgenus']:  # != 'no rank':
                lite_refseq_df.loc[index, 'NCBI-{}'.format(rank_name)] = lineage_name


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
virus_report_fp = os.path.join(refseq_dir, 'viruses.txt')  # Organism
virus_report_df = build_ncbi_virus_report(virus_report_fp)  # Organism/Name

# Minor MANUAL adjustments
taxid_adjustments = {
    '745177': '166921',  # 745177 merged to 166921 on Nov 24, 2017
    '867696': '1987993'  # 867696 points to 1987993 on NCBI taxonomy page
    }
virus_report_df['TaxID'] = virus_report_df['TaxID'].replace(taxid_adjustments)

lite_virus_report_df = virus_report_df.loc[:, ['Organism/Name', 'Group', 'SubGroup', 'Host']].copy()
lite_virus_report_df.drop_duplicates(subset=['Organism/Name'], keep='first', inplace=True)
lite_virus_report_df.set_index('Organism/Name', inplace=True)

# Merge does work here, but need to "double" merge as RefSeq and NCBI do not always agree on taxonomy
for n, (index, series) in enumerate(lite_refseq_df.iterrows(), start=1):
    organism = series['Organism']
    ncbi_name = series['NCBI-species']
    sys.stdout.write('\rProcessing {}/{} ({} of {})'.format(organism, ncbi_name, n, total_organisms))

    try:
        lite_refseq_df.loc[index, 'Host'] = lite_virus_report_df.loc[organism, 'Host']
        lite_refseq_df.loc[index, 'Report-Name'] = organism
        lite_refseq_df.loc[index, 'Report-Family'] = lite_virus_report_df.loc[organism, 'SubGroup']
    except KeyError:  # Organism name doesnt exist
        try:
            if not pd.isnull(ncbi_name):
                lite_refseq_df.loc[index, 'Host'] = lite_virus_report_df.loc[ncbi_name, 'Host']
                lite_refseq_df.loc[index, 'Report-Name'] = ncbi_name
                lite_refseq_df.loc[index, 'Report-Family'] = lite_virus_report_df.loc[ncbi_name, 'SubGroup']
        except KeyError:
            continue


def build_ncbi_taxonomy(names_fp):
    print('Processing taxonomy information from {}'.format(names_fp))

    # Really NCBI?
    names_df = pd.read_csv(names_fp, delimiter='\t\|\t', header=None, engine='python',
                           names=['tax_id', 'name_txt', 'unique name', 'unique class'])
    names_df['tax_id'] = names_df['tax_id'].astype('str')

    # Really NCBI... X2!
    names_df['unique class'] = names_df['unique class'].str.replace('\t|', '')
    names_df['unique class'] = names_df['unique class'].str.replace('|', '')

    names_taxid_adjustments = {
        '1914209': '197310'  # 1914209 points to nothing on NCBI taxonomy page
    }
    names_df['tax_id'] = names_df['tax_id'].replace(names_taxid_adjustments)

    names_df.set_index('name_txt', inplace=True)

    return names_df


# Read NCBI taxonomy file
# This helps with taxonomic identification, but more importantly, gives SYNONYMS
# names_fi = os.path.join(refseq_dir, 'names.dmp')
# names_df = build_ncbi_taxonomy(names_fi)

def build_ictv_report(ictv_report):
    print('\nReading ICTV taxonomy file {}'.format(ictv_report))
    # Sort, Order, Family, Subfamily, Genus, Species, Type Species? Exemplar Accession Number, Exemplar Isolate,
    # Genome Composition, Last Change, MSL of Last Change, Proposal, Taxon History URL
    report_df = pd.read_csv(ictv_report_fp, delimiter=',', header=0, quotechar='"', dtype=object)

    return report_df


# Load ICTV taxonomy
ictv_report_fp = os.path.join(refseq_dir, 'ICTV_Master_Species_List_2018a_v1.csv')
ictv_report_df = build_ictv_report(ictv_report_fp)

lite_ictv_report_df = ictv_report_df.loc[:, ['Order', 'Family', 'Subfamily', 'Genus', 'Species']].copy()
lite_ictv_report_df.drop_duplicates(subset=['Species'], keep='first', inplace=True)
lite_ictv_report_df.set_index('Species', inplace=True)

for n, (index, series) in enumerate(lite_refseq_df.iterrows(), start=1):
    organism = series['Organism']
    ncbi_name = series['NCBI-species']
    sys.stdout.write('\rProcessing {}/{} ({} of {})'.format(organism, ncbi_name, n, total_organisms))

    try:
        lite_refseq_df.loc[index, 'ICTV-Species'] = organism
        for level in lite_ictv_report_df.columns:
            lite_refseq_df.loc[index, 'ICTV-{}'.format(level)] = lite_ictv_report_df.loc[organism, level]
    except KeyError:  # Organism name doesnt exist
        try:
            if not pd.isnull(ncbi_name):
                lite_refseq_df.loc[index, 'ICTV-Species'] = ncbi_name
                for level in lite_ictv_report_df.columns:
                    lite_refseq_df.loc[index, 'ICTV-{}'.format(level)] = lite_ictv_report_df.loc[ncbi_name, level]
        except KeyError:
            continue

virus_host_db_fp = os.path.join(refseq_dir, 'virushostdb.tsv')
virus_host_db_df = pd.read_csv(virus_host_db_fp, delimiter='\t', header=0, dtype={'virus tax id': str,
                                                                                  'host tax id': str})
virus_host_db_df.drop_duplicates(subset=['virus tax id'], keep='first', inplace=True)
virus_host_db_df.set_index('virus tax id', inplace=True)

for n, (index, series) in enumerate(lite_refseq_df.iterrows(), start=1):
    organism = series['Organism']
    ncbi_name = series['NCBI-species']
    taxID = series['TaxID']
    sys.stdout.write('\rProcessing {}/{} ({} of {})'.format(organism, ncbi_name, n, total_organisms))

    host = series['Host']
    if (not pd.isnull(taxID)) and (taxID in virus_host_db_df.index):  # (pd.isnull(host)) - doesn't matter, can filter
        try:
            lite_refseq_df.loc[index, 'Virus-Host-DB Host'] = virus_host_db_df.loc[taxID, 'host lineage'].split(';')[0]
        except AttributeError:
            continue

lite_refseq_fp = os.path.join(refseq_dir, 'lite_refseq_df.csv')
lite_refseq_df.to_csv(lite_refseq_fp, sep=',', quotechar='"')

# Gather genera that are found in non-prokaryotes - helps when filtering hosts for non-ICTV-examined
# UPDATE After adding the Virus-Host DB, only 41 genomes from RefSeq do not have EITHER a RefSeq host or V-H DB host
# Of those 41, 29 actually have a TaxID (12 do not). And of ALL 41, only 2 are phage - Salmonella phage c341 and
# Enterobacteria phage phiK
# non_prokaryotes = ['algae', 'diatom', 'environment', 'eukaryotes', 'fungi', 'human', 'invertebrates', 'plants',
#                    'vertebrates', 'plants', 'protozoa', '-']

# non_prokaryote_df = lite_refseq_df[lite_refseq_df['Host'].str.contains('|'.join(non_prokaryotes), na=False)]
# non_prokaryote_genera = set(non_prokaryote_df['ICTV-Genus'].unique().tolist() + non_prokaryote_df['NCBI-genus'].unique().tolist())

prokaryote_df = lite_refseq_df[(lite_refseq_df['Host'].str.contains('bacteria|archaea', na=False)) | (lite_refseq_df['Virus-Host-DB Host'].str.contains('Bacteria|Archaea', na=False))].copy()
# Remove those that we KNOW are not prokaryotic
# prokaryote_df = prokaryote_df[~(prokaryote_df['ICTV-Genus'].isin(non_prokaryote_genera) & prokaryote_df['NCBI-genus'].isin(non_prokaryote_genera))]

# Haloarcula californiae icosahedral virus 1 changed from 'environment' to 'archaea' as host

# Primary annotation
# db_dir = '/Users/bolduc.10/Research/Bioinformatics/Repositories/vcontact2/Databases'
# fst_ref = os.path.join(db_dir, '2017_ICTV_update_HB.csv')  # Name
# fst_ref_df = pd.read_csv(fst_ref, header=0, delimiter=',', engine='python')  # C is having issues with file
# fst_ref_df['Name'] = fst_ref_df['Name'].str.replace('_', ' ')

# A byproduct of speed and not parsing NCBI's entire synonyms database, it misses 34 of the 2K+ viruses
aka = {
    'Acinetobacter phage vB AbaP CEB1': 'Acinetobacter phage vB_AbaP_CEB1',
    'Actinomyces phage Av-1': 'Actinomyces virus Av1',
    'Bacillus phage BMBtpLA': 'Bacillus phage vB_BtS_BMBtp3',
    'Bacillus phage G': 'Bacillus virus G',
    'Bacillus virus Gamma': 'Bacillus phage Gamma',  # !
    'Bacillus virus GIL16c': 'Bacillus phage GIL16c',  # !
    'Bordetella phage BPP-1': 'Bordetella virus BPP1',
    'Cronobacter phage ESP2949-1': 'Cronobacter virus ESP29491',
    'Endosymbiont phage APSE-1': 'Hamiltonella virus APSE1',
    'Enterobacteria phage JS98': 'Escherichia phage JS98',
    'Enterobacteria phage P22': 'Salmonella virus P22',
    'Enterobacteria phage PRD1': 'Salmonella virus PRD1',
    'Enterobacteria phage PsP3': 'Salmonella virus PsP3',
    'Enterococcus phage vB IME195': 'Enterococcus phage vB_EfaP_IME195',
    'Enterococcus phage vB IME196': 'Enterococcus phage vB_EfaS_IME196',
    'Enterococcus phage vB IME197': 'Enterococcus phage vB_EfaS_IME197',
    'Enterococcus phage vB IME198': 'Enterococcus phage vB_EfaS_IME198',  # Exists, but "_" vB_IME198
    'Escherichia phage Jk06': 'Escherichia virus KP26',
    'Klebsiella phage vB KP1': 'Klebsiella phage vB_Kp1',
    'Mesorhizobium phagevB MloP Lo5R7ANS': 'Mesorhizobium phage vB_MloP_Lo5R7ANS',
    'Mycobacterium phage Lockley': 'Mycobacterium virus Lockley',
    'Mycobacterium phage Rumpelstiltskin': 'Mycobacterium virus Rumpelstiltskin',
    'Pseudoalteromonas Phage H103': 'Pseudoalteromonas phage H103',
    'Pseudomonas phage PA1phi': 'Pseudomonas virus PA1KOR',
    'Proteus phage pPM 1': 'Proteus phage pPM_01',
    'Rhodococcus phage E3 ': 'Rhodococcus phage E3',
    ' Rhodococcus phage E3': 'Rhodococcus phage E3',
    'Streptococcus pyogenes phage 315.1': 'Streptococcus phage 315.1',
    'Streptococcus pyogenes phage 315.2': 'Streptococcus phage 315.2',
    'Streptococcus pyogenes phage 315.3': 'Streptococcus phage 315.3',
    'Streptococcus pyogenes phage 315.4': 'Streptococcus phage 315.4',
    'Streptococcus pyogenes phage 315.5': 'Streptococcus phage 315.5',
    'Streptococcus pyogenes phage 315.6': 'Streptococcus phage 315.6',
    'uncultured phage crAssphage': 'uncultured crAssphage'
}
# Acinetobacter phage vB_AbaP_CEB1
# Campylobacter phage IBB35
# Chlamydia phage 3 (2 and 4 are in there?!)
# Salmonella phage g341c

reverse_replacer = {
    'Acinetobacter phage IME AB3': 'Acinetobacter phage IME_AB3',
 'Acinetobacter phage vB AbaM Acibel004': 'Acinetobacter phage vB_AbaM_Acibel004',
 'Acinetobacter phage vB AbaM IME200': 'Acinetobacter phage vB_AbaM_IME200',
 'Acinetobacter phage vB AbaM phiAbaA1': 'Acinetobacter phage vB_AbaM_phiAbaA1',
 'Acinetobacter phage vB AbaP Acibel007': 'Acinetobacter phage vB_AbaP_Acibel007',
 'Acinetobacter phage vB AbaP PD-6A3': 'Acinetobacter phage vB_AbaP_PD-6A3',
 'Acinetobacter phage vB AbaP PD-AB9': 'Acinetobacter phage vB_AbaP_PD-AB9',
 'Acinetobacter phage vB AbaS TRS1': 'Acinetobacter phage vB_AbaS_TRS1',
 'Aeromonas phage vB AsaM-56': 'Aeromonas phage vB_AsaM-56',
 'Alteromonas phage vB AmaP AD45-P1': 'Alteromonas phage vB_AmaP_AD45-P1',
 'Arthrobacter phage vB ArS-ArV2': 'Arthrobacter phage vB_ArS-ArV2',
 'Arthrobacter phage vB ArtM-ArV1': 'Arthrobacter phage vB_ArtM-ArV1',
 'Bacillus phage vB BanS-Tsamsa': 'Bacillus phage vB_BanS-Tsamsa',
 'Bacillus phage vB BceM Bc431v3': 'Bacillus phage vB_BceM_Bc431v3',
 'Bacillus phage vB BhaS-171': 'Bacillus phage vB_BhaS-171',
 'Bacillus phage vB BtS BMBtp3': 'Bacillus phage vB_BtS_BMBtp3',
 'Citrobacter phage vB CfrM CfP1': 'Citrobacter phage vB_CfrM_CfP1',
 'Clostridium phage vB CpeS-CP51': 'Clostridium phage vB_CpeS-CP51',
 'Cronobacter phage vB CsaM GAP161': 'Cronobacter phage vB_CsaM_GAP161',
 'Cronobacter phage vB CsaM GAP31': 'Cronobacter phage vB_CsaM_GAP31',
 'Cronobacter phage vB CsaM GAP32': 'Cronobacter phage vB_CsaM_GAP32',
 'Cronobacter phage vB CsaP GAP52': 'Cronobacter phage vB_CsaP_GAP52',
 'Cronobacter phage vB CskP GAP227': 'Cronobacter phage vB_CskP_GAP227',
 'Enterobacteria phage UAB Phi20': 'Enterobacteria phage UAB_Phi20',
 'Enterobacteria phage UAB Phi78': 'Enterobacteria phage UAB_Phi78',
 'Enterobacteria phage VT2phi 272': 'Enterobacteria phage VT2phi_272',
 'Enterobacteria phage vB EcoM VR5': 'Enterobacteria phage vB_EcoM_VR5',
 'Enterobacteria phage vB EcoP ACG-C91': 'Enterobacteria phage vB_EcoP_ACG-C91',
 'Enterobacteria phage vB EcoS NBD2': 'Enterobacteria phage vB_EcoS_NBD2',
 'Enterobacteria phage vB EcoS Rogue1': 'Enterobacteria phage vB_EcoS_Rogue1',
 'Enterobacteria phage vB KleM-RaK2': 'Enterobacteria phage vB_KleM-RaK2',
 'Enterobacteriaphage UAB Phi87': 'Enterobacteriaphage UAB_Phi87',
 'Enterococcus phage IME EF3': 'Enterococcus phage IME_EF3',
 'Enterococcus phage vB EfaP IME195': 'Enterococcus phage vB_EfaP_IME195',
 'Enterococcus phage vB EfaS IME196': 'Enterococcus phage vB_EfaS_IME196',
 'Enterococcus phage vB EfaS IME197': 'Enterococcus phage vB_EfaS_IME197',
 'Enterococcus phage vB EfaS IME198': 'Enterococcus phage vB_EfaS_IME198',
 'Enterococcus phage vB Efae230P-4': 'Enterococcus phage vB_Efae230P-4',
 'Erwinia phage vB EamM Asesino': 'Erwinia phage vB_EamM_Asesino',
 'Erwinia phage vB EamM Caitlin': 'Erwinia phage vB_EamM_Caitlin',
 'Erwinia phage vB EamM ChrisDB': 'Erwinia phage vB_EamM_ChrisDB',
 'Erwinia phage vB EamM EarlPhillipIV': 'Erwinia phage vB_EamM_EarlPhillipIV',
 'Erwinia phage vB EamM Huxley': 'Erwinia phage vB_EamM_Huxley',
 'Erwinia phage vB EamM Kwan': 'Erwinia phage vB_EamM_Kwan',
 'Erwinia phage vB EamM Phobos': 'Erwinia phage vB_EamM_Phobos',
 'Erwinia phage vB EamM-Y2': 'Erwinia phage vB_EamM-Y2',
 'Erwinia phage vB EamP Frozen': 'Erwinia phage vB_EamP_Frozen',
 'Erwinia phage vB EamP-L1': 'Erwinia phage vB_EamP-L1',
 'Erwinia phage vB EamP-S6': 'Erwinia phage vB_EamP-S6',
 'Escherichia phage 64795 ec1': 'Escherichia phage 64795_ec1',
 'Escherichia phage LM33 P1': 'Escherichia phage LM33_P1',
 'Escherichia phage vB Eco ACG-M12': 'Escherichia phage vB_Eco_ACG-M12',
 'Escherichia phage vB EcoM 112': 'Escherichia phage vB_EcoM_112',
 'Escherichia phage vB EcoM ACG-C40': 'Escherichia phage vB_EcoM_ACG-C40',
 'Escherichia phage vB EcoM AYO145A': 'Escherichia phage vB_EcoM_AYO145A',
 'Escherichia phage vB EcoM Alf5': 'Escherichia phage vB_EcoM_Alf5',
 'Escherichia phage vB EcoM ECO1230-10': 'Escherichia phage vB_EcoM_ECO1230-10',
 'Escherichia phage vB EcoM JS09': 'Escherichia phage vB_EcoM_JS09',
 'Escherichia phage vB EcoM PhAPEC2': 'Escherichia phage vB_EcoM_PhAPEC2',
 'Escherichia phage vB EcoM VR20': 'Escherichia phage vB_EcoM_VR20',
 'Escherichia phage vB EcoM VR25': 'Escherichia phage vB_EcoM_VR25',
 'Escherichia phage vB EcoM VR26': 'Escherichia phage vB_EcoM_VR26',
 'Escherichia phage vB EcoM VR7': 'Escherichia phage vB_EcoM_VR7',
 'Escherichia phage vB EcoM-UFV13': 'Escherichia phage vB_EcoM-UFV13',
 'Escherichia phage vB EcoM-VpaE1': 'Escherichia phage vB_EcoM-VpaE1',
 'Escherichia phage vB EcoM-ep3': 'Escherichia phage vB_EcoM-ep3',
 'Escherichia phage vB EcoP 24B': 'Escherichia phage vB_EcoP_24B',
 'Escherichia phage vB EcoP G7C': 'Escherichia phage vB_EcoP_G7C',
 'Escherichia phage vB EcoP GA2A': 'Escherichia phage vB_EcoP_GA2A',
 'Escherichia phage vB EcoP PhAPEC5': 'Escherichia phage vB_EcoP_PhAPEC5',
 'Escherichia phage vB EcoP PhAPEC7': 'Escherichia phage vB_EcoP_PhAPEC7',
 'Escherichia phage vB EcoP SU10': 'Escherichia phage vB_EcoP_SU10',
 'Escherichia phage vB EcoS AHP42': 'Escherichia phage vB_EcoS_AHP42',
 'Escherichia phage vB EcoS AHS24': 'Escherichia phage vB_EcoS_AHS24',
 'Escherichia phage vB EcoS AKS96': 'Escherichia phage vB_EcoS_AKS96',
 'Escherichia phage vB EcoS FFH1': 'Escherichia phage vB_EcoS_FFH1',
 'Gokushovirinae Bog1183 53': 'Gokushovirinae Bog1183_53',
 'Gokushovirinae Bog5712 52': 'Gokushovirinae Bog5712_52',
 'Gokushovirinae Bog8989 22': 'Gokushovirinae Bog8989_22',
 'Gokushovirinae Fen672 31': 'Gokushovirinae Fen672_31',
 'Gokushovirinae Fen7875 21': 'Gokushovirinae Fen7875_21',
 'Klebsiella phage vB Kp1': 'Klebsiella phage vB_Kp1',
 'Klebsiella phage vB KpnM KB57': 'Klebsiella phage vB_KpnM_KB57',
 'Klebsiella phage vB KpnM KpV477': 'Klebsiella phage vB_KpnM_KpV477',
 'Klebsiella phage vB KpnP KpV289': 'Klebsiella phage vB_KpnP_KpV289',
 'Klebsiella phage vB KpnP SU503': 'Klebsiella phage vB_KpnP_SU503',
 'Klebsiella phage vB KpnP SU552A': 'Klebsiella phage vB_KpnP_SU552A',
 'Listeria phage vB LmoM AG20': 'Listeria phage vB_LmoM_AG20',
 'Listeria phage vB LmoS 188': 'Listeria phage vB_LmoS_188',
 'Listeria phage vB LmoS 293': 'Listeria phage vB_LmoS_293',
 'Mannheimia phage vB MhM 3927AP2': 'Mannheimia phage vB_MhM_3927AP2',
 'Mannheimia phage vB MhM 587AP1': 'Mannheimia phage vB_MhM_587AP1',
 'Mannheimia phage vB MhS 1152AP2': 'Mannheimia phage vB_MhS_1152AP2',
 'Mannheimia phage vB MhS 535AP2': 'Mannheimia phage vB_MhS_535AP2',
 'Mannheimia phage vB MhS 587AP2': 'Mannheimia phage vB_MhS_587AP2',
 'Mesorhizobium phage vB MloP Lo5R7ANS': 'Mesorhizobium phage vB_MloP_Lo5R7ANS',
 'Microbacterium phage vB MoxS-ISF9': 'Microbacterium phage vB_MoxS-ISF9',
 'Microviridae Bog1249 12': 'Microviridae Bog1249_12',
 'Microviridae Bog5275 51': 'Microviridae Bog5275_51',
 'Microviridae Bog9017 22': 'Microviridae Bog9017_22',
 'Microviridae Fen2266 11': 'Microviridae Fen2266_11',
 'Microviridae Fen418 41': 'Microviridae Fen418_41',
 'Microviridae Fen4707 41': 'Microviridae Fen4707_41',
 'Microviridae Fen685 11': 'Microviridae Fen685_11',
 'Microviridae Fen7786 21': 'Microviridae Fen7786_21',
 'Microviridae Fen7895 21': 'Microviridae Fen7895_21',
 'Microviridae Fen7918 21': 'Microviridae Fen7918_21',
 'Microviridae Fen7940 21': 'Microviridae Fen7940_21',
 'Morganella phage vB MmoM MP1': 'Morganella phage vB_MmoM_MP1',
 'Morganella phage vB MmoP MP2': 'Morganella phage vB_MmoP_MP2',
 'Mycobacterium phage vB MapS FF47': 'Mycobacterium phage vB_MapS_FF47',
 'Paenibacillus phage phiIBB Pl23': 'Paenibacillus phage phiIBB_Pl23',
 'Paracoccus phage vB PmaS IMEP1': 'Paracoccus phage vB_PmaS_IMEP1',
 'Propionibacterium phage ATCC29399B C': 'Propionibacterium phage ATCC29399B_C',
 'Propionibacterium phage ATCC29399B T': 'Propionibacterium phage ATCC29399B_T',
 'Propionibacterium phage P100 1': 'Propionibacterium phage P100_1',
 'Propionibacterium phage P100 A': 'Propionibacterium phage P100_A',
 'Proteus phage pPM 01': 'Proteus phage pPM_01',
 'Proteus phage vB PmiM Pm5461': 'Proteus phage vB_PmiM_Pm5461',
 'Proteus phage vB PmiP Pm5460': 'Proteus phage vB_PmiP_Pm5460',
 'Pseudomonas phage CHA P1': 'Pseudomonas phage CHA_P1',
 'Pseudomonas phage PAK P1': 'Pseudomonas phage PAK_P1',
 'Pseudomonas phage PAK P2': 'Pseudomonas phage PAK_P2',
 'Pseudomonas phage PAK P3': 'Pseudomonas phage PAK_P3',
 'Pseudomonas phage PAK P4': 'Pseudomonas phage PAK_P4',
 'Pseudomonas phage PAK P5': 'Pseudomonas phage PAK_P5',
 'Pseudomonas phage YMC11/06/C171 PPU BP': 'Pseudomonas phage YMC11/06/C171_PPU_BP',
 'Pseudomonas phage YMC11/07/P54 PAE BP': 'Pseudomonas phage YMC11/07/P54_PAE_BP',
 'Pseudomonas phage vB Pae PS44': 'Pseudomonas phage vB_Pae_PS44',
 'Pseudomonas phage vB Pae-Kakheti25': 'Pseudomonas phage vB_Pae-Kakheti25',
 'Pseudomonas phage vB Pae-TbilisiM32': 'Pseudomonas phage vB_Pae-TbilisiM32',
 'Pseudomonas phage vB PaeM C2-10 Ab1': 'Pseudomonas phage vB_PaeM_C2-10_Ab1',
 'Pseudomonas phage vB PaeM MAG1': 'Pseudomonas phage vB_PaeM_MAG1',
 'Pseudomonas phage vB PaeM PAO1 Ab03': 'Pseudomonas phage vB_PaeM_PAO1_Ab03',
 'Pseudomonas phage vB PaeM PAO1 Ab27': 'Pseudomonas phage vB_PaeM_PAO1_Ab27',
 'Pseudomonas phage vB PaeM PS24': 'Pseudomonas phage vB_PaeM_PS24',
 'Pseudomonas phage vB PaeP C2-10 Ab09': 'Pseudomonas phage vB_PaeP_C2-10_Ab09',
 'Pseudomonas phage vB PaeP C2-10 Ab22': 'Pseudomonas phage vB_PaeP_C2-10_Ab22',
 'Pseudomonas phage vB PaeP MAG4': 'Pseudomonas phage vB_PaeP_MAG4',
 'Pseudomonas phage vB PaeP PAO1 Ab05': 'Pseudomonas phage vB_PaeP_PAO1_Ab05',
 'Pseudomonas phage vB PaeP PPA-ABTNL': 'Pseudomonas phage vB_PaeP_PPA-ABTNL',
 'Pseudomonas phage vB PaeP Tr60 Ab31': 'Pseudomonas phage vB_PaeP_Tr60_Ab31',
 'Pseudomonas phage vB PaeP p2-10 Or1': 'Pseudomonas phage vB_PaeP_p2-10_Or1',
 'Pseudomonas phage vB PaeS PAO1 Ab18': 'Pseudomonas phage vB_PaeS_PAO1_Ab18',
 'Pseudomonas phage vB PaeS PAO1 Ab30': 'Pseudomonas phage vB_PaeS_PAO1_Ab30',
 'Pseudomonas phage vB PaeS PM105': 'Pseudomonas phage vB_PaeS_PM105',
 'Pseudomonas phage vB PaeS SCH Ab26': 'Pseudomonas phage vB_PaeS_SCH_Ab26',
 'Pseudomonas phage vB PsyM KIL1': 'Pseudomonas phage vB_PsyM_KIL1',
 'Rhizobium phage vB RglS P106B': 'Rhizobium phage vB_RglS_P106B',
 'Rhizobium phage vB RleM P10VF': 'Rhizobium phage vB_RleM_P10VF',
 'Rhizobium phage vB RleM PPF1': 'Rhizobium phage vB_RleM_PPF1',
 'Rhizobium phage vB RleS L338C': 'Rhizobium phage vB_RleS_L338C',
 'Rhodovulum phage vB RhkS P1': 'Rhodovulum phage vB_RhkS_P1',
 'Salmonella phage 100268 sal2': 'Salmonella phage 100268_sal2',
 'Salmonella phage 103203 sal5': 'Salmonella phage 103203_sal5',
 'Salmonella phage 118970 sal1': 'Salmonella phage 118970_sal1',
 'Salmonella phage 118970 sal2': 'Salmonella phage 118970_sal2',
 'Salmonella phage 118970 sal3': 'Salmonella phage 118970_sal3',
 'Salmonella phage 118970 sal4': 'Salmonella phage 118970_sal4',
 'Salmonella phage 64795 sal3': 'Salmonella phage 64795_sal3',
 'Salmonella phage vB SPuM SP116': 'Salmonella phage vB_SPuM_SP116',
 'Salmonella phage vB SalM PM10': 'Salmonella phage vB_SalM_PM10',
 'Salmonella phage vB SalM SJ2': 'Salmonella phage vB_SalM_SJ2',
 'Salmonella phage vB SalM SJ3': 'Salmonella phage vB_SalM_SJ3',
 'Salmonella phage vB SemP Emek': 'Salmonella phage vB_SemP_Emek',
 'Salmonella phage vB SenMS16': 'Salmonella phage vB_SenMS16',
 'Salmonella phage vB SenS-Ent2': 'Salmonella phage vB_SenS-Ent2',
 'Salmonella phage vB SenS-Ent3': 'Salmonella phage vB_SenS-Ent3',
 'Salmonella phage vB SnwM CGG4-1': 'Salmonella phage vB_SnwM_CGG4-1',
 'Salmonella phage vB SosS Oslo': 'Salmonella phage vB_SosS_Oslo',
 'Staphylococcus phage vB SauM Remus': 'Staphylococcus phage vB_SauM_Remus',
 'Staphylococcus phage vB SauM Romulus': 'Staphylococcus phage vB_SauM_Romulus',
 'Staphylococcus phage vB SauS phi2': 'Staphylococcus phage vB_SauS_phi2',
 'Staphylococcus phage vB SepS SEP9': 'Staphylococcus phage vB_SepS_SEP9',
 'Stenotrophomonas phage vB SmaS-DLP 2': 'Stenotrophomonas phage vB_SmaS-DLP_2',
 'Synechococcus phage S-RIM2 R1 1999': 'Synechococcus phage S-RIM2 R1_1999',
 'Vibrio phage ICP2 2013 A Haiti': 'Vibrio phage ICP2_2013_A_Haiti',
 'Vibrio phage vB VchM-138': 'Vibrio phage vB_VchM-138',
 'Vibrio phage vB VpaM MAR': 'Vibrio phage vB_VpaM_MAR',
 'Vibrio phage vB VpaS MAR10': 'Vibrio phage vB_VpaS_MAR10',
 'Xanthomonas phage vB XveM DIBBI': 'Xanthomonas phage vB_XveM_DIBBI',
 'Yersinia phage vB YenM TG1': 'Yersinia phage vB_YenM_TG1',
 'Yersinia phage vB YenP AP10': 'Yersinia phage vB_YenP_AP10',
 'Yersinia phage vB YenP AP5': 'Yersinia phage vB_YenP_AP5',
 'Yersinia phage vB YenP ISAO8': 'Yersinia phage vB_YenP_ISAO8'}

replacer = {**aka, **reverse_replacer}

# for to_replace, replace_with in replacer.items():
#     fst_ref_df['Name'] = fst_ref_df['Name'].str.replace(to_replace, replace_with)
# underscored = fst_ref_df[fst_ref_df['Name'].str.contains('_')]
# Double-checked "_" items and 0 (ZERO) are not present (i.e. all that have _ in name STILL have it in the official)

# Check to see if NEW taxonomy matches the original manuscript
# AND
# counts = 0
# for n, (index, series) in enumerate(prokaryote_df.iterrows(), start=1):
#     organism = series['Organism']
#     ncbi_name = series['NCBI-species']
#     ictv_name = series['ICTV-Species']
#     sys.stdout.write('\rProcessing {}/{} ({} of {})'.format(organism, ncbi_name, n, total_organisms))
#
#     if organism in fst_ref_df['Name'].tolist():
#         fst_ref_df.loc[fst_ref_df['Name'] == organism, 'Identified'] = 'True'
#         prokaryote_df.loc[index, '2017-PeerJ'] = organism
#         prokaryote_df.loc[index, '2017-Order'] = fst_ref_df.loc[fst_ref_df['Name'] == organism, 'Order_2017'].values[0]
#         prokaryote_df.loc[index, '2017-Family'] = fst_ref_df.loc[fst_ref_df['Name'] == organism, 'Family_2017'].values[0]
#         prokaryote_df.loc[index, '2017-Genus'] = fst_ref_df.loc[fst_ref_df['Name'] == organism, 'Genus_2017'].values[0]
#         counts += 1
#     elif organism in fst_ref_df['ICTV_species_2017'].tolist():
#         fst_ref_df.loc[fst_ref_df['ICTV_species_2017'] == organism, 'Identified'] = 'True'
#         prokaryote_df.loc[index, '2017-PeerJ'] = organism
#         prokaryote_df.loc[index, '2017-Order'] = fst_ref_df.loc[fst_ref_df['ICTV_species_2017'] == organism, 'Order_2017'].values[0]
#         prokaryote_df.loc[index, '2017-Family'] = fst_ref_df.loc[fst_ref_df['ICTV_species_2017'] == organism, 'Family_2017'].values[0]
#         prokaryote_df.loc[index, '2017-Genus'] = fst_ref_df.loc[fst_ref_df['ICTV_species_2017'] == organism, 'Genus_2017'].values[0]
#         counts += 1
#
#     elif ncbi_name in fst_ref_df['ICTV_species_2017'].tolist():
#         fst_ref_df.loc[fst_ref_df['ICTV_species_2017'] == ncbi_name, 'Identified'] = 'True'
#         prokaryote_df.loc[index, '2017-PeerJ'] = ncbi_name
#         prokaryote_df.loc[index, '2017-Order'] = fst_ref_df.loc[fst_ref_df['ICTV_species_2017'] == ncbi_name, 'Order_2017'].values[0]
#         prokaryote_df.loc[index, '2017-Family'] = fst_ref_df.loc[fst_ref_df['ICTV_species_2017'] == ncbi_name, 'Family_2017'].values[0]
#         prokaryote_df.loc[index, '2017-Genus'] = fst_ref_df.loc[fst_ref_df['ICTV_species_2017'] == ncbi_name, 'Genus_2017'].values[0]
#         counts += 1
#     elif ncbi_name in fst_ref_df['Name'].tolist():
#         fst_ref_df.loc[fst_ref_df['Name'] == ncbi_name, 'Identified'] = 'True'
#         prokaryote_df.loc[index, '2017-PeerJ'] = ncbi_name
#         prokaryote_df.loc[index, '2017-Order'] = fst_ref_df.loc[fst_ref_df['Name'] == ncbi_name, 'Order_2017'].values[0]
#         prokaryote_df.loc[index, '2017-Family'] = fst_ref_df.loc[fst_ref_df['Name'] == ncbi_name, 'Family_2017'].values[0]
#         prokaryote_df.loc[index, '2017-Genus'] = fst_ref_df.loc[fst_ref_df['Name'] == ncbi_name, 'Genus_2017'].values[0]
#         counts += 1

prokaryote_fp = os.path.join(refseq_dir, 'lite_refseq_prokaryotes_df.csv')
prokaryote_df.to_csv(prokaryote_fp, sep=',', quotechar='"')
#
# fst_ref_df.set_index('Name')
# fst_ref_out_fp = os.path.join(db_dir, '2017_ICTV_update_HB_BB.csv')
# fst_ref_df.to_csv(fst_ref_out_fp)

# final_refseq = refseq_df.merge(, left_on='Organism', right_on='Organism', how='right')
# WRITE OUT ['Organism/Name', 'origin', 'order', 'family', 'genus'] for TAXONOMY
taxonomy_dict = {}
for n, (index, series) in enumerate(prokaryote_df.iterrows(), start=1):
    organism = series['Organism']
    # Actually doesn't matter the name, all that matters is what its taxonomy is
    sys.stdout.write('\rProcessing {}/{} ({} of {})'.format(organism, ncbi_name, n, total_organisms))

    origin = 'RefSeq-88'
    order = np.nan
    family = np.nan
    genus = np.nan

    # if pd.isnull(series['2017-Order']) and (series['2017-Order'] != 'unassigned'):
    if pd.isnull(series['ICTV-Order']) and (series['ICTV-Order'] != 'unassigned'):
        if pd.isnull(series['NCBI-order']) and (series['NCBI-order'] != 'unassigned'):
            order = 'unassigned'
        else:
            order = series['NCBI-order']
    else:
        order = series['ICTV-Order']
    # else:
    #     order = series['2017-Order']

    # if pd.isnull(series['2017-Family']) and (series['2017-Family'] != 'unassigned'):
    if pd.isnull(series['ICTV-Family']) and (series['ICTV-Family'] != 'unassigned'):
        if pd.isnull(series['NCBI-family']) and (series['NCBI-family'] != 'unassigned'):
            family = 'unassigned'
        else:
            family = series['NCBI-family']
    else:
        family = series['ICTV-Family']
    # else:
    #     family = series['2017-Family']

    # if pd.isnull(series['2017-Subfamily']) and (series['2017-subfamily'] != 'unassigned'):
    if pd.isnull(series['ICTV-Subfamily']) and (series['ICTV-Subfamily'] != 'unassigned'):
        if pd.isnull(series['NCBI-subfamily']) and (series['NCBI-subfamily'] != 'unassigned'):
            subfamily = 'unassigned'
        else:
            subfamily = series['NCBI-subfamily']
    else:
        subfamily = series['ICTV-Subfamily']
    # else:
    #     subfamily = series['2017-Family']

    # if pd.isnull(series['2017-Genus']) and (series['2017-Genus'] != 'unassigned'):
    if pd.isnull(series['ICTV-Genus']) and (series['ICTV-Genus'] != 'unassigned'):
        if pd.isnull(series['NCBI-genus']) and (series['NCBI-genus'] != 'unassigned'):
            genus = 'unassigned'
        else:
            genus = series['NCBI-genus']
    else:
        genus = series['ICTV-Genus']
    # else:
    #     genus = series['2017-Genus']

    taxonomy_dict[organism] = {
        'Organism/Name': organism,
        'origin': origin,
        'order': order,
        'family': family,
        'subfamily': subfamily,
        'genus': genus
    }

taxonomy_df = pd.DataFrame.from_dict(taxonomy_dict, orient='index')
taxonomy_df.sort_values(by='Organism/Name', inplace=True)
taxonomy_df.replace('unassigned', np.nan, inplace=True)
taxonomy_fp = os.path.join(refseq_dir, 'ViralRefSeq-prokaryotes-v88.reference.csv')
taxonomy_df.to_csv(taxonomy_fp, index=False)

# Filter RefSeq proteins by prokaryotes-only
final_refseq = refseq_df[refseq_df['Organism'].isin(prokaryote_df['Organism'].tolist())]

# WRITE OUT VIRAL SEQUENCES
proteins_to_write = []
protein2contig_dict = {}
for i, series in final_refseq.iterrows():

    sequence = series['Seq']
    organism = series['Organism']
    accession = series['Accession']
    protein = series['Protein']
    description = '{} [{}]'.format(protein, organism)

    ncbi_record = SeqRecord(Seq(sequence, IUPAC.protein), id=accession, name=accession, description=description)

    proteins_to_write.append(ncbi_record)

    protein2contig_dict[accession] = {
        'protein_id': accession,
        'contig_id': organism,
        'keywords': protein
    }

print(len(protein2contig_dict))
print(len(proteins_to_write))

viral_refseq_faa_fp = os.path.join(refseq_dir, 'ViralRefSeq-prokaryotes-v88.faa')
with open(viral_refseq_faa_fp, 'w') as viral_refseq_faa_fh:
    SeqIO.write(proteins_to_write, viral_refseq_faa_fh, 'fasta')

# WRITE OUT GENE-TO-CONTIGS FILE - for both database AND benchmarking inputs
protein2contig_df = pd.DataFrame.from_dict(protein2contig_dict, orient='index')
protein2contig_df.sort_values(by='contig_id', inplace=True)
protein2contig_fp = os.path.join(refseq_dir, 'ViralRefSeq-prokaryotes-v88.protein2contig.csv')
protein2contig_df.to_csv(protein2contig_fp, index=False)

print('Program Complete')
