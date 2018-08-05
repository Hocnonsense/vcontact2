"""General options"""
import os
import subprocess 
import logging

COLORS = {"refseq":"#B31B1B",
          "tara":"#32C6A6",
          "Myoviridae": "#4daf4a",
          "Siphoviridae" :"#984ea3",
          "Podoviridae":"#ff7f00",
          "red":"#e41a1c",
          "blue":"#377eb8",
          "green":"#4daf4a",
          "violet":"#984ea3",
          "orange":"#ff7f00",
          "yellow":"#ffff33",
          "brown":"#a65628",
          "pink":"#f781bf",
          "grey":"#999999",
          "Non affiliated":"#999999",
          "Non-affiliated":"#999999",
          "Inoviridae":"#f781bf",}

KEYWORDS =  ['ATPase', 'wedge', 'junction', 'assembly',  'activator',  'baseplate',  'capsid',  'chaperone',  'diphosphate',  'endonuclease',  'exonuclease',  'fiber',  'head',  'helicase',  'helix',  'homing',  'hydrolase',  'inhibitor',  'injection',  'integrase',  'kinase',  'ligase',  'lysis',  'lysozyme',  'membrane',  'methylase',  'methyltransferase',  'neck',  'nuclease',  'polymerase',  'portal',  'plate',  'scaffold',  'primase',  'prohead',  'protease',  'recombination',  'recombinase',  'transposase',  'reductase',  'repair',  'regulator',  'replication',  'repressor',  'ribonucleoside',  'ribonucleotide',  'structural',  'synthase',  'tail',  'tube',  'tRNA',  'terminase',  'transcription',  'transcriptional',  'transferase',  'virion']
PELAGIPHAGES = ["NC_020481","NC_020482","NC_020483","NC_020484"]
PELAGIPHAGES_NAMES = ["HTVC010P", "HTVC011P", "HTVC019P", "HTVC008M"]

