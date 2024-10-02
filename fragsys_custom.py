### IMPORTS ###

import os
import re
import Bio
import sys
import math
import scipy
import pickle
import shutil
import logging
import argparse
import Bio.SeqIO
import importlib
import prointvar
import statistics
import subprocess
import Bio.AlignIO
import numpy as np
import configparser
import pandas as pd
import seaborn as sns
from Bio.Seq import Seq
from Bio import SeqUtils
from Bio import pairwise2
from scipy import cluster
import varalign.alignments
import scipy.stats as stats
import varalign.align_variants
import matplotlib.pyplot as plt
from prointvar.pdbx import PDBXreader
from prointvar.pdbx import PDBXwriter
import matplotlib.patches as mpatches
from prointvar.sifts import SIFTSreader
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from prointvar.config import config as cfg
from prointvar.dssp import DSSPrunner, DSSPreader
from scipy.spatial.distance import squareform, pdist
from prointvar.fetchers import download_sifts_from_ebi
from prointvar.fetchers import download_structure_from_pdbe

### DICTIONARIES AND LISTS

# arpeggio_suffixes = [
#     "atomtypes", "bs_contacts", "contacts", "specific.sift", 
#     "sift","specific.siftmatch", "siftmatch", "specific.polarmatch",
#     "polarmatch", "ri", "rings", "ari", "amri", "amam", "residue_sifts"
# ]

pdb_clean_suffixes = ["break_residues", "breaks"]

simple_ions = [
    "ZN", "MN", "CL", "MG", "CD", "NI", "NA", "IOD", "CA", "BR", "XE"
]

acidic_ions = [
    "PO4", "ACT", "SO4", "MLI", "CIT", "ACY", "VO4"
]

non_relevant_ligs_manual = [
    "DMS", "EDO", "HOH", "TRS", "GOL", "OGA", "FMN", "PG4", "PGR",
    "MPD", "TPP", "MES", "PLP", "HYP", "CSO", "UNX", "EPE", "PEG",
    "PGE", "DOD", "SUI"
]

non_relevant = non_relevant_ligs_manual + simple_ions + acidic_ions

pdb_resnames = [
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",       ### ONE OF THESE MIGHT BE ENOUGH ###
    "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"
]

aas = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",                   ### ONE OF THESE MIGHT BE ENOUGH ###
    "THR", "TRP", "TYR", "VAL", "GLX", "GLI", "NLE", "CYC"
]

aas_1l= [
    "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V",
    "-"
]

bbone = ["N", "CA", "C", "O"]

aa_code = {
    "ALA" : 'A', "CYS" : 'C', "ASP" : 'D', "GLU" : 'E',
    "PHE" : 'F', "GLY" : 'G', "HIS" : 'H', "ILE" : 'I',
    "LYS" : 'K', "LEU" : 'L', "MET" : 'M', "ASN" : 'N',
    "PRO" : 'P', "GLN" : 'Q', "ARG" : 'R', "SER" : 'S',
    "THR" : 'T', "VAL" : 'V', "TRP" : 'W', "TYR" : 'Y',
    "PYL" : 'O', "SEC" : 'U', "HYP" : 'P', "CSO" : 'C', # WEIRD ONES
    "SUI" : 'D',
}

sample_colors = ["#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe", "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080", "#ffffff", "#000000"]
#sample_colors = list(itertools.islice(rgbs(), 200)) # new_colours

consvar_class_colours = [
    "royalblue", "green", "grey", "firebrick", "orange"
]

interaction_to_color = { # following Arpeggio's colour scheme
    'clash': '#000000',
    'covalent':'#999999',
    'vdw_clash': '#999999',
    'vdw': '#999999',
    'proximal': '#999999',
    'hbond': '#f04646',
    'weak_hbond': '#fc7600',
    'xbond': '#3977db', #halogen bond
    'ionic': '#e3e159',
    'metal_complex': '#800080',
    'aromatic': '#00ccff',
    'hydrophobic': '#006633',
    'carbonyl': '#ff007f',
    'polar': '#f04646',
    'weak_polar': '#fc7600',
}

wd = os.getcwd()

### CONFIG FILE READING AND VARIABLE SAVING

config = configparser.ConfigParser()
config.read("fragsys_config.txt")

dssp_bin = config["binaries"].get("dssp_bin")
stamp_bin = config["binaries"].get("stamp_bin")
transform_bin = config["binaries"].get("transform_bin")
clean_pdb_python_bin = config["binaries"].get("clean_pdb_python_bin")
clean_pdb_bin = config["binaries"].get("clean_pdb_bin")
arpeggio_python_bin = config["binaries"].get("arpeggio_python_bin")
arpeggio_bin = config["binaries"].get("arpeggio_bin")

gnomad_vcf = config["dbs"].get("gnomad_vcf")
swissprot_pkl = config["dbs"].get("swissprot_pkl")
swissprot_path = config["dbs"].get("swissprot")
ensembl_sqlite_path = config["dbs"].get("ensembl_sqlite")

stampdir = config["other"].get("stampdir")
lig_clust_metric = config["other"].get("lig_clust_metric")
lig_clust_method = config["other"].get("lig_clust_method")

msa_fmt = config["formats"].get("msa_fmt")

stamp_ROI = config["other"].get("stamp_roi")

#arpeggio_dist = float(config["thresholds"].get("arpeggio_dist"))            # this should not be needed anymore after we switch to pdbe-arpeggio
CONS_t_h = float(config["thresholds"].get("CONS_t_h"))                      # conservation score upper threshold to consider position highly divergent. Currently only working with Shenkin divergence score.
CONS_t_l = float(config["thresholds"].get("CONS_t_l"))                      # conservation score lower threshold to consider position highly conserved. Currenyly only working with Shenkin divergence score.
CONS_ts = [CONS_t_l, CONS_t_h]
JACKHMMER_n_it = int(config["thresholds"].get("JACKHMMER_n_it"))
lig_clust_dist = float(config["thresholds"].get("lig_clust_dist"))
MES_t = float(config["thresholds"].get("MES_t"))  
MES_sig_t = float(config["thresholds"].get("MES_sig_t"))                          # Missense Enrichment Score threshold to consider a position missense-depleted, or enriched.

### FUNCTIONS

## UTILITIES

def dump_pickle(data, pickle_out):
    """
    Dumps pickle.
    """
    with open(pickle_out, "wb") as f:
        pickle.dump(data, f)

def load_pickle(pickle_in):
    """
    Loads pickle.
    """
    with open(pickle_in, "rb") as f:
        data = pickle.load(f)
    return data

## SETUP FUNCTIONS

def setup_dirs(dirs):
    """
    Creates directories if they do not exist.
    """
    for dirr in dirs:
        if os.path.isdir(dirr):
            continue
        else:
            os.mkdir(dirr)

## STAMPING FUNCTIONS

def generate_STAMP_domains(pdbs_dir, domains_out, roi = stamp_ROI):
    """
    Genereates domains file, needed to run STAMP.
    """
    pdb_files = [file for file in os.listdir(pdbs_dir) if file.endswith(".pdb")]
    with open(domains_out, "w+") as fh:
        for pdb in pdb_files:
            pdb_root, pdb_ext = os.path.splitext(pdb)
            if pdb_ext == ".pdb":
                fh.write("{} {} {{{}}}\n".format(os.path.join(pdbs_dir, pdb), pdb_root + "_" + roi, roi))

def stamp(domains, prefix, out):
    """
    Runs STAMP using the domains file as input.
    """
    if "STAMPDIR" not in os.environ:
        os.environ["STAMPDIR"] = stampdir
    args = [
        stamp_bin, "-l", domains, "-rough", "-n",
        str(2), "-prefix", prefix, ">", out
    ]
    cmd = " ".join(args)
    exit_code = os.system(cmd)
    return cmd, exit_code

def transform(matrix):
    """
    Runs TRANSFORM to obtain set of transformed coordinates.
    """
    if "STAMPDIR" not in os.environ:
        os.environ["STAMPDIR"] = stampdir
    args = [transform_bin, "-f", matrix, "-het"]
    cmd = " ".join(args)
    exit_code = os.system(cmd)
    return cmd, exit_code

def fnames_from_domains(domains_out):
    """
    Returns a list of the pdb file names from the STAMP domains file.
    """
    with open(domains_out) as f:
        lines = f.readlines()
    fnames = []
    for line in lines:
        fnames.append(line.split()[1] + ".pdb")
    return fnames

def move_supp_files(struc_files, supp_pdbs_dir, cwd):
    """
    Moves set of supperimposed coordinate files to appropriate directory.
    """
    for file in struc_files:
        if os.path.isfile(os.path.join(cwd, file)):
            shutil.move(
                os.path.join(cwd, file),
                os.path.join(supp_pdbs_dir, file)
            )

def move_stamp_output(prefix, stamp_out_dir, cwd):
    """
    Moves STAMP output files to appropriate directory.
    """
    stamp_files = sorted([file for file in os.listdir(cwd) if prefix in file]) + ["stamp_rough.trans"]
    for file in stamp_files:
        filepath = os.path.join(cwd, file)
        if os.path.isfile(filepath):
            shutil.move(filepath, os.path.join(stamp_out_dir, file))
    out_from = os.path.join(cwd, prefix + ".out")
    out_to = os.path.join(stamp_out_dir, prefix + ".out")
    doms_from = os.path.join(cwd, prefix + ".domains")
    doms_to = os.path.join(stamp_out_dir, prefix + ".domains")
    if os.path.isfile(out_from):
        shutil.move(out_from, out_to)
    if os.path.isfile(doms_from):
        shutil.move(doms_from, doms_to)

def simplify_pdb(supp_file, simple_file, struc_fmt = "mmcif"):
    """
    Simplifies pdb file by removing all non-ATOM records.
    """
    df = PDBXreader(inputfile = supp_file).atoms(format_type = struc_fmt, excluded=())
    df_hetatm = df.query('group_PDB != "ATOM"')
    if df_hetatm.empty:
        return # no HETATM records, no simple file is written
    else:
        w = PDBXwriter(outputfile = simple_file)
        w.run(df_hetatm, format_type = struc_fmt, category = "auth")

## LIGAND FUNCTIONS

def get_lig_data(supp_pdbs_dir, ligs_df_path, struc_fmt = "mmcif"):
    """
    From a directory containing a set of structurally superimposed pdbs,
    writes a .pkl file indicating the name, chain and residue number of the
    ligand(s) of interest in every pdb.
    """
    ligs_df = pd.DataFrame([])
    supp_pdb_files = [file for file in os.listdir(supp_pdbs_dir) if file.endswith(".pdb")]
    for struc in supp_pdb_files:
        struc_path = os.path.join(supp_pdbs_dir, struc)
        df = PDBXreader(inputfile = struc_path).atoms(format_type = struc_fmt, excluded=())
        hetatm_df = df.query('group_PDB == "HETATM"')
        ligs = hetatm_df.label_comp_id.unique().tolist()
        #lois = [lig for lig in ligs if lig not in non_relevant]
        lois = ligs #currently taking all ligands
        for loi in lois:
            loi_df = hetatm_df.query('label_comp_id == @loi')
            #lois_df_un = loi_df.drop_duplicates(["label_comp_id", "label_asym_id"])[["label_comp_id", "label_asym_id", "auth_seq_id"]]
            lois_df_un = loi_df.drop_duplicates(["label_comp_id", "auth_asym_id", "auth_seq_id"])[["label_comp_id", "auth_asym_id", "auth_seq_id"]] # changing this to auth, cause it seems pdbe-arpeggio uses auth (except for _com_id).
            lois_df_un["struc_name"] = struc
            ligs_df = ligs_df.append(lois_df_un)
    ligs_df = ligs_df[["struc_name", "label_comp_id", "auth_asym_id", "auth_seq_id"]]
    ligs_df.to_pickle(ligs_df_path)
    return ligs_df

## SIFTS FUNCTIONS

def get_swissprot(): 
    """
    Retrieves sequences and their data from Swiss-Prot

    :param db: absolute path to a fasta file containing sequences, Swiss-Prot database by default
    :type db: str
    :returns: dictionary containing the sequence id, description and sequence for all proteins in Swiss-Prot
    :rtpe: dict
    """
    swissprot_dict = Bio.SeqIO.parse(swissprot_path, "fasta")
    proteins = {}
    for protein in swissprot_dict:
        acc = protein.id.split("|")[1]
        proteins[acc] = {}
        proteins[acc]["id"] = protein.id
        proteins[acc]["desc"] = protein.description
        proteins[acc]["seq"] = protein.seq
    return proteins

def retrieve_mapping_from_struc(struc, uniprot_id, struc_dir, mappings_dir, swissprot, struc_fmt = "mmcif"):
    """
    Retrieves the mapping between the UniProt sequence and the PDB sequence by doing an alignment.
    """
    input_struct = os.path.join(struc_dir, struc)
    pdb_structure = PDBXreader(inputfile = input_struct).atoms(format_type = struc_fmt, excluded=()) # ProIntVar reads the local file
    
    seq_record = str(swissprot[uniprot_id]["seq"])
    pps = pdb_structure.query('group_PDB == "ATOM"')[['label_comp_id', 'auth_asym_id', 'auth_seq_id']].drop_duplicates().groupby('auth_asym_id')  # groupby chain
    pdb_chain_seqs = [(chain, SeqUtils.seq1(''.join(seq['label_comp_id'].values)), seq['auth_seq_id'].values) for chain, seq in pps] # list of tuples like: [(chain_id, chain_seq, [chain resnums])]
    alignments = [pairwise2.align.globalxs(str(seq_record),chain_seq[1], -5, -1) for chain_seq in pdb_chain_seqs] # list of lists of tuples containing SwissProt seq - PDB chain seq pairwise alignment
    
    maps = []
    for pdb_chain_seq, alignment in zip(pdb_chain_seqs, alignments):
        PDB_UniProt_map = pd.DataFrame(
            [(i, x) for i, x in enumerate(alignment[0][1], start=1)],  # create aligned PDB sequences to dataframe
            columns=['UniProt_ResNum', 'PDB_ResName']
        )
        PDB_UniProt_map = PDB_UniProt_map.assign(UniProt_ResName = list(alignment[0][0]))
        PDB_index = PDB_UniProt_map.query('PDB_ResName != "-"').index
        PDB_UniProt_map = PDB_UniProt_map.assign(PDB_ResNum = pd.Series(pdb_chain_seq[2], index = PDB_index)) # adds PDB_ResNum column
        PDB_UniProt_map = PDB_UniProt_map.assign(PDB_ChainID = pd.Series(pdb_chain_seq[0], index = PDB_index)) # adds PDB_ChainId column
        maps.append(PDB_UniProt_map)
    prointvar_mapping = pd.concat(maps)
    prointvar_mapping = prointvar_mapping[['UniProt_ResNum','UniProt_ResName','PDB_ResName','PDB_ResNum','PDB_ChainID']]
    prointvar_mapping = prointvar_mapping[~prointvar_mapping.PDB_ResNum.isnull()]
    prointvar_mapping.PDB_ResNum = prointvar_mapping.PDB_ResNum.astype(int)
    struc_root, _ = os.path.splitext(struc) 
    prointvar_mapping_csv = os.path.join(mappings_dir, struc_root + "_mapping.csv")
    prointvar_mapping.PDB_ResNum = prointvar_mapping.PDB_ResNum.astype(str) # PDB_ResNum is a string, not an integer
    prointvar_mapping.to_csv(prointvar_mapping_csv, index = False)
    return prointvar_mapping

## DSSP FUNCTIONS

def run_dssp(struc, supp_pdbs_dir, dssp_dir):
    """
    Runs DSSP, saves and return resulting output dataframe
    """
    struc_root, _ = os.path.splitext(struc)
    dssp_csv = os.path.join(dssp_dir, struc_root + ".csv") # output csv filepath
    dssp_out = os.path.join(dssp_dir, struc_root + ".dssp")
    struc_in = os.path.join(supp_pdbs_dir, struc)
    DSSPrunner(inputfile = struc_in, outputfile = dssp_out).write()            # runs DSSP
    dssp_data = DSSPreader(inputfile = dssp_out).read()            # reads DSSP output
    dssp_data = dssp_data.rename(index = str, columns = {"RES": "PDB_ResNum"})
    dssp_data.PDB_ResNum = dssp_data.PDB_ResNum.astype(str)
    dssp_cols = ["PDB_ResNum", "SS", "ACC", "KAPPA", "ALPHA", "PHI", "PSI", "RSA"]    # selects subset of columns
    dssp_data.to_csv(dssp_csv, index = False)
    return dssp_data[dssp_cols]

## ARPEGGIO FUNCTIONS

def run_clean_pdb(pdb_path):
    """
    Runs pdb_clean.py to prepare files for Arpeggio.
    """
    args = [
        clean_pdb_python_bin, clean_pdb_bin, pdb_path
    ]
    cmd = " ".join(args)
    exit_code = os.system(cmd)
    return cmd, exit_code

def run_arpeggio(pdb_path, lig_sel, out_dir):
    """
    runs Arpeggio
    """
    args = [
        arpeggio_python_bin, arpeggio_bin, pdb_path,
        "-s", lig_sel, "-o", out_dir, "--mute"
    ]
    cmd = " ".join(args)
    exit_code = os.system(cmd)
    return exit_code, cmd

def switch_columns(df, names):
    """
    Switches columns in Arpeggio DataFrame, so that ligand atoms and protein
    atoms are always on the same column.
    """
    # Columns to switch
    columns_to_switch = [
        'auth_asym_id', 'auth_atom_id', 'auth_seq_id', 'label_comp_id'
    ]

    # Iterate through the DataFrame and switch columns where necessary
    for index, row in df.iterrows():
        if row['label_comp_id_end'] in names:
            for col in columns_to_switch:
                bgn_col = f"{col}_bgn"
                end_col = f"{col}_end"
                df.at[index, bgn_col], df.at[index, end_col] = df.at[index, end_col], df.at[index, bgn_col]

    return df

def map_values(row, pdb2up):
    """
    maps UniProt ResNums from SIFTS dictionary from CIF file to Arpeggio dataframe.
    """
    try:
        return pdb2up[row.auth_asym_id_end][row.auth_seq_id_end]
    except KeyError:
        log.debug(f'Residue {row.auth_seq_id_end} chain {row.auth_asym_id_end} has no mapping to UniProt')
        return np.nan # if there is no mapping, return NaN

def create_resnum_mapping_dict(df):
    """
    Creates a dictionary with the mapping between PDB_ResNum and UniProt_ResNum
    froma  previously created dataframe.
    """
    chain_mapping = {}
    
    # Iterate over the dataframe
    for index, row in df.iterrows():
        chain_id = row['PDB_ChainID']
        pdb_resnum = row['PDB_ResNum']
        uniprot_resnum = row['UniProt_ResNum']
        
        # Initialize dictionary for the chain ID if it doesn't exist
        if chain_id not in chain_mapping:
            chain_mapping[chain_id] = {}
        
        # Add PDB_ResNum as key and UniProt_ResNum as value for the current chain
        chain_mapping[chain_id][str(pdb_resnum)] = uniprot_resnum
    
    return chain_mapping

def process_arpeggio_df(arp_df, pdb_id, ligand_names, pdb2up):
    """
    Process Arpeggio Df to put in appropriate
    format to extract fingerprings. Also filter out
    non-relevant interactions.
    """
    
    arp_df_end_expanded = arp_df['end'].apply(pd.Series)
    arp_df_bgn_expanded = arp_df['bgn'].apply(pd.Series)

    arp_df = arp_df.join(arp_df_end_expanded).drop(labels='end', axis=1)
    arp_df = arp_df.join(arp_df_bgn_expanded, lsuffix = "_end", rsuffix = "_bgn").drop(labels='bgn', axis = 1)

    arp_df.auth_seq_id_bgn = arp_df.auth_seq_id_bgn.astype(str)
    arp_df.auth_seq_id_end = arp_df.auth_seq_id_end.astype(str)

    inter_df = arp_df.query('interacting_entities == "INTER" & type == "atom-atom"').copy().reset_index(drop = True)

    inter_df = inter_df[inter_df['contact'].apply(lambda x: 'clash' not in x)].copy().reset_index(drop = True) # filtering out clashes
    
    inter_df = inter_df.query('label_comp_id_bgn in @pdb_resnames or label_comp_id_end in @pdb_resnames').copy().reset_index(drop = True) # filtering out ligand-ligand interactions
    
    if inter_df.empty:
        log.warning("No protein-ligand interaction  for {}".format(pdb_id))
        return inter_df, "no-PL-inters"
    
    inter_df = inter_df.query('label_comp_id_bgn in @ligand_names or label_comp_id_end in @ligand_names').copy().reset_index(drop = True) # filtering out non-LOI interactions (only to avoid re-running Arpeggio, once it has been run with wrong selection)
    
    switched_df = switch_columns(inter_df, ligand_names)

    switched_df = switched_df.query('label_comp_id_end in @pdb_resnames').copy() # filtering out non-protein-ligand interactions

    # Apply the function and create a new column
    switched_df["UniProt_ResNum_end"] = switched_df.apply(lambda row: map_values(row, pdb2up), axis=1)
    
    switched_df = switched_df.sort_values(by=["auth_asym_id_end", "UniProt_ResNum_end", "auth_atom_id_end"]).reset_index(drop = True)
    
    return switched_df, "OK"

def process_arpeggio(struc, all_ligs, clean_pdbs_dir, arpeggio_dir, mappings_dir, struc_fmt = "mmcif"): 
    """
    Processes arpeggio output to generate the two tables that will
    be used later in the analysis.
    """
    lig_cons_splits = []
    arpeggio_lig_conss = []

    struc_root, _ = os.path.splitext(struc)
        
    for lig in all_ligs:
        arpeggio_out_path_lig = os.path.join(arpeggio_dir, "{}.clean_{}.bs_contacts".format(struc_root, lig))
        
        all_cons_lig = pd.read_table(arpeggio_out_path_lig, header = None)

        lig_cons_split_lig = reformat_arpeggio(all_cons_lig)

        lig_cons_split_lig = add_resnames_to_arpeggio_table(struc, clean_pdbs_dir, lig_cons_split_lig, struc_fmt)

        lig_cons_split_lig = ligand_to_atom2(lig_cons_split_lig, lig)
        
        lig_cons_split_lig["contact_type"] = lig_cons_split_lig.apply(lambda row: contact_type(row), axis = 1)

        arpeggio_lig_cons_lig = lig_cons_split_lig.sort_values(by = ["Chain (Atom1)", "ResNum (Atom1)"])
        arpeggio_lig_cons_lig = arpeggio_lig_cons_lig[["ResNum (Atom1)","Chain (Atom1)", 'ResName (Atom1)']] 
        arpeggio_lig_cons_lig = arpeggio_lig_cons_lig.drop_duplicates(subset = ["ResNum (Atom1)", "Chain (Atom1)"])
        arpeggio_lig_cons_lig = arpeggio_lig_cons_lig.rename(index = str, columns = {"ResNum (Atom1)": "PDB_ResNum", "Chain (Atom1)": "PDB_ChainID"}) 
        arpeggio_lig_cons_lig = arpeggio_lig_cons_lig.astype({"PDB_ResNum": int})
        
        lig_cons_splits.append(lig_cons_split_lig)
        arpeggio_lig_conss.append(arpeggio_lig_cons_lig)
    
    lig_cons_split = pd.concat(lig_cons_splits)
    arpeggio_lig_cons = pd.concat(arpeggio_lig_conss)
    
    ########################### ADDED TO AVOID INCORRECT BS DEFINITION DUE TO PDB RESNUMS ###########################

    struc_mapping = pd.read_csv(os.path.join(mappings_dir, "{}.mapping".format(struc_root)))

    mapping_dict = {}
    for chain, chain_df in struc_mapping.groupby("PDB_ChainID"):
        for i, row in chain_df.iterrows():
            mapping_dict[(chain, str(row["PDB_ResNum"]))] = row["UniProt_ResNum"]
    
    lig_cons_split['ResNum (Atom1)'] = lig_cons_split['ResNum (Atom1)'].astype(str)
    arpeggio_lig_cons['PDB_ResNum'] = arpeggio_lig_cons['PDB_ResNum'].astype(str)
    lig_cons_split["UniProt_Resnum"] = lig_cons_split.set_index(['Chain (Atom1)', 'ResNum (Atom1)']).index.map(mapping_dict.get)
    arpeggio_lig_cons["UniProt_Resnum"] = arpeggio_lig_cons.set_index(['PDB_ChainID', 'PDB_ResNum']).index.map(mapping_dict.get)

    ########################### ADDED TO AVOID INCORRECT BS DEFINITION DUE TO PDB RESNUMS ###########################

    old_len1 = len(arpeggio_lig_cons)
    old_len2 = len(lig_cons_split)
    lig_cons_split = lig_cons_split[~lig_cons_split.UniProt_Resnum.isnull()]
    arpeggio_lig_cons = arpeggio_lig_cons[~arpeggio_lig_cons.UniProt_Resnum.isnull()]
    new_len1 = len(arpeggio_lig_cons)
    new_len2 = len(lig_cons_split)
    if new_len1 != old_len1:
        log.warning("{} residues lacked mapping from PDB to UniProt at arpeggio_lig_cons for {}".format(new_len1 - old_len1, struc_root))
    if new_len2 != old_len2:
        log.warning("{} residues lacked mapping from PDB to UniProt at lig_cons_split for {}".format(new_len2 - old_len2, struc_root))

    ########################### ADDED TO AVOID INCORRECT BS DEFINITION DUE TO LACKING MAPPING TO PDB RESNUMS ###########################

    lig_cons_split.to_csv(os.path.join(arpeggio_dir,  "arpeggio_all_cons_split_{}.csv".format(struc_root)), index = False)
    arpeggio_lig_cons.to_csv(os.path.join(arpeggio_dir,  "arpeggio_lig_cons_{}.csv".format(struc_root)), index = False)
    
    return lig_cons_split, arpeggio_lig_cons

def reformat_arpeggio(arpeggio_df):
    """
    Starts formatting arpeggio table.
    """
    arpeggio_df.columns = [
        'Atom_1', 'Atom_2', 'Clash', 'Covalent', 'VdW Clash', 'Vdw', 'Proximal', 'Hydrogen Bond',
        'Weak Hydrogen Bond', 'Halogen bond',  'Ionic', 'Metal Complex', 'Aromatic', 'Hydrophobic',
        'Carbonyl', 'Polar', 'Weak Polar', 'Atom proximity', 'Vdw proximity', 'Interacting entities'
    ]
    lig_cons = arpeggio_df.loc[arpeggio_df["Interacting entities"] == "INTER"] # Selecting only the interactions between our specified atoms (i.e ligands) and other selections
    lig_cons = lig_cons.sort_values(by = ["Atom_1"])
    
    split_atom1 = lig_cons.Atom_1.str.split("/", expand = True) # splits atom1 column into three new columns: chain, resnum and atom
    split_atom1.columns = ["Chain (Atom1)", "ResNum (Atom1)", "Atom (Atom1)"]
    split_atom2 = lig_cons.Atom_2.str.split("/", expand = True)  # splits atom2 column into three new columns: chain, resnum and atom
    split_atom2.columns = ["Chain (Atom2)", "ResNum (Atom2)", "Atom (Atom2)"]
    
    lig_cons_split = pd.merge(split_atom2, lig_cons, left_index = True, right_index = True) # Making a table of the contacts, but with the atom identifier split into chain, resnum, and atom
    lig_cons_split = pd.merge(split_atom1, lig_cons_split, left_index = True, right_index = True) # Making a table of the contacts, but with the atom identifier split into chain, resnum, and atom
    lig_cons_split = lig_cons_split.drop(axis = 1, labels = ["Atom_1", "Atom_2"])
    lig_cons_split["ResNum (Atom1)"] = lig_cons_split["ResNum (Atom1)"].astype(int)
    lig_cons_split["ResNum (Atom2)"] = lig_cons_split["ResNum (Atom2)"].astype(int)
    return lig_cons_split

def add_resnames_to_arpeggio_table(structure, clean_pdbs_dir, arpeggio_cons_split, struc_fmt = "mmcif"):
    """
    adds residue names to arpeggio table, needed for later table mergings
    """
    structure_root, _ = os.path.splitext(structure) 
    structure_path = os.path.join(clean_pdbs_dir, "{}.clean.pdb".format(structure_root))
    pdb_structure = PDBXreader(inputfile = structure_path).atoms(format_type = struc_fmt, excluded=())
    resnames_dict = {(row.auth_asym_id, int(row.auth_seq_id)): row.label_comp_id for index, row in pdb_structure.drop_duplicates(['auth_asym_id', 'auth_seq_id']).iterrows()}
    arpeggio_cons_split["ResName (Atom1)"] = arpeggio_cons_split.set_index(["Chain (Atom1)", "ResNum (Atom1)"]).index.map(resnames_dict.get)
    arpeggio_cons_split["ResName (Atom2)"] = arpeggio_cons_split.set_index(["Chain (Atom2)", "ResNum (Atom2)"]).index.map(resnames_dict.get)
    return arpeggio_cons_split

def ligand_to_atom2(lig_cons_split, lig):
    """
    formats arpeggio table so that the ligand atoms are always Atom2
    """
    ordered_cols = [
        'Chain (Atom1)', 'ResNum (Atom1)', 'ResName (Atom1)', 'Atom (Atom1)',
        'Chain (Atom2)', 'ResNum (Atom2)', 'ResName (Atom2)', 'Atom (Atom2)',
        'Clash', 'Covalent', 'VdW Clash', 'Vdw', 'Proximal', 'Hydrogen Bond',
        'Weak Hydrogen Bond', 'Halogen bond', 'Ionic', 'Metal Complex', 'Aromatic',
        'Hydrophobic', 'Carbonyl', 'Polar', 'Weak Polar','Atom proximity',
        'Vdw proximity', 'Interacting entities'
    ]
    lig_is_atom1 = lig_cons_split[lig_cons_split["ResName (Atom1)"] == lig].sort_values("ResNum (Atom2)")
    lig_is_atom2 = lig_cons_split[lig_cons_split["ResName (Atom2)"] == lig].sort_values("ResNum (Atom1)")

    lig_is_atom1.rename(columns = {
        'Chain (Atom1)': 'Chain (Atom2)', 'Chain (Atom2)': 'Chain (Atom1)',
        'ResNum (Atom1)': 'ResNum (Atom2)', 'ResNum (Atom2)': 'ResNum (Atom1)',
        'Atom (Atom1)': 'Atom (Atom2)', 'Atom (Atom2)': 'Atom (Atom1)',
        'ResName (Atom1)': 'ResName (Atom2)', 'ResName (Atom2)': 'ResName (Atom1)'
    }, inplace = True)

    lig_cons_split_rf = pd.concat([lig_is_atom1[ordered_cols], lig_is_atom2[ordered_cols]]) # new dataframe so ligand is always atom2
    lig_cons_split_rf = lig_cons_split_rf[lig_cons_split_rf["ResName (Atom1)"].isin(aas)]
    return lig_cons_split_rf

def contact_type(row):
    """
    determines whether a row is backbone or sidechain contact
    """
    if row["Atom (Atom1)"] in bbone:
        return "backbone"
    else:
        return "sidechain"

def generate_dictionary(mmcif_file):
    """
    Generates coordinate dictionary from a mmCIF file.
    """
    # Parse the mmCIF file
    mmcif_dict = MMCIF2Dict(mmcif_file)

    # Initialise the result dictionary
    result = {}

    # Iterate through the atoms and populate the dictionary
    for i, auth_asym_id in enumerate(mmcif_dict["_atom_site.auth_asym_id"]):
        label_comp_id_end = mmcif_dict["_atom_site.label_comp_id"][i]
        auth_seq_id = mmcif_dict["_atom_site.auth_seq_id"][i]
        auth_atom_id_end = mmcif_dict["_atom_site.auth_atom_id"][i]
        x = mmcif_dict["_atom_site.Cartn_x"][i]
        y = mmcif_dict["_atom_site.Cartn_y"][i]
        z = mmcif_dict["_atom_site.Cartn_z"][i]

        # Dictionary creation
        result[auth_asym_id, label_comp_id_end, auth_seq_id, auth_atom_id_end] = [x, y, z]

    return result

def determine_width(interactions):
    """
    Generates cylinder width for 3DMol.js interaction
    representation depending on Arpeggio contact
    fingerprint.
    """
    return 0.125 if 'vdw_clash' in interactions else 0.0625

def determine_color(interactions):
    """
    Generates cylinder colour for 3DMol.js interaction
    representation depending on Arpeggio contact
    fingerprint.
    """
    undef = ['covalent', 'vdw', 'vdw_clash', 'proximal']
    if len(interactions) == 1 and interactions[0] in undef:
        return '#999999'
    else:
        colors = [interaction_to_color[interaction] for interaction in interactions if interaction in interaction_to_color and interaction not in undef]
        if colors:
            return colors[0]
        else:
            log.critical("No color found for {}".format(interactions))
            return None  # Return the first color found, or None if no match

### FUNCTIONS FOR SITE DEFINITION

def generate_ligs_res_df(arpeggio_dir):
    """
    Given a directory containing processed Arpeggio output files,
    returns a dataset containing information about the ligands of
    interest binding the protein and their labels.
    """
    p = re.compile("""arpeggio_all_cons_split_(.+).csv""")
    lig_files = sorted([file for file in os.listdir(arpeggio_dir) if file.startswith("arpeggio_all_cons_split")])
    file_ids, ligs, resnums, chains, ligs_ress = [[], [], [], [], []]
    for file in lig_files:
        m = p.match(file)
        if m:
            file_id = m.group(1)
        else:
            log.critical("Failed getting file ID for {}".format(file))
        #comps = file.split("_")
        #file_id = "{}_{}".format(comps[-2], comps[-1].split(".")[0]) # CAREFUL WITH THIS ID
        df = pd.read_csv(os.path.join(arpeggio_dir, file))
        for lig, lig_df in df.groupby(["ResName (Atom2)", "ResNum (Atom2)", "Chain (Atom2)"]):
            lig_ress = lig_df["UniProt_Resnum"].unique().tolist()
            file_ids.append(file_id)
            ligs.append(lig[0])
            resnums.append(lig[1])
            chains.append(lig[2])
            ligs_ress.append(lig_ress)
    lig_data_df = pd.DataFrame(list(zip(file_ids, ligs, resnums, chains, ligs_ress)), columns = ["pdb_id", "lig_name", "lig_resnum", "lig_chain", "binding_res"])
    labs = [file_ids[i] + "_" + str(ligs[i]) + "_" + str(resnums[i]) + "_" + str(chains[i]) for i in range(len(ligs))]
    lig_data_df["lab"] = labs
    return lig_data_df

def get_dis_file(lig_data_df, out):
    """
    Creates .dis file to be fed to OC.
    """
    lig_res = lig_data_df.binding_res.tolist() #this is a list of lists, each list contains residue numbers interacting with ligand
    labs = lig_data_df.lab.tolist()
    intersect_dict = get_intersect_rel_matrix(lig_res)
    n_ligs = len(lig_res)
    with open(out, "w+") as fh:
        fh.write(str(n_ligs) + "\n")
        for lab in labs: #labs contain unique identifier for a ligand
            fh.write(lab + "\n")
        for i in range(n_ligs):
            for j in range(i+1, n_ligs):
                fh.write(str(intersect_dict[i][j]) + "\n")

def get_intersect_rel_matrix(binding_ress):
    """
    Given a set of ligand binding residues, calcualtes a
    similarity matrix between all the different sets of ligand
    binding residues.
    """
    inters = {i: {} for i in range(len(binding_ress))}
    for i in range(len(binding_ress)):
        inters[i][i] = intersection_rel(binding_ress[i], binding_ress[i])
        for j in range(i+1, len(binding_ress)):
            inters[i][j] = intersection_rel(binding_ress[i], binding_ress[j])
            inters[j][i] = inters[i][j]
    return inters

def intersection_rel(l1, l2):
    """
    Calculates relative intersection.
    """
    len1 = len(l1)
    len2 = len(l2)
    I_max = min([len1, len2])
    I = len(list(set(l1).intersection(l2)))
    return I/I_max

def write_bs_files(frag_mean_coords, bs_def_out, attr_out, chimera_script_out):
    """
    Writes files for binding site definition.
    """
    frag_mean_coords = frag_mean_coords.dropna()
    frag_mean_coords.binding_site = frag_mean_coords.binding_site.astype(int)
    frag_mean_coords.lig_resnum = frag_mean_coords.lig_resnum.astype(int)
    chimera_atom_spec = (
        ':'+ frag_mean_coords.lig_resnum.astype(str) +
        '.'+ frag_mean_coords.lig_chain +
        '&#/name==' + frag_mean_coords.pdb_path
        )
    frag_mean_coords = frag_mean_coords.assign(chimera_atom_spec = chimera_atom_spec)  
    frag_mean_coords.to_pickle(bs_def_out) # saves table to pickle
    write_bs_attribute_file(frag_mean_coords, attr_out)
    bs_labs = frag_mean_coords.binding_site.unique().tolist()
    write_chimera_script(chimera_script_out, bs_labs)

def write_bs_attribute_file(clustered_fragments, attr_out):
    """
    writes Chimera attribute file to later colour ligands
    according to the binding site they bind to
    """
    with open(attr_out, "w") as out:
        out.write("attribute: binding_site\n")
        out.write("match mode: 1-to-1\n")
        out.write("recipient: residues\n")
        out.write("\n".join("\t" + clustered_fragments.chimera_atom_spec.values + "\t" + clustered_fragments.binding_site.astype(str)))

def write_chimera_script(chimera_script_out, bs_labels):
    """
    writes Chimera script that will format the superimposed structures
    as well as colour ligands according to their binding site
    """
    
    #chimera_rgb_string = [','.join(map(str, rgb)) for rgb in sample_colors]
    cmds = [
        "~rib", "rib #0", "ksdssp", "set silhouette", "set silhouettewidth 3",
        "background solid white", "~dis", "sel ~@/color=white", "dis sel", "namesel lois",
        "~sel"
    ]
    with open(chimera_script_out, 'w') as out:
        out.write('# neutral colour for everything not assigned a cluster\n')
        out.write('colour white\n')
    
        out.write('# colour each binding site\n')
        for i in range(0, len(bs_labels)):
            #out.write('colour {} :/binding_site=={}\n'.format(','.join(list(map(str, list(sample_colors[i])))), i))
            out.write('colour {} :/binding_site=={}\n'.format(sample_colors[i], i))
        out.write("### SOME FORMATTING ###\n")
        out.write("\n".join(cmds))
    log.info("Chimera script successfully created")

def get_residue_bs_membership(cluster_ress):
    """
    Returns a dictionary indicating to which ligand binding
    site each ligand binding residue is found in. A residue
    might contribute to more than one adjacent binding site.
    """
    all_bs_ress = []
    for v in cluster_ress.values():
        all_bs_ress.extend(v)
    all_bs_ress = sorted(list(set(all_bs_ress)))
    
    bs_ress_membership_dict = {}
    for bs_res in all_bs_ress:
        bs_ress_membership_dict[bs_res] = []
        for k, v in cluster_ress.items():
            if bs_res in v:
                bs_ress_membership_dict[bs_res].append(k) # which binding site each residue belongs to
    
    return bs_ress_membership_dict

### CONSERVATION + VARIATION FUNCTIONS

def create_alignment_from_struc(example_struc, fasta_path, struc_fmt = "mmcif", n_it = JACKHMMER_n_it, seqdb = swissprot_path):
    """
    Given an example structure, creates and reformats an MSA.
    """
    create_fasta_from_seq(get_seq_from_pdb(example_struc, struc_fmt = struc_fmt), fasta_path) # CREATES FASTA FILE FROM PDB FILE
    fasta_root, _ = os.path.splitext(fasta_path)    
    hits_out = "{}.out".format(fasta_root)
    hits_aln = "{}.sto".format(fasta_root)
    hits_aln_rf = "{}_rf.sto".format(fasta_root)    
    jackhmmer(fasta_path, hits_out, hits_aln , n_it = n_it, seqdb = seqdb,) # RUNS JACKHAMMER USING AS INPUT THE SEQUENCE FROM THE PDB AND GENERATES ALIGNMENT
    add_acc2msa(hits_aln, hits_aln_rf)

def get_seq_from_pdb(pdb_path, struc_fmt = "mmcif"): # MIGHT NEED TO BE MORE SELECTIVE WITH CHAIN, ETC
    """
    Generates aa sequence string from a pdb coordinates file.
    """
    struc = PDBXreader(pdb_path).atoms(format_type = struc_fmt, excluded=())
    return "".join([aa_code[aa] for aa in struc.query('group_PDB == "ATOM"').drop_duplicates(["auth_seq_id"]).label_comp_id.tolist()])

def create_fasta_from_seq(seq, out):
    """
    Saves input sequence to fasta file to use as input for jackhmmer.
    """
    with open(out, "w+") as fh:
        fh.write(">query\n{}\n".format(seq))

def jackhmmer(seq, hits_out, hits_aln, n_it = JACKHMMER_n_it, seqdb = swissprot_path):
    """
    Runs jackhmmer on an input seq for a number of iterations and returns exit code, should be 0 if all is OK.
    """
    args = ["jackhmmer", "-N", str(n_it), "-o", hits_out, "-A", hits_aln, seq, seqdb]
    cmd = " ".join(args)
    exit_code = os.system(cmd)
    return cmd, exit_code

def add_acc2msa(aln_in, aln_out, fmt_in = msa_fmt):
    """
    Modifies AC field of jackhmmer alignment in stockholm format.
    
    :param aln_in: path of input alignment
    :type aln_in: str, required
    :param aln_out: path of output alignment
    :type aln_in: str, required
    :param fmt_in: input and output MSA format
    :type aln_in: str, defaults to stockholm
    """
    aln = Bio.SeqIO.parse(aln_in, fmt_in)
    recs = []
    for rec in aln:
        if rec.id == "query":
            continue
        else:
            rec.annotations["accession"] = rec.id.split("|")[1]
            recs.append(rec)
    Bio.SeqIO.write(recs, aln_out, fmt_in)

def get_target_prot_cols(msa_in, msa_fmt = msa_fmt): 
    """
    Returns list of MSA col idx that are popualted on the protein target.
    """
    seqs = [str(rec.seq) for rec in Bio.SeqIO.parse(msa_in, msa_fmt) if "query" in rec.id]
    occupied_cols = [i+1 for seq in seqs for i, el in enumerate(seq) if el != "-"]
    return sorted(list(set(occupied_cols)))

def calculate_shenkin(aln_in, aln_fmt, out = None):
    """
    Given an MSA, calculates Shenkin ans occupancy, gap
    percentage for all columns.
    """
    cols = in_columns(aln_in, aln_fmt)
    scores = []
    occ = []
    gaps = []
    occ_pct = []
    gaps_pct = []
    for k, v in cols.items():
        scores.append(get_shenkin(k, v))
        stats = (get_stats(v))
        occ.append(stats[0])
        gaps.append(stats[1])
        occ_pct.append(stats[2])
        gaps_pct.append(stats[3])
    df = pd.DataFrame(list(zip(list(range(1,len(scores)+1)),scores, occ,gaps, occ_pct, gaps_pct)), columns = ["col", "shenkin", "occ", "gaps", "occ_pct", "gaps_pct"])
    if out != None:
        df.to_pickle(out)#, index = False)
    return df

def get_stats(col):
    """
    for a given MSA column, calculates some basic statistics
    such as column residue occupancy ang gaps
    """
    n_seqs = len(col)
    gaps = col.count("-")
    occ = n_seqs - gaps
    occ_pct = round(100*(occ/n_seqs), 2)
    gaps_pct = round(100-occ_pct, 2)
    return occ, gaps, occ_pct, gaps_pct

def in_columns(aln_in, infmt):
    """
    Returns dictionary in which column idx are the key
    and a list containing all aas aligned to that column
    is the value.
    """
    aln = Bio.AlignIO.read(aln_in, infmt)
    n_cols = len(aln[0])
    cols = {}
    for col in range(1,n_cols+1):
        cols[col] = []
    for row in aln:
        seq = str(row.seq)
        for i in range(0, len(seq)):
            cols[i+1].append(seq[i])
    return cols

def get_shenkin(i_col, col):
    """
    Calculates Shenkin score for an MSA column.
    """
    S = get_entropy(get_freqs(i_col, col))
    return round((2**S)*6,2)

def get_freqs(i_col, col):
    """
    Calculates amino acid frequences for a given MSA column.
    """

    abs_freqs = {aa: 0 for aa in aas_1l}
    non_standard_aas = {}
    for aa in col:
        aa = aa.upper()
        if col.count("-") == len(col):
            abs_freqs["-"] = 1
            return abs_freqs
        if aa in aas_1l:
            abs_freqs[aa] += 1
        else:
            if aa not in non_standard_aas:
                non_standard_aas[aa] = 0
            non_standard_aas[aa] += 1
    all_ns_aas = sum(non_standard_aas.values())
    if all_ns_aas != 0:
        log.warning("Column {} presents non-standard AAs: {}".format(str(i_col), non_standard_aas))
    rel_freqs = {k: v/(len(col) - all_ns_aas) for k, v in abs_freqs.items()}
    return rel_freqs

def get_entropy(freqs):
    """
    Calculates Shannon's entropy from a set of aa frequencies.
    """
    S = 0
    for f in freqs.values():
        if f != 0:
            S += f*math.log2(f)
    return -S

def format_shenkin(shenkin, prot_cols, out = None):
    """
    Formats conservation dataframe and also
    calculates two normalised versions of it.
    """
    shenkin_filt = shenkin[shenkin.col.isin(prot_cols)].copy()
    shenkin_filt.index = range(1, len(shenkin_filt) + 1) # CONTAINS SHENKIN SCORE, OCCUPANCY/GAP PROPORTION OF CONSENSUS COLUMNS
    min_shenkin = min(shenkin_filt.shenkin)
    max_shenkin = max(shenkin_filt.shenkin)
    shenkin_filt.loc[:, "rel_norm_shenkin"] = round(100*(shenkin_filt.shenkin - min_shenkin)/(max_shenkin - min_shenkin), 2) # ADDING NEW COLUMNS WITH DIFFERENT NORMALISED SCORES
    shenkin_filt.loc[:, "abs_norm_shenkin"] = round(100*(shenkin_filt.shenkin - 6)/(120 - 6), 2)
    if out != None:
        shenkin_filt.to_pickle(out)#, index = False)
    return shenkin_filt

def get_human_subset_msa(aln_in, human_msa_out, fmt_in = msa_fmt):
    """
    Creates a subset MSA containing only human sequences.
    """
    msa = Bio.AlignIO.read(aln_in, fmt_in)
    human_recs = []
    for rec in msa:
        if "HUMAN" in rec.name:
            human_recs.append(rec)
    Bio.SeqIO.write(human_recs, human_msa_out, fmt_in)

def cp_sqlite(wd, og_path = ensembl_sqlite_path):
    """
    Copies ensembl_cache.sqlite
    to execution directory.
    """
    hidden_var_dir = os.path.join(wd, ".varalign")
    sqlite_name = os.path.basename(og_path)
    if not os.path.isdir(hidden_var_dir):
        os.mkdir(hidden_var_dir)
    else:
        pass
    cp_path = os.path.join(hidden_var_dir, sqlite_name)
    shutil.copy(og_path, cp_path)
    return cp_path

def rm_sqlite(cp_path):
    """
    Removes ensembl_cache.sqlite
    from execution directory.
    """
    hidden_var_dir = os.path.dirname(cp_path)
    os.remove(cp_path)
    os.rmdir(hidden_var_dir)

def format_variant_table(df, col_mask, vep_mask = ["missense_variant"], tab_format = "gnomad"):
    """
    Formats variant table, by gettint rid of empty rows that are not human sequences,
    changning column names and only keeping those variants that are of interest and
    are present in columns of interest.
    """
    df_filt = df.copy(deep = True)
    df_filt.reset_index(inplace = True)
    if tab_format == "gnomad":
        df_filt.columns = [" ".join(col).strip() for col in df_filt.columns.tolist()]
    df_filt.columns = [col.lower().replace(" ", "_") for col in df_filt.columns.tolist()]
    df_filt = df_filt[df_filt.source_id.str.contains("HUMAN")]
    df_filt = df_filt.dropna(subset = ["vep_consequence"])
    df_filt = df_filt[df_filt.vep_consequence.isin(vep_mask)]
    df_filt = df_filt[df_filt.alignment_column.isin(col_mask)]
    return df_filt

def get_missense_df(aln_in, variants_df, shenkin_aln, prot_cols, aln_out, aln_fmt = msa_fmt, get_or = True):
    """
    Generates a dataframe for the subset of human sequences with variants
    mapping to them. Calculates shenkin, and occupancy data, and then
    enrichment in variants.
    """
    variants_aln = generate_subset_aln(aln_in, aln_fmt, variants_df, aln_out)
    if variants_aln == "":
        return pd.DataFrame()
    variants_aln_info = calculate_shenkin(variants_aln, aln_fmt)
    variants_aln_info = variants_aln_info[variants_aln_info.col.isin(prot_cols)]
    vars_df = pd.DataFrame(variants_df.alignment_column.value_counts().reindex(prot_cols, fill_value = 0).sort_index()).reset_index()
    vars_df.index = range(1, len(prot_cols) + 1)
    vars_df.columns = ["col", "variants"]
    merged = pd.merge(variants_aln_info, vars_df, on = "col", how = "left")
    merged.index = range(1, len(vars_df) + 1)
    merged["shenkin"] = shenkin_aln["shenkin"]
    merged["rel_norm_shenkin"] = shenkin_aln["rel_norm_shenkin"] 
    merged["abs_norm_shenkin"] = shenkin_aln["abs_norm_shenkin"]
    if get_or == True:
        merged_or = get_OR(merged)
        return merged_or
    else:
        return merged

def generate_subset_aln(aln_in, aln_fmt, df, aln_out = None):
    """
    Creates a subset MSA containing only human sequences that present
    missense variants and returns the path of such MSA.
    """
    seqs_ids = df.source_id.unique().tolist()
    aln = Bio.SeqIO.parse(aln_in, aln_fmt)
    variant_seqs = [rec for rec in aln if rec.id in seqs_ids]
    n_variant_seqs = len(variant_seqs)
    if n_variant_seqs == 0:
        return ""
    else:
        log.info("There are {} sequences with variants for {}".format(str(n_variant_seqs), aln_in))
    if aln_out == None:
        aln_root, aln_ext = os.path.splitext(aln_in)
        aln_out =  "{}_variant_seqs{}".format(aln_root, aln_ext)
    Bio.SeqIO.write(variant_seqs, aln_out, aln_fmt)
    return aln_out

def get_OR(df, variant_col = "variants"):
    """
    Calculates OR, and associated p-value and CI,
    given a missense dataframe with variants and occupancy.
    """
    tot_occ = sum(df.occ)
    tot_vars = sum(df[variant_col])
    idx = df.index.tolist()
    for i in idx:
        i_occ = df.loc[i,"occ"]
        i_vars = df.loc[i,variant_col]
        rest_occ = tot_occ - i_occ
        rest_vars = tot_vars - i_vars
        if i_occ == 0:
            oddsr = np.nan 
            pval = np.nan
            se_or = np.nan
            log.debug("0 occupancy. Returning np.nan")
        else:
            if i_vars == 0:
                i_occ += 0.5
                i_vars += 0.5
                rest_occ += 0.5
                rest_vars += 0.5
                log.debug("0 variants. Adding 0.5 to each cell")
            oddsr, pval = stats.fisher_exact([[i_vars, rest_vars], [i_occ, rest_occ]])
            vals = [i_vars, rest_vars, i_occ, rest_occ]
            se_or = 1.96*(math.sqrt(sum(list(map((lambda x: 1/x), vals)))))
        df.loc[i, "oddsratio"] = round(oddsr, 2)
        df.loc[i, "pvalue"] = round(pval, 2)
        df.loc[i, "se_OR"] = round(se_or, 2)
    return df

def add_miss_class(df, miss_df_out = None, cons_col = "shenkin", MES_t = MES_t, cons_ts = CONS_ts, colours = consvar_class_colours):
    """
    Adds two columns to missense dataframe. These columns will put columns
    into classes according to their divergence and missense enrichment score.
    It also adds a column that will help colour MSA columns according to their
    classifications.
    """
    for i in df.index:
        if df.loc[i, cons_col] <= cons_ts[0] and df.loc[i, "oddsratio"] < MES_t:
            df.loc[i, "miss_class"] = "CMD"
        elif df.loc[i, cons_col] <= cons_ts[0] and df.loc[i, "oddsratio"] > MES_t:
            df.loc[i, "miss_class"] = "CME"
        elif df.loc[i, cons_col] >= cons_ts[1] and df.loc[i, "oddsratio"] < MES_t:
            df.loc[i, "miss_class"] = "UMD"
        elif df.loc[i, cons_col] >= cons_ts[1] and df.loc[i, "oddsratio"] > MES_t:
            df.loc[i, "miss_class"] = "UME"
        else:
            df.loc[i, "miss_class"] = "None"
    coloring = {
        "CMD": colours[0],
        "CME": colours[1],
        "UMD": colours[3],
        "UME": colours[4],
        "None": colours[2]
    }
    df["miss_color"] =  df.miss_class.map(coloring)
    
    if miss_df_out != None:
        df.to_pickle(miss_df_out)#, index = False)
    return df

def merge_shenkin_df_and_mapping(shenkin_df, mapping_df, aln_ids):
    """
    merges conservation, and variation table with MSA-UniProt
    mapping table, so conservation and variation data
    are mapped to UniProt residues.
    """
    shenkin_df = shenkin_df.rename(index = str, columns = {"col": "alignment_column"}) # renaming columns to be consistent with other StruVarPi dataframes
    prot_mapping = mapping_df.copy(deep = True).loc[aln_ids]
    prot_mapping.columns = prot_mapping.columns.droplevel(1)
    prot_mapping.reset_index(inplace = True)
    prot_mapping = prot_mapping.rename_axis(None, axis = "columns")
    prot_mapping.rename(index = None, columns = {"Alignment": "MSA_column", "Protein_position": "UniProt_ResNum"}, inplace = True)
    # Merging the VarAlign data to the Pfam alignment to gain conservation and variation data for the whole family...
    mapped_data = pd.merge(
        prot_mapping[["MSA_column", "UniProt_ResNum"]], shenkin_df,
        left_on = "MSA_column", right_on = "alignment_column"
    ).drop("MSA_column", axis = 1)
    return mapped_data

### SETTING UP LOG

logging.basicConfig(filename = "fragsys_REVAMPED.log", format = '%(asctime)s %(name)s [%(levelname)-8s] - %(message)s', level = logging.INFO)

log = logging.getLogger("FRAGSYS_REVAMPED")

### MAIN FUNCTION

def main(args):
    """
    Main function of the script. Calls all other functions.
    """

    ### PARSING ARGUMENTS

    log.info("Logging initiated")

    for arg, value in sorted(vars(args).items()):
        log.info("Argument %s: %r", arg, value)

    input_dir = args.input_dir
    uniprot_id = args.uniprot_id
    struc_fmt = args.struc_fmt
    override = args.override
    override_variants = args.override_variants
    run_variants = args.variants

    ### SETTING UP DIRECTORIES

    input_id = os.path.normpath(input_dir).split(os.sep)[-1]
    output_dir = os.path.join(wd, "output", input_id)
    results_dir = os.path.join(output_dir, "results")
    stamp_out_dir = os.path.join(output_dir, "stamp_out")
    supp_pdbs_dir = os.path.join(output_dir, "supp_pdbs")
    supp_cifs_dir = os.path.join(output_dir, "supp_cifs")
    simple_pdbs_dir = os.path.join(output_dir, "simple_pdbs")
    clean_pdbs_dir = os.path.join(output_dir, "clean_pdbs")
    pdb_clean_dir = os.path.join(output_dir, "pdb_clean")
    mappings_dir = os.path.join(output_dir, "mappings")
    dssp_dir = os.path.join(output_dir, "dssp")
    arpeggio_dir = os.path.join(output_dir, "arpeggio")
    varalign_dir = os.path.join(output_dir, "varalign")
    
    dirs = [
        output_dir, results_dir, stamp_out_dir,
        supp_pdbs_dir, simple_pdbs_dir, clean_pdbs_dir,
        pdb_clean_dir, mappings_dir, dssp_dir,
        arpeggio_dir, varalign_dir, supp_cifs_dir,
    ]
    
    setup_dirs(dirs) # creates directory /output/input_dir_name and then all results subdirectories: results, stamp_out, clean_pdbs, etc...

    log.info("Directories created")

    ### CHECKING IF FINAL RESULTS TABLE ALREADY EXISTS

    final_table_out = os.path.join(results_dir, "{}_results_table.pkl".format(input_id))
    if os.path.isfile(final_table_out) and not override:
        log.info("Final results table already exists")
        log.info("Skipping to the end")
        sys.exit(0)
    if struc_fmt == "mmcif":
        strucs = [f for f in os.listdir(input_dir) if f.endswith(".cif") and "clean" not in f] # different extension depending on structure format
    elif struc_fmt == "pdb":
        strucs = [f for f in os.listdir(input_dir) if f.endswith(".pdb") and "clean" not in f]
    n_strucs = len(strucs)
    log.info("Number of structures: {}".format(n_strucs))

    ### CLEANING FILES

    ## if structure format is mmcif, we need to convert to pdb since clean_pdb.py only works with pdb files

    for struc in strucs:
        struc_root, struc_ext = os.path.splitext(struc)  # "structure_name", ".pdb"
        pdb_path = os.path.join(input_dir, struc)        # /input_dir/structure_name.pdb
        pdb_path_root, _ =  os.path.splitext(pdb_path)   # /input_dir/structure_name
        clean_struc_path_from = f'{pdb_path_root}.clean{struc_ext}' # /input_dir/structure_name.clean.pdb
        clean_struc_path_to = os.path.join(clean_pdbs_dir, f'{struc_root}.clean{struc_ext}') # /output_dir/clean_pdbs/structure_name.clean.pdb
        if os.path.isfile(clean_struc_path_to): # if clean file already exists in clean files dir
            log.debug("PDB {} already cleaned".format(struc))
            pass
        else: # if clean file does not exist in clean files dir
            cmd, ec = run_clean_pdb(pdb_path) # run pdb_clean.py on pdb file
            if ec == 0:
                log.debug("{} cleaned".format(struc))
                if os.path.isfile(clean_struc_path_from): 
                    shutil.move(clean_struc_path_from, clean_struc_path_to) # move clean file to clean files dir
                for pdb_clean_suff in pdb_clean_suffixes: # for each suffix in pdb_clean_suffixes
                    pdb_clean_file_from = os.path.join(input_dir, "{}.{}".format(struc, pdb_clean_suff))
                    pdb_clean_file_to = os.path.join(pdb_clean_dir, "{}.{}".format(struc, pdb_clean_suff))
                    if os.path.isfile(pdb_clean_file_from):
                        shutil.move(pdb_clean_file_from, pdb_clean_file_to)
            else:
                log.critical("pdb_clean.py failed with {}".format(cmd))

    ### STAMP SECTION

    domains_out = os.path.join(results_dir, "{}_stamp.domains".format(input_id))
    if os.path.isfile(domains_out):
        log.debug("STAMP domains file already exists")
        pass
    else:
        generate_STAMP_domains(clean_pdbs_dir, domains_out)
        log.info("STAMP domains file generated")

    prefix = "{}_stamp".format(input_id)
    n_domains = len(open(domains_out).readlines())
    log.info("Number of domains: {}".format(str(n_domains)))
    matrix_file = "{}.{}".format(prefix, str(n_domains-1))
    last_matrix_path = os.path.join(stamp_out_dir, matrix_file)

    if os.path.isfile(last_matrix_path):
        log.debug("STAMP matrix files already exist")
        pass
    else:
        cmd, ec = stamp(
            domains_out,
            prefix, os.path.join(results_dir, "{}.out".format(prefix))
        )
        if ec == 0:
            log.info("STAMP matrix files generated")
        else:
            log.critical("STAMP failed with {}".format(cmd))
        
    fnames = fnames_from_domains(domains_out)

    c = 0 # counting the number of superposed pdbs
    for file in fnames:
        if os.path.isfile(os.path.join(supp_pdbs_dir, file)): # only when they alaready have been transformed
            c += 1
    if c == n_domains:
        log.debug("All structure domains already are superposed")
        pass
    else:
        if not os.path.isfile(matrix_file): # RUNNING TRANSFORM ONCE STAMP OUTPUT HAS BEEN MOVED TO STAMP_OUT_DIR
            matrix_file = os.path.join(stamp_out_dir, matrix_file) # needs to change to accommodate for how many domains in .domains file 
        cmd, ec = transform(matrix_file) #running transform with matrix on cwd
        if ec == 0:
            log.info("Structures transformed")
        else:
            log.critical("TRANSFORM failed with {}".format(cmd))
    
    log.info("STAMP and TRANSFORM completed")
    
    move_supp_files(fnames, supp_pdbs_dir, wd)

    move_stamp_output(prefix, stamp_out_dir, wd)

    simple_pdbs = [f for f in os.listdir(simple_pdbs_dir) if f.endswith(".pdb")]
    n_simple_pdbs = len(simple_pdbs) # number of simplified pdbs, will be 0 the first time thiis command is executed

    if n_simple_pdbs == n_domains:
        log.debug("Structure domains already simplified")
        pass
    else:
        supp_files = sorted([f for f in os.listdir(supp_pdbs_dir) if f.endswith(".pdb")])
        shutil.copy(os.path.join(supp_pdbs_dir, supp_files[0]), simple_pdbs_dir) # copy first chain as is
        print(f'Keeping protein atoms for {supp_files[0]}')
        for file in supp_files[1:]: # we keep protei atoms for first chain
            if file.endswith(".pdb"):
                supp_file = os.path.join(supp_pdbs_dir, file)
                simple_file = os.path.join(simple_pdbs_dir, file) 
                if os.path.isfile(simple_file):
                    log.debug("Simple pdb file already exists")
                    pass
                else:
                    #remove_extra_ligs(supp_file, simple_file)
                    simplify_pdb(supp_file, simple_file, struc_fmt)
        log.info("All structure domains have been simplified") # what we want to do here is seimply keep protein atoms for first chain, this is to make visualisation quicker and simpler
    
    log.info("PDB simplification completed")

    ### GET LIGAND DATA

    lig_data_path = os.path.join(results_dir, "{}_lig_data.pkl".format(input_id))
    if override or not os.path.isfile(lig_data_path):
        ligs_df = get_lig_data(simple_pdbs_dir, lig_data_path, struc_fmt) # this now uses auth fields as it seems that is what pdbe-arpeggio uses.
        log.info("Saved ligand data")
    else:
        ligs_df = pd.read_pickle(lig_data_path)
        log.debug("Loaded ligand data")

    ### UNIPROT MAPPING SECTION

    swissprot = load_pickle(swissprot_pkl)
    log.debug("Swissprot loaded")

    for struc in fnames: #fnames are now the files of the STAMPED PDB files, not the original ones
        struc_root, _ =  os.path.splitext(struc)
        struc_mapping_path = os.path.join(mappings_dir, "{}_mapping.csv".format(struc_root))
        struc_mapping_dict_path = os.path.join(mappings_dir, "{}_mapping.pkl".format(struc_root))
        if override or not os.path.isfile(struc_mapping_path):
            mapping = retrieve_mapping_from_struc(struc, uniprot_id, supp_pdbs_dir, mappings_dir, swissprot, struc_fmt = struc_fmt) # giving supp, here, instead of simple because we want them all
            mapping_dict = create_resnum_mapping_dict(mapping)
            dump_pickle(mapping_dict, struc_mapping_dict_path)
            log.info("Mapping files for {} generated".format(struc_root))
        else:
            mapping = pd.read_csv(struc_mapping_path)
            mapping_dict = load_pickle(struc_mapping_dict_path)
            log.debug("Mapping files for {} already exists".format(struc_root))

    log.info("UniProt mapping section completed")
            
    ### DSSP SECTION   

    dssp_mapped_out = os.path.join(results_dir, "{}_dssp_mapped.pkl".format(input_id))

    if override or not os.path.isfile(dssp_mapped_out):

        mapped_dssps = []
        for struc in fnames: #fnames are now the files of the STAMPED PDB files, not the original ones
            ## DSSP
            struc_root, _ =  os.path.splitext(struc)
            dssp_csv = os.path.join(dssp_dir, "{}.csv".format(struc_root))

            if override or not os.path.isfile(dssp_csv):
                dssp_data = run_dssp(struc, supp_pdbs_dir, dssp_dir)
                log.info("DSSP run successfully on {}".format(struc_root))
            else:
                dssp_data = pd.read_csv(dssp_csv)
                log.debug("DSSP data already exists")
            
            ## UNIPROT MAPPING

            dssp_data.PDB_ResNum = dssp_data.PDB_ResNum.astype(str) 

            mapping = pd.merge(mapping, dssp_data, left_on = "PDB_ResNum", right_on = "PDB_ResNum") # don't think this merging worked well

            mapped_dssps.append(mapping)

            mapped_dssp_df = pd.concat(mapped_dssps)

            mapped_dssp_df.to_pickle(os.path.join(results_dir, "{}_dssp_mapped.pkl".format(input_id)))
    else:
        mapped_dssp_df = pd.read_pickle(dssp_mapped_out)
        log.debug("Loaded mapped DSSP data")

    log.info("DSSP section completed")

    ### ARPEGGIO PART ###

    struc2ligs = {}
    for struc in fnames:
        struc_root, _ =  os.path.splitext(struc)
        struc2ligs[struc] = []
        struc_df = ligs_df.query('struc_name == @struc')

        pdb_path = os.path.join(supp_pdbs_dir, struc)
        pdb_path_root, _ =  os.path.splitext(pdb_path)

        pdb_df = PDBXreader(pdb_path).atoms(format_type = struc_fmt, excluded=())
        cif_out = os.path.join(supp_cifs_dir, "{}.cif".format(struc_root))

        pdb_df["label_alt_id"] = "."
        # cif_df["pdbx_PDB_ins_code"] = "?"
        pdb_df["pdbx_formal_charge"] = "?"

        cif_cols_order = [
            "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id", "label_asym_id",
            # "label_entity_id",
            "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
            "B_iso_or_equiv", "pdbx_formal_charge", "auth_seq_id", "auth_comp_id", "auth_asym_id", "auth_atom_id", "pdbx_PDB_model_num"
        ]
        w = PDBXwriter(outputfile = cif_out)
        w.run(pdb_df[cif_cols_order], format_type = "mmcif")
        
        if struc_df.empty:
            log.warning("No ligands in {}".format(struc))
            continue

        lig_sel = " ".join(["/{}/{}/".format(row.auth_asym_id, row.auth_seq_id) for _, row in  struc_df.iterrows()])

        ligand_names = list(set([row.label_comp_id for _, row in struc_df.iterrows()]))

        arpeggio_default_json_name = os.path.basename(struc).split(".")[0]
        arpeggio_default_out = os.path.join(arpeggio_dir, f"{arpeggio_default_json_name}.json")  # this is how arpeggio names the file (splits by "." and takes the first part)

        arpeggio_out = os.path.join(arpeggio_dir, struc_root + ".json")
        arpeggio_proc_df_out = os.path.join(arpeggio_dir, struc_root + "_proc.pkl")

        if override or not os.path.isfile(arpeggio_proc_df_out):

            if override or not os.path.isfile(arpeggio_out):

                ec, cmd = run_arpeggio(cif_out, lig_sel, arpeggio_dir)
                if ec != 0:
                    log.error("Arpeggio failed for {} with {}".format(struc, cmd))
                    continue

                shutil.move(arpeggio_default_out, arpeggio_out)

            arp_df = pd.read_json(arpeggio_out)

            pdb2up = load_pickle(os.path.join(mappings_dir, "{}_mapping.pkl".format(struc_root)))

            proc_inters, fp_stat = process_arpeggio_df(arp_df, struc_root, ligand_names, pdb2up)

            coords_dict = generate_dictionary(cif_out)

            proc_inters["coords_end"] = proc_inters.set_index(["auth_asym_id_end", "label_comp_id_end", "auth_seq_id_end", "auth_atom_id_end"]).index.map(coords_dict.get)
            proc_inters["coords_bgn"] = proc_inters.set_index(["auth_asym_id_bgn", "label_comp_id_bgn", "auth_seq_id_bgn", "auth_atom_id_bgn"]).index.map(coords_dict.get)

            proc_inters["width"] = proc_inters["contact"].apply(determine_width)
            proc_inters["color"] = proc_inters["contact"].apply(determine_color)
            proc_inters.to_pickle(arpeggio_proc_df_out)
            log.debug("Arpeggio processed data saved")
    else:
        proc_inters = pd.read_pickle(arpeggio_proc_df_out)
        log.debug("Loaded arpeggio processed data")

    print("Finishing here!")
    sys.exit(0)

    ligand_contact_list = []
    for struc in fnames:
        struc_root, _ =  os.path.splitext(struc)
        all_ligs = ligs_df.query('struc_name == @struc').label_comp_id.unique().tolist()
        arpeggio_out1 = os.path.join(arpeggio_dir, "arpeggio_all_cons_split_{}.csv".format(struc_root)) # output file 1
        arpeggio_out2 = os.path.join(arpeggio_dir,  "arpeggio_lig_cons_{}.csv".format(struc_root)) # output file 2
        if os.path.isfile(arpeggio_out1) and os.path.isfile(arpeggio_out2):
            lig_cons_split = pd.read_csv(arpeggio_out1)
            arpeggio_lig_cons = pd.read_csv(arpeggio_out2)
        else:
            if len(all_ligs) == 0:
                log.warning("No LOIs in {}".format(struc_root))
                continue
            else:
                lig_cons_split, arpeggio_lig_cons = process_arpeggio(struc, all_ligs, clean_pdbs_dir, arpeggio_dir, mappings_dir) ### NOW PROCESSES ALL LIGANDS ###
                log.debug("Arpeggio output processed for {}!".format(struc_root))
        ligand_contact = arpeggio_lig_cons["PDB_ResNum"].astype(str)
        ligand_contact_list.append(ligand_contact)

    log.info("ARPEGGIO section completed")

    ### BINDING SITE DEFINITION SECTION

    pdb_paths = [os.path.join(clean_pdbs_dir, file) for file in os.listdir(clean_pdbs_dir)]

    ligs = ligs_df.label_comp_id.unique().tolist()
    string_name = "{}_BS_def_{}_{}_{}".format(input_id, lig_clust_method, lig_clust_metric, lig_clust_dist)
    bs_def_out = os.path.join(results_dir, "{}.pkl".format(string_name))
    attr_out = os.path.join(results_dir, "{}.attr".format(string_name))
    chimera_script_out = os.path.join(results_dir, "{}.com".format(string_name))

    ## DEFINE SITES WITHOUT OC

    # GET DISTANCE MATRIX

    rel_dist_out = os.path.join(results_dir, "{}_rel_dist.pkl".format(input_id))
    lig_inters_out = os.path.join(results_dir, "{}_lig_inters.pkl".format(input_id))
    if override or not os.path.isfile(lig_inters_out):
        lig_data_df = generate_ligs_res_df(arpeggio_dir)
        lig_data_df.to_pickle(lig_inters_out)
        log.info("Saved ligand interactions dataframe")
    else:
        lig_data_df = pd.read_pickle(lig_inters_out)
        log.info("Loaded ligand interactions dataframe")

    if override or not os.path.isfile(rel_dist_out):
        lig_res = lig_data_df.binding_res.tolist()
        intersect_dict = get_intersect_rel_matrix(lig_res)
        irel_df = pd.DataFrame(intersect_dict).round(3)
        dist_df = 1 - irel_df # distance matrix in pd.Dataframe() format
        dist_df.to_pickle(rel_dist_out) 
        log.info("Saved relative distance matrix")
    else:
        dist_df = pd.read_pickle(rel_dist_out)
        log.info("Loaded relative distance matrix")

    if os.path.isfile(bs_def_out) and os.path.isfile(attr_out) and os.path.isfile(chimera_script_out):
        lig_data_df = pd.read_pickle(bs_def_out)
        log.debug("Loaded binding site definition")
        pass
    else:
        labs = lig_data_df.lab.tolist()
        condensed_dist_mat = scipy.spatial.distance.squareform(dist_df) # condensed distance matrix to be used for clustering
        linkage = scipy.cluster.hierarchy.linkage(condensed_dist_mat, method = lig_clust_method, optimal_ordering = True)
        cut_tree = scipy.cluster.hierarchy.cut_tree(linkage, height = lig_clust_dist)
        cluster_ids = [int(cut) for cut in cut_tree]
        cluster_id_dict = {labs[i]: cluster_ids[i] for i in range(len(labs))} #dictionary indicating membership for each lig

        pdb_paths = [os.path.join(clean_pdbs_dir, file) for file in os.listdir(clean_pdbs_dir)]
        lig_data_df["binding_site"] = lig_data_df.lab.map(cluster_id_dict)
        #pdb_files_dict = {f.split("/")[-1].split(".")[0]: f.split("/")[-1] for f in pdb_paths}
        pdb_files_dict = {}
        for f in pdb_paths:
            file_name = os.path.basename(f)
            file_root, _ = os.path.splitext(os.path.splitext(file_name)[0])
            if file_root not in pdb_files_dict.keys():
                pdb_files_dict[file_root] = file_name
        
        print(pdb_files_dict)

        lig_data_df["pdb_path"] = lig_data_df.pdb_id.map(pdb_files_dict)

        #print(pdb_files_dict)

        #print(lig_data_df.head())

        write_bs_files(lig_data_df, bs_def_out, attr_out, chimera_script_out)

        log.info("Binding site were defined")

    log.info("Binding site definition completed")

    lig_bs_out = os.path.join(results_dir, "{}_lig_bs.pkl".format(input_id))
    lig_ress_out = os.path.join(results_dir, "{}_lig_ress.pkl".format(input_id))
    site_ress_out = os.path.join(results_dir, "{}_site_ress.pkl".format(input_id))
    res_bss_out = os.path.join(results_dir, "{}_res_bss.pkl".format(input_id))

    if override or not os.path.isfile(lig_bs_out):
        bs_lig_dict = dict(zip(lig_data_df.lab,lig_data_df.binding_site)) # {ligand_id: binding_site_id}
        dump_pickle(bs_lig_dict, lig_bs_out)
        log.info("Saved ligand --> binding site dictionary!")
    else:
        bs_lig_dict = load_pickle(lig_bs_out)
        log.debug("Loaded ligand --> binding site dictionary!")

    if override or not os.path.isfile(lig_ress_out):
        ligand_ress_dict = dict(zip(lig_data_df.lab,lig_data_df.binding_res)) # {ligand_id: [binding_residues]]}
        dump_pickle(ligand_ress_dict, lig_ress_out)
        log.info("Saved ligand --> [binding residues] dictionary!")
    else:
        ligand_ress_dict = load_pickle(lig_ress_out)
        log.debug("Loaded ligand --> [binding residues] dictionary!")

    if override or not os.path.isfile(site_ress_out):
        site_ress_dict = {}
        for site_id, site_rows in lig_data_df.groupby("binding_site"):
            site_ress_dict[site_id] = []
            for _, site_row in site_rows.iterrows():
                site_ress_dict[site_id].extend(site_row.binding_res)
        site_ress_dict = {k: sorted(list(set(v))) for k, v in site_ress_dict.items()} # {binding_site_id: [binding_residues]}
        dump_pickle(site_ress_dict, site_ress_out)
        log.info("Saved binding site --> [binding residues] dictionary!")
    else:
        site_ress_dict = load_pickle(site_ress_out)
        log.debug("Loaded binding site --> [binding residues] dictionary!")

    if override or not os.path.isfile(res_bss_out):
        ress_bss_dict = get_residue_bs_membership(site_ress_dict) #  {residue: [binding_site_ids]}
        dump_pickle(ress_bss_dict, res_bss_out)
        log.info("Saved residue --> [binding site ids] dictionary!")
    else:
        ress_bss_dict = load_pickle(res_bss_out)
        log.debug("Loaded residue --> [binding site ids] dictionary!")

    ### DSSP DATA ANALYSIS
    dsspd_filt = mapped_dssp_df.query('UniProt_ResNum == UniProt_ResNum and PDB_ResName != "X" and RSA == RSA').copy()
    dsspd_filt.SS = dsspd_filt.SS.fillna("C")
    dsspd_filt.SS = dsspd_filt.SS.replace("", "C")
    dsspd_filt.UniProt_ResNum = dsspd_filt.UniProt_ResNum.astype(int)

    AA_dict_out = os.path.join(results_dir, "{}_ress_AA.pkl".format(input_id))
    RSA_dict_out = os.path.join(results_dir, "{}_ress_RSA.pkl".format(input_id))
    SS_dict_out = os.path.join(results_dir, "{}_ress_SS.pkl".format(input_id))
    rsa_profs_out = os.path.join(results_dir, "{}_bss_RSA_profiles.pkl".format(input_id))
    ss_profs_out = os.path.join(results_dir, "{}_bss_SS_profiles.pkl".format(input_id))
    aa_profs_out = os.path.join(results_dir, "{}_bss_AA_profiles.pkl".format(input_id))

    if override or not os.path.isfile(AA_dict_out):
        ress_AA_dict = {
            up_resnum: dsspd_filt.query('UniProt_ResNum == @up_resnum').PDB_ResName.mode()[0] # gets dict per UP residue and more frequent AA.
            for up_resnum in dsspd_filt.UniProt_ResNum.unique().tolist()
        }
        dump_pickle(ress_AA_dict, AA_dict_out)
        log.info("Saved AA dict")
    else:
        ress_AA_dict = load_pickle(AA_dict_out)
        log.debug("Loaded AA dict")
    
    if override or not os.path.isfile(RSA_dict_out):
        ress_RSA_dict = {
            up_resnum: round(dsspd_filt.query('UniProt_ResNum == @up_resnum').RSA.mean(), 2) # gets dict per UP residue and average RSA.
            for up_resnum in dsspd_filt.UniProt_ResNum.unique().tolist()
        }
        dump_pickle(ress_RSA_dict, RSA_dict_out)
        log.info("Saved RSA dict")
    else:
        ress_RSA_dict = load_pickle(RSA_dict_out)
        log.debug("Loaded RSA dict")

    if override or not os.path.isfile(SS_dict_out):
        ress_SS_dict = {
            up_resnum: dsspd_filt.query('UniProt_ResNum == @up_resnum').SS.mode()[0] # gets dict per UP residue and more frequent SS.
            for up_resnum in dsspd_filt.UniProt_ResNum.unique().tolist()
        }
        dump_pickle(ress_SS_dict, SS_dict_out)
        log.info("Saved SS dict")
    else:
        ress_SS_dict = load_pickle(SS_dict_out)
        log.debug("Loaded SS dict")

    if override or not os.path.isfile(rsa_profs_out):
        rsa_profiles = {}
        for k, v in site_ress_dict.items():
            rsa_profiles[k] = []
            for v2 in v:
                if v2 in ress_RSA_dict:
                    rsa_profiles[k].append(ress_RSA_dict[v2])
                else:
                    log.warning("Cannot find RSA data for UP residue {} in {}".format(str(v2), input_id))
        dump_pickle(rsa_profiles, rsa_profs_out)
        log.info("Saved RSA profiles")
    else:
        rsa_profiles = load_pickle(rsa_profs_out)
        log.debug("Loaded RSA profiles")
    
    if override or not os.path.isfile(ss_profs_out):
        ss_profiles = {}
        for k, v in site_ress_dict.items():
            ss_profiles[k] = []
            for v2 in v:
                if v2 in ress_SS_dict:
                    ss_profiles[k].append(ress_SS_dict[v2])
                else:
                    log.warning("Cannot find SS data for UP residue {} in {}".format(str(v2), input_id))
        dump_pickle(ss_profiles, ss_profs_out)
        log.info("Saved SS profiles")
    else:
        ss_profiles = load_pickle(ss_profs_out)
        log.debug("Loaded SS profiles")

    if override or not os.path.isfile(aa_profs_out):
        aa_profiles = {}
        for k, v in site_ress_dict.items():
            aa_profiles[k] = []
            for v2 in v:
                if v2 in ress_AA_dict:
                    aa_profiles[k].append(ress_AA_dict[v2])
                else:
                    log.warning("Cannot find AA data for UP residue {} in {}".format(str(v2), input_id))
        dump_pickle(aa_profiles, aa_profs_out)
        log.info("Saved AA profiles")
    else:
        aa_profiles = load_pickle(aa_profs_out)
        log.debug("Loaded AA profiles")

    log.info("DSSP analysis completed")

    ### VARIATION SECTION

    if run_variants:

        example_struc = os.path.join(clean_pdbs_dir, os.listdir(clean_pdbs_dir)[0])
        fasta_path = os.path.join(varalign_dir, "{}.fa".format(input_id))
        fasta_root, _ = os.path.splitext(fasta_path)    
        hits_aln = "{}.sto".format(fasta_root)  
        hits_aln_rf = "{}_rf.sto".format(fasta_root)
        shenkin_out = os.path.join(varalign_dir, "{}_shenkin.csv".format(input_id))
        shenkin_filt_out = os.path.join(varalign_dir, "{}_shenkin_filt.csv".format(input_id))

        if override_variants or not os.path.isfile(hits_aln_rf):
            create_alignment_from_struc(example_struc, fasta_path)
            log.info("Saved MSA to file")
            pass
        else:
            log.debug("Loaded MSA from file")

        ### CONSERVATION ANALYSIS

        prot_cols = prot_cols = get_target_prot_cols(hits_aln)
        
        shenkin_out = os.path.join(varalign_dir, "{}_rf_shenkin.pkl".format(input_id))
        if override_variants or not os.path.isfile(shenkin_out):
            shenkin = calculate_shenkin(hits_aln_rf, "stockholm", shenkin_out)
            log.info("Saved conservation data table")
        else:
            shenkin = pd.read_pickle(shenkin_out)
            log.debug("Loaded conservation data table")
        
        shenkin_filt_out = os.path.join(varalign_dir, "{}_rf_shenkin_filt.pkl".format(input_id))
        if override_variants or not os.path.isfile(shenkin_filt_out):
            shenkin_filt = format_shenkin(shenkin, prot_cols, shenkin_filt_out)
            log.info("Saved filtered conservation data table")
        else:
            shenkin_filt = pd.read_pickle(shenkin_filt_out)
            log.debug("Loaded filtered conservation data table")

        aln_obj = Bio.AlignIO.read(hits_aln_rf, msa_fmt) #crashes if target protein is not human!
        aln_info_path = os.path.join(varalign_dir, "{}_rf_info_table.p.gz".format(input_id))
        if override_variants or not os.path.isfile(aln_info_path):
            aln_info = varalign.alignments.alignment_info_table(aln_obj)
            aln_info.to_pickle(aln_info_path)
            log.info("Saved MSA info table")
        else:
            aln_info = pd.read_pickle(aln_info_path)
            log.debug("Loaded MSA info table")
        
        log.info("There are {} sequences in MSA".format(len(aln_info)))
        
        indexed_mapping_path = os.path.join(varalign_dir, "{}_rf_mappings.p.gz".format(input_id))
        if override_variants or not os.path.isfile(indexed_mapping_path):
            indexed_mapping_table = varalign.align_variants._mapping_table(aln_info) # now contains all species
            indexed_mapping_table.to_pickle(indexed_mapping_path) # important for merging later on
            log.info("Saved MSA mapping table")
        else:
            indexed_mapping_table = pd.read_pickle(indexed_mapping_path)
            log.debug("Loaded MSA mapping table")    

        aln_info_human = aln_info.query('species == "HUMAN"')

        if len(aln_info_human) > 0:
            log.info("There are {} HUMAN sequences in the MSA".format(len(aln_info_human)))

            human_hits_msa = os.path.join(varalign_dir, "{}_rf_human.sto".format(input_id))

            if override_variants or not os.path.isfile(human_hits_msa):
                get_human_subset_msa(hits_aln_rf, human_hits_msa)
            else:
                pass
            ### copy ensemble SQLite to directory where this is being executed
            cp_path = cp_sqlite(wd)
            log.debug("ENSEMBL_CACHE SQLite copied correctly")

            variant_table_path = os.path.join(varalign_dir, "{}_rf_human_variants.p.gz".format(input_id))
            if override_variants or not os.path.isfile(variant_table_path):
                try:
                    variants_table = varalign.align_variants.align_variants(aln_info_human, path_to_vcf = gnomad_vcf,  include_other_info = False, write_vcf_out = False)     
                except ValueError as e:
                    variants_table = pd.DataFrame()
                    log.warning("No variants were retrieved")

                variants_table.to_pickle(variant_table_path)

            else:
                variants_table = pd.read_pickle(variant_table_path)

            ### remove ensembl SQLite from directory where this is being executed
            rm_sqlite(cp_path)
            log.debug("ENSEMBL_CACHE SQLite removed correctly")

            if variants_table.empty: # variant table is empty. E.g., P03915. Only 3 human sequences. They are all mitochondrial (not in gnomAD)
                pass

            else:
                # in order to be able to read the vcf and parse the DB, the ensemble.cache.sqlite file must be in the ./.varalign directory

                human_miss_vars = format_variant_table(variants_table, prot_cols) # GET ONLY MISSENSE VARIANTS ROWS
                human_miss_vars_msa_out = os.path.join(varalign_dir, "{}_rf_human_missense_variants_seqs.sto".format(input_id))

                miss_df_out = os.path.join(results_dir, "{}_missense_df.pkl".format(input_id))
                
                if override or not os.path.isfile(miss_df_out): # we leave it as override and not override_variants to fix the wrong pseudocounts
                    missense_variants_df = get_missense_df(
                        hits_aln_rf, human_miss_vars,
                        shenkin_filt, prot_cols, human_miss_vars_msa_out
                    )

                    if missense_variants_df.empty:
                        log.warning("No missense variants found for MSA")
                        pass

                    else:
                        missense_variants_df = add_miss_class(
                            missense_variants_df, miss_df_out,
                            cons_col = "abs_norm_shenkin",
                        )
                        log.info("Saved missense dataframe")
                else:
                    missense_variants_df = pd.read_pickle(miss_df_out)
                    log.debug("Loaded missense dataframe")

                if missense_variants_df.empty:
                    log.warning("No missense variants found for MSA of {}".format(input_id))
                    pass
                else:
                    # ADDS COLUMNS FROM MISSENSE DF TO SHENKIN FILT DF, CONSERVATION AND VARIATION DATA ABOUT HUMAN VARIANT SUB MSA
                    shenkin_filt.loc[:, "human_shenkin"] = missense_variants_df.shenkin
                    shenkin_filt.loc[:, "human_occ"] = missense_variants_df.occ
                    shenkin_filt.loc[:, "human_gaps"] = missense_variants_df.gaps
                    shenkin_filt.loc[:, "human_occ_pct"] = missense_variants_df.occ_pct
                    shenkin_filt.loc[:, "human_gaps_pct"] = missense_variants_df.gaps_pct
                    shenkin_filt.loc[:, "variants"] = missense_variants_df.variants
                    shenkin_filt.loc[:, "oddsratio"] = missense_variants_df.oddsratio
                    shenkin_filt.loc[:, "pvalue"] = missense_variants_df.pvalue
                    shenkin_filt.loc[:, "se_OR"] = missense_variants_df.se_OR

        else:
            log.warning("No human sequences in MSA")
            pass

        shenkin_mapped_out = os.path.join(results_dir, "{}_ress_consvar.pkl".format(input_id))
        if override or not os.path.isfile(shenkin_mapped_out): # we leave it as override and not override_variants to fix the wrong pseudocounts
            aln_ids = list(set([seqid[0] for seqid in indexed_mapping_table.index.tolist() if uniprot_id in seqid[0]])) # THIS IS EMPTY IF QUERY SEQUENCE IS NOT FOUND
            n_aln_ids = len(aln_ids)
            if n_aln_ids != 1:
                log.warning("There are {} sequences matching input protein accession".format(str(n_aln_ids)))
            mapped_data = merge_shenkin_df_and_mapping(shenkin_filt, indexed_mapping_table, aln_ids)
            mapped_data.to_pickle(shenkin_mapped_out)
        else:
            mapped_data = pd.read_pickle(shenkin_mapped_out)
        log.info("Saved conservation + variant data")

        if not mapped_dssp_df.empty:
            mapped_data["AA"] = mapped_data.UniProt_ResNum.map(ress_AA_dict)
            mapped_data["RSA"] = mapped_data.UniProt_ResNum.map(ress_RSA_dict)
            mapped_data["SS"] = mapped_data.UniProt_ResNum.map(ress_SS_dict)
        else:
            log.warning("No DSSP data available")
            pass

        mapped_data["binding_sites"] = mapped_data.UniProt_ResNum.map(ress_bss_dict)
        mapped_data.to_pickle(final_table_out)
        log.info("Saved final table")
    else:
        log.info("Not running variants")

    log.info("Variants section completed")

    log.info("THE END")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Clusters ligands, define, and characterise binding sites.")
    parser.add_argument("input_dir", type = str, help = "Path to directory containing input structures")
    parser.add_argument("uniprot_id", type = str, help = "Uniprot ID of the protein")
    parser.add_argument("--struc_fmt", type = str, choices=["pdb", "mmcif"], default = "mmcif", help="Format of the input structures (must be 'pdb' or 'mmcif')")
    parser.add_argument("--override", help = "Override any previously generated files.", action = "store_true")
    parser.add_argument("--override_variants", help = "Override any previously generated files (ONLY VARIANTS SECTION).", action = "store_true")
    parser.add_argument("--variants", help = "Retrieves Human variants form MSA and generates tables.", action = "store_true")

    args = parser.parse_args()

    main(args)
