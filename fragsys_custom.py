### IMPORTS ###

import os
import re
import Bio
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
from prointvar.config import config as cfg
from prointvar.dssp import DSSPrunner, DSSPreader
from scipy.spatial.distance import squareform, pdist
from prointvar.fetchers import download_sifts_from_ebi
from prointvar.fetchers import download_structure_from_pdbe

### DICTIONARIES AND LISTS

arpeggio_suffixes = [
        "atomtypes", "bs_contacts", "contacts", "specific.sift", 
        "sift","specific.siftmatch", "siftmatch", "specific.polarmatch",
        "polarmatch", "ri", "rings", "ari", "amri", "amam", "residue_sifts"
]

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
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
        'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',                   ### ONE OF THESE MIGHT BE ENOUGH ###
        'THR', 'TRP', 'TYR', 'VAL', 'GLX', 'GLI', 'NLE', 'CYC'
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
swissprot_path = config["dbs"].get("swissprot")
stampdir = config["other"].get("stampdir")
oc_dist = float(config["clustering"].get("oc_dist"))
oc_metric = config["clustering"].get("oc_metric")
oc_method = config["clustering"].get("oc_method")
mes_sig_t = float(config["other"].get("mes_sig_t"))
msa_fmt = config["other"].get("msa_fmt")
stamp_ROI = config["other"].get("stamp_roi")
arpeggio_dist = float(config["other"].get("arpeggio_dist"))

### FUNCTIONS

## UTILITIES

def dump_pickle(data, f_out):
    """
    dumps pickle
    """
    with open(f_out, "wb") as f:
        pickle.dump(data, f)
        
def load_pickle(f_in):
    """
    loads pickle
    """
    with open(f_in, "rb") as f:
        data = pickle.load(f)
    return data

## SETUP FUNCTIONS

def setup_dirs(dirs):
    """
    Creates directories if they don't exist.
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
    with open(domains_out, "w+") as fh:
        for pdb in os.listdir(pdbs_dir):
            if pdb[-4:] == ".pdb":
                fh.write("{} {} {{{}}}\n".format(os.path.join(pdbs_dir, pdb), pdb[:-4], roi))

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
    exit_code = os.system(" ".join(args))
    if exit_code != 0:
        print(" ".join(args))
    return exit_code

def transform(matrix):
    """
    Runs TRANSFORM to obtain set of transformed coordinates.
    """
    if "STAMPDIR" not in os.environ:
        os.environ["STAMPDIR"] = stampdir
    args = [transform_bin, "-f", matrix, "-het"]
    exit_code = os.system(" ".join(args))
    if exit_code != 0:
        print(" ".join(args))
    return exit_code

def move_supp_files(unsupp_pdbs_dir, supp_pdbs_dir):
    """
    Moves set of supperimposed coordinate files to appropriate directory.
    """
    struc_files = os.listdir(unsupp_pdbs_dir)
    cwd = os.getcwd()
    for file in struc_files:
        if os.path.isfile(os.path.join(cwd, file)):
            shutil.move(
                os.path.join(cwd, file),
                os.path.join(supp_pdbs_dir, file)
            )

def move_stamp_output(wd, prefix, stamp_out_dir):
    """
    Moves STAMP output files to appropriate directory.
    """
    cwd = os.getcwd()
    stamp_files = sorted([file for file in os.listdir(cwd) if prefix in file]) + ["stamp_rough.trans"]
    for file in stamp_files:
        filepath = os.path.join(cwd, file)
        if os.path.isfile(filepath):
            shutil.move(filepath, os.path.join(stamp_out_dir, file))
    out_from = os.path.join(wd, prefix + ".out")
    out_to = os.path.join(stamp_out_dir, prefix + ".out")
    doms_from = os.path.join(wd, prefix + ".domains")
    doms_to = os.path.join(stamp_out_dir, prefix + ".domains")
    if os.path.isfile(out_from):
        shutil.move(out_from, out_to)
    if os.path.isfile(doms_from):
        shutil.move(doms_from, doms_to)
        
## LIGAND FUNCTIONS

def get_lig_data(supp_pdbs_dir, ligs_df_path):
    """
    From a directory containing a set of structurally superimposed pdbs,
    writes a csv file indicating the name, chain and residue number of the
    ligand(s) of interest in every pdb.
    """
    ligs_df = pd.DataFrame([])
    for struc in os.listdir(supp_pdbs_dir):
        struc_path = os.path.join(supp_pdbs_dir, struc)
        df = PDBXreader(inputfile = struc_path).atoms(format_type = "pdb", excluded=())
        hetatm_df = df.query('group_PDB == "HETATM"')
        ligs = hetatm_df.label_comp_id.unique().tolist()
        lois = [lig for lig in ligs if lig not in non_relevant]
        for loi in lois:
            loi_df = hetatm_df.query('label_comp_id == @loi')
            lois_df_un = loi_df.drop_duplicates(["label_comp_id", "label_asym_id"])[["label_comp_id", "label_asym_id", "auth_seq_id"]]
            lois_df_un["struc_name"] = struc
            ligs_df = ligs_df.append(lois_df_un)
    ligs_df = ligs_df[["struc_name","label_comp_id", "label_asym_id", "auth_seq_id"]]
    ligs_df.to_csv(ligs_df_path, index = False)
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

def retrieve_mapping_from_struc(struc, uniprot_id, struc_dir, sifts_dir, swissprot):
    """
    Retrieves the mapping between the UniProt sequence and the PDB sequence by doing an alignment.
    """
    input_struct = os.path.join(struc_dir, struc)
    pdb_structure = PDBXreader(inputfile = input_struct).atoms(format_type = "pdb") # ProIntVar reads the local file
    
    seq_record = str(swissprot[uniprot_id]["seq"])
    pps = pdb_structure.query('group_PDB == "ATOM"')[['label_comp_id', 'label_asym_id', 'label_seq_id_full']].drop_duplicates().groupby('label_asym_id')  # groupby chain
    pdb_chain_seqs = [(chain, SeqUtils.seq1(''.join(seq['label_comp_id'].values)), seq['label_seq_id_full'].values) for chain, seq in pps] # list of tuples like: [(chain_id, chain_seq, [chain resnums])]
    alignments = [pairwise2.align.globalxs(str(seq_record),chain_seq[1], -5, -1) for chain_seq in pdb_chain_seqs] # list of lists of tuples containing SwissProt seq - PDB chain seq pairwise alignment
    
    maps = []
    for pdb_chain_seq, alignment in zip(pdb_chain_seqs, alignments):
        PDB_UniProt_map = pd.DataFrame(
            [(i, x) for i, x in enumerate(alignment[0][1], start=1)],  # create aligned PDB sequences to dataframe
            columns=['UniProt_ResNum', 'PDB_ResName']
        )
        PDB_UniProt_map = PDB_UniProt_map.assign(UniProt_ResName = list(alignment[0][0]))
        PDB_index = PDB_UniProt_map.query('PDB_ResName != "-"').index
        PDB_UniProt_map = PDB_UniProt_map.assign(PDB_ResNum = pd.Series(pdb_chain_seq[2], index=PDB_index)) # adds PDB_ResNum column
        PDB_UniProt_map = PDB_UniProt_map.assign(PDB_ChainID = pd.Series(pdb_chain_seq[0], index=PDB_index)) # adds PDB_ChainId column
        maps.append(PDB_UniProt_map)
    prointvar_mapping = pd.concat(maps)
    prointvar_mapping = prointvar_mapping[['UniProt_ResNum','UniProt_ResName','PDB_ResName','PDB_ResNum','PDB_ChainID']]
    prointvar_mapping = prointvar_mapping[~prointvar_mapping.PDB_ResNum.isnull()]
    prointvar_mapping.PDB_ResNum = prointvar_mapping.PDB_ResNum.astype(int)
    prointvar_mapping_csv = os.path.join(sifts_dir, struc.replace(".pdb", ".mapping"))
    prointvar_mapping.to_csv(prointvar_mapping_csv, index = False)
    return prointvar_mapping

## DSSP FUNCTIONS

def run_dssp(struc, supp_pdbs_dir, dssp_dir):
    """
    Runs DSSP, saves and return resulting output dataframe
    """
    dssp_csv = os.path.join(dssp_dir, "dssp_" + struc.replace("pdb", "csv")) # output csv filepath
    dssp_out = os.path.join(dssp_dir, struc.replace("pdb", "dssp"))
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
    exit_code = os.system(" ".join(args))
    if exit_code != 0:
        print(" ".join(args))
    return exit_code

def run_arpeggio(pdb_path, lig_name, dist = arpeggio_dist):
    """
    Runs Arpeggio.
    """
    args = [
        arpeggio_python_bin, arpeggio_bin, pdb_path, "-s",
        "RESNAME:{}".format(lig_name), "-i", str(dist),# "-v" removing verbose
    ]
    exit_code = os.system(" ".join(args))
    if exit_code != 0:
        print(" ".join(args))
    return exit_code

def move_arpeggio_output(supp_pdbs_dir, clean_pdbs_dir, arpeggio_dir, pdb_clean_dir, strucs, struc2ligs):
    """
    Moves pdb_clean.py and Arpeggio output files to appropriate directory.
    """  
    for struc in strucs:
        clean_struc = struc.replace(".pdb", ".clean.pdb")
        clean_struc_path_from = os.path.join(supp_pdbs_dir, clean_struc)
        clean_struc_path_to = os.path.join(clean_pdbs_dir, clean_struc)
        if os.path.isfile(clean_struc_path_from):
            shutil.move(clean_struc_path_from, clean_struc_path_to)
        for arpeggio_suff in arpeggio_suffixes:
            for the_lig in struc2ligs[struc]:
                arpeggio_file_from_supp = os.path.join(supp_pdbs_dir, struc[:-3] + "clean_{}".format(the_lig) + "." + arpeggio_suff)
                arpeggio_file_to = os.path.join(arpeggio_dir, struc[:-3] + "clean_{}".format(the_lig) + "." + arpeggio_suff)
                arpeggio_file_from_clean = os.path.join(clean_pdbs_dir, struc[:-3] + "clean_{}".format(the_lig) + "." + arpeggio_suff)
                if os.path.isfile(arpeggio_file_from_supp):
                    shutil.move(arpeggio_file_from_supp, arpeggio_file_to)
                elif os.path.isfile(arpeggio_file_from_clean):
                    shutil.move(arpeggio_file_from_clean, arpeggio_file_to)
        for pdb_clean_suff in pdb_clean_suffixes:
            pdb_clean_file_from = os.path.join(supp_pdbs_dir, struc + "." + pdb_clean_suff)
            pdb_clean_file_to = os.path.join(pdb_clean_dir, struc + "." + pdb_clean_suff)
            if os.path.isfile(pdb_clean_file_from):
                shutil.move(pdb_clean_file_from, pdb_clean_file_to)
            elif os.path.isfile(os.path.join(pdb_clean_dir, struc + "." + pdb_clean_suff)):
                shutil.move(os.path.join(pdb_clean_dir, struc + "." + pdb_clean_suff), pdb_clean_file_to)

def process_arpeggio(struc, all_ligs, clean_pdbs_dir, arpeggio_dir, sifts_dir):
    """
    processes arpeggio output to generate the two tables that will
    be used later in the analysis
    """
    lig_cons_splits = []
    arpeggio_lig_conss = []
        
    for lig in all_ligs:
        arpeggio_out_path_lig = os.path.join(arpeggio_dir, struc.replace("pdb", "clean_{}.bs_contacts".format(lig)))
        
        all_cons_lig = pd.read_table(arpeggio_out_path_lig, header = None)

        lig_cons_split_lig = reformat_arpeggio(all_cons_lig)

        lig_cons_split_lig = add_resnames_to_arpeggio_table(struc, clean_pdbs_dir, lig_cons_split_lig)

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

    struc_mapping = pd.read_csv(os.path.join(sifts_dir, struc.replace(".pdb", ".mapping")))

    mapping_dict = {}
    for chain, chain_df in struc_mapping.groupby("PDB_ChainID"):
        for i, row in chain_df.iterrows():
            mapping_dict[(chain, str(row["PDB_ResNum"]))] = row["UniProt_ResNum"]
    
    #sifts_dict = extend_sifts_mapping_dict(sifts_dict, bio2asym_chain_dict, cif2pdb_chain_dict)
    
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
        log.warning("{} residues lacked mapping from PDB to UniProt at arpeggio_lig_cons for {}".format(new_len1 - old_len1, struc))
    if new_len2 != old_len2:
        log.warning("{} residues lacked mapping from PDB to UniProt at lig_cons_split for {}".format(new_len2 - old_len2, struc))

    ########################### ADDED TO AVOID INCORRECT BS DEFINITION DUE TO LACKING MAPPING TO PDB RESNUMS ###########################

    lig_cons_split.to_csv(os.path.join(arpeggio_dir,  "arpeggio_all_cons_split_" + struc.replace(".pdb", ".csv")), index = False)
    arpeggio_lig_cons.to_csv(os.path.join(arpeggio_dir,  "arpeggio_lig_cons_" + struc.replace(".pdb", ".csv")), index = False)
    
    return lig_cons_split, arpeggio_lig_cons

def reformat_arpeggio(arpeggio_df):
    """
    starts formatting arpeggio table
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

def add_resnames_to_arpeggio_table(structure, clean_pdbs_dir, arpeggio_cons_split):
    """
    adds residue names to arpeggio table, needed for later table mergings
    """
    structure_path = os.path.join(clean_pdbs_dir, structure.replace("pdb", "clean.pdb"))
    pdb_structure = PDBXreader(inputfile = structure_path).atoms(format_type = "pdb", excluded=())
    resnames_dict = {(row.label_asym_id, int(row.label_seq_id)): row.label_comp_id for index, row in pdb_structure.drop_duplicates(['label_asym_id', 'label_seq_id']).iterrows()}
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
    
### FUNCTIONS FOR SITE DEFINITION

def def_bs_oc(results_dir, pdb_files, prot, bs_def_out, attr_out, chimera_script_out, arpeggio_dir, metric = oc_metric, dist = oc_dist, method = oc_method):
    """
    given a set of pdb structures, and other arguments, clusters ligands in space,
    defines binding sites and writes chimera attribute files and chimera script to
    format the superimposed structures to facilitate visualisation
    
    alt_fmt is a boolean I added so it works with slightly different input. Still PDB
    files, but some which coordinates were transformed using PDBe-KB transformation
    matrices. They have different nomenclature and therefore indices to get pdb_files_dict
    must be different
    """
    lig_data_df, labs = generate_ligs_res_df(arpeggio_dir)

    if len(lig_data_df) == 1: #should only happen in the case of only one LOI
        lig_data_df["binding_site"] = 0
    else:
        dis_out = os.path.join(results_dir, "{}_{}.dis".format(prot, metric))
        get_dis_file(lig_data_df, labs, dis_out, metric = metric)
        
        ocout, ec = oc(dis_out, method = method, cut_t = dist)
        oc_dict = oc2dict(ocout)
        cluster_id_dict = {}
        for k, v in oc_dict.items():
            for member in v["members"]:
                cluster_id_dict[member] = v["new_id"]
                
        lig_data_df["lab"] = labs 
        lig_data_df["binding_site"] = lig_data_df.lab.map(cluster_id_dict)
        #print(pdb_files)
        pdb_files_dict = {f.split("/")[-1].split(".")[0]: f.split("/")[-1] for f in pdb_files}
        lig_data_df["pdb_path"] = lig_data_df.pdb_id.map(pdb_files_dict)

    print(lig_data_df)
    
    write_bs_files(lig_data_df, bs_def_out, attr_out, chimera_script_out)
    
    return lig_data_df, labs

def generate_ligs_res_df(arpeggio_dir):
    """
    Given a directory containing processed Arpeggio output files,
    returns a dataset containing information about the ligands of
    interest binding the protein and their labels.
    """
    lig_files = sorted([file for file in os.listdir(arpeggio_dir) if file.startswith("arpeggio_all_cons_split")])
    file_ids, ligs, resnums, chains, ligs_ress = [[], [], [], [], []]
    for file in lig_files:
        file_id = file.split("_")[-1].split(".")[0]
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
    return lig_data_df, labs

def get_dis_file(lig_data_df, labs, out, metric = oc_metric):
    """
    Creates .dis file to be fed to OC.
    """
    lig_res = lig_data_df.binding_res.tolist() #this is a list of lists, each list contains residue numbers interacting with ligand
    if metric == "i_rel":
        intersect_dict = get_intersect_rel_matrix(lig_res)
    elif metric == "i_abs":
        intersect_dict = get_intersect_matrix(lig_res)
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

def get_intersect_matrix(binding_ress):
    """
    Given a set of ligand binding residues, calcualtes a
    similarity matrix between all the different sets of ligand
    binding residues.
    """
    inters = {i: {} for i in range(len(binding_ress))}
    for i in range(len(binding_ress)):
        inters[i][i] = intersection(binding_ress[i], binding_ress[i])
        for j in range(i+1, len(binding_ress)):
            inters[i][j] = intersection(binding_ress[i], binding_ress[j])
            inters[j][i] = inters[i][j]
    return inters

def intersection(l1, l2):
    """
    Calculates intersection.
    """
    I = len(list(set(l1).intersection(l2)))
    return I

def oc(oc_in, type_mat = "sim", method = oc_method, cut_t = oc_dist):
    """
    Runs OC and returns exit code, should be 0 if all is OK.
    """
    oc_out = oc_in.replace(".dis", "_{}_{}_cut_{}.ocout".format(type_mat, method, cut_t))
    args = [
        "/homes/2394007/oc", type_mat, method, "id", "cut", str(cut_t),
        "ps", oc_out[:-6], "<", oc_in, ">", oc_out
    ]
    exit_code = os.system(" ".join(args))
    return oc_out, exit_code

def oc2dict(ocout):
    """
    Parses OC output to generate dict.
    """
    oc_dict = {}
    re_cluster = re.compile("""##\s+(\d+)\s(\d+\.*\d*)\s(\d+)\s*""")
    with open (ocout, "r") as fh:
        s = 0
        n_singleton = 0
        n_cluster = 0
        for line in fh:
            if s == 0:
                if line.startswith("##"):
                    m = re_cluster.match(line)
                    cluster_id, score, cluster_size = m.groups()
                    oc_dict[cluster_id] = {"score": float(score), "size": int(cluster_size), "new_id": n_cluster}
                    n_cluster += 1
                elif line.startswith(" "):
                    oc_dict[cluster_id]["members"] = line.strip().split(" ")
                else:
                    if line.strip() == "":
                        pass
                    elif line.strip() == "UNCLUSTERED ENTITIES":
                        s = 1
            else:
                oc_dict["singleton_{}".format(n_singleton)] = {"score": "", "size": 1, "members": [line.strip(),], "new_id": n_cluster}
                n_singleton += 1
                n_cluster += 1
    return oc_dict

def write_bs_files(frag_mean_coords, bs_def_out, attr_out, chimera_script_out):
    """
    Writes files for binding site definition.
    """
    frag_mean_coords = frag_mean_coords.dropna()
    frag_mean_coords.binding_site = frag_mean_coords.binding_site.astype(int)
    frag_mean_coords.lig_resnum = frag_mean_coords.lig_resnum.astype(int)
    chimera_atom_spec = (':'+ frag_mean_coords.lig_resnum.astype(str) +
                     '.'+ frag_mean_coords.lig_chain +
                     '&#/name==' + frag_mean_coords.pdb_path)
    frag_mean_coords = frag_mean_coords.assign(chimera_atom_spec = chimera_atom_spec)  
    frag_mean_coords.to_csv(bs_def_out, index = False) # saves table to csv
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
            out.write('colour {} :/binding_site=={}\n'.format(','.join(list(map(str, list(sample_colors[i])))), i))
        out.write("### SOME FORMATTING ###\n")
        out.write("\n".join(cmds))
    print("Chimera script successfully created!")

### SETTING UP LOG

logging.basicConfig(filename = "fragsys_custom.log", format = '%(asctime)s %(name)s [%(levelname)-8s] - %(message)s', level = logging.INFO)

log = logging.getLogger("FRAGSYS_CUSTOM")

### MAIN FUNCTION

def main(args):
    """
    Main function of the script. Calls all other functions.
    """
    log.info("Logging initiated")

    input_dir = args.input_dir
    uniprot_id = args.uniprot_id
    override = args.override
    override_variants = args.override_variants
    run_variants = args.variants

    input_id = input_dir.split("/")[-1]
    output_dir = "./output/{}".format(input_id)

    supp_pdbs_dir = os.path.join(output_dir, "supp_pdbs")
    clean_pdbs_dir = os.path.join(output_dir, "clean_pdbs")
    stamp_out_dir = os.path.join(output_dir, "stamp_out")
    results_dir = os.path.join(output_dir, "results")
    sifts_dir = os.path.join(output_dir, "sifts")
    dssp_dir = os.path.join(output_dir, "dssp")
    pdb_clean_dir = os.path.join(output_dir, "pdb_clean")
    arpeggio_dir = os.path.join(output_dir, "arpeggio")
    varalign_dir = os.path.join(output_dir, "varalign")
    
    dirs = [
        #unsupp_cifs_dir, unsupp_pdbs_dir,
        output_dir, supp_pdbs_dir,
        clean_pdbs_dir, stamp_out_dir, results_dir, pdb_clean_dir,
        arpeggio_dir, sifts_dir, dssp_dir, varalign_dir,
        #figs_dir
        ]
    
    setup_dirs(dirs)

    log.info("Directories created")
    strucs = [f for f in os.listdir(input_dir) if f.endswith(".pdb") and "clean" not in f]
    n_strucs = len(strucs)
    log.info("Number of structures: {}".format(n_strucs))

    ### STAMP SECTION

    domains_out = os.path.join(results_dir, "{}_stamp.domains".format(input_id))
    if os.path.isfile(domains_out):
        #log.info("STAMP domains file already exists")
        pass
    else:
        generate_STAMP_domains(input_dir, domains_out)
        log.info("STAMP domains file generated")

    prefix = "{}_stamp".format(input_id)
    matrix_file = prefix + "." + str(n_strucs-1)

    last_matrix_path = os.path.join(stamp_out_dir, matrix_file)

    if os.path.isfile(last_matrix_path):
        #log.info("STAMP matrix files already exist")
        pass
    else:
        ec = stamp(
            domains_out,
            prefix, os.path.join(results_dir, prefix + ".out")
        )
        if ec == 0:
            log.info("STAMP matrix files generated")
            pass
        else:
            log.error("Something went wrong with STAMP")
        

    c = 0 # counting the number of superseded pdbs
    for file in strucs:
        if os.path.isfile(os.path.join(supp_pdbs_dir, file)): # only when they alaready have been transformed
            c += 1
    if c == n_strucs:
        log.info("All structures already are superposed")
        pass
    else:
        if not os.path.isfile(matrix_file): # RUNNING TRANSFORM ONCE STAMP OUTPUT HAS BEEN MOVED TO STAMP_OUT_DIR
            matrix_file = os.path.join(stamp_out_dir, prefix + "." + str(n_strucs-1))
        ec = transform(matrix_file) #running transform with matrix on cwd
        if ec == 0:
            pass
        else:
            log.error("Something went wrong with TRANSFORM")
    
    log.info("STAMP and TRANSFORM completed")
            
    move_supp_files(input_dir, supp_pdbs_dir)

    wd = os.getcwd()

    move_stamp_output(wd, prefix, stamp_out_dir)

    ### GET LIGAND DATA

    lig_data_path = os.path.join(results_dir, "{}_lig_data.csv".format(input_id))
    if os.path.isfile(lig_data_path):
        ligs_df = pd.read_csv(lig_data_path)
    else:
        ligs_df = get_lig_data(supp_pdbs_dir, lig_data_path)

    log.info("Ligand data retrieved")

    ### UNIPROT MAPPING SECTION

    swissprot = load_pickle(swissprot_path)
    log.info("Swissprot loaded")

    pdb_mappings = []
    for struc in strucs:
         ## DSSP
        dssp_csv = os.path.join(dssp_dir, "dssp_" + struc.replace("pdb", "csv"))
        if os.path.isfile(dssp_csv):
            dssp_data = pd.read_csv(dssp_csv)
            #log.info("DSSP data already exists")
            pass
        else:
            dssp_data = run_dssp(struc, supp_pdbs_dir, dssp_dir)
            log.info("DSSP run successfully on {}".format(struc))

        ## UNIPROT MAPPING
        
        struc_mapping_path = os.path.join(sifts_dir, struc.replace(".pdb", ".mapping"))
        if os.path.isfile(struc_mapping_path):
            mapping = pd.read_csv(struc_mapping_path)
            #log.info("Mapping file for {} already exists".format(struc))
            pass
        else:
            mapping = retrieve_mapping_from_struc(struc, uniprot_id, supp_pdbs_dir, sifts_dir, swissprot)
            log.info("Mapping file for {} generated".format(struc))
        
        mapping = pd.merge(mapping, dssp_data, left_on = "PDB_ResNum", right_on = "PDB_ResNum")
        pdb_mappings.append(mapping)

    ### ARPEGGIO PART ###

    struc2ligs = {}
    for struc in strucs:
        struc2ligs[struc] = []
        struc_df = ligs_df.query('struc_name == @struc')
        pdb_path = os.path.join(supp_pdbs_dir, struc)
        clean_pdb_path = pdb_path.replace(".pdb", ".clean.pdb")
        pdb_clean_1 = os.path.join(clean_pdbs_dir, struc.replace(".pdb", ".clean.pdb"))
        if os.path.isfile(pdb_clean_1):
            #("PDB {} already cleaned".format(struc))
            pass
        else:
            ec = run_clean_pdb(pdb_path)
            if ec == 0:
                log.info("PDB {} cleaned".format(struc))
                pass
            else:
                log.error("Something went wrong when cleaning {} :(".format(struc[:4]))
                pass

        ligs = struc_df.label_comp_id.unique().tolist()

        for the_lig in ligs: # RUNs ARPEGGIO ONCE FOR EACH LIGAND
            struc2ligs[struc].append(the_lig)

            if not os.path.isfile(clean_pdb_path):
                clean_pdb_path = pdb_clean_1

            lig_contacts_1 = os.path.join(arpeggio_dir, struc[:-3] + "clean_{}.bs_contacts".format(the_lig))
            lig_contacts_2 = os.path.join(supp_pdbs_dir, struc[:-3] + "clean_{}.bs_contacts".format(the_lig))
            
            if os.path.isfile(lig_contacts_1) or os.path.isfile(lig_contacts_2):
                #("Arpeggio already ran for {} in {}!".format(the_lig, struc))
                continue

            ec = run_arpeggio(clean_pdb_path, the_lig)
            if ec == 0:
                log.info("Arpeggio ran sucessfully for {} in {}!".format(the_lig, struc))
                for arpeggio_suff in arpeggio_suffixes: # CHANGES ARPEGGIO OUTPUT FILENAMES SO THEY INCLUDE LIGAND NAME
                    arpeggio_file_old_name_supp = os.path.join(supp_pdbs_dir, struc[:-3] + "clean" + "." + arpeggio_suff)
                    arpeggio_file_new_name_supp = os.path.join(supp_pdbs_dir, struc[:-3] + "clean_{}".format(the_lig) + "." + arpeggio_suff)
                    arpeggio_file_old_name_clean = os.path.join(clean_pdbs_dir, struc[:-3] + "clean" + "." + arpeggio_suff)
                    arpeggio_file_new_name_clean = os.path.join(clean_pdbs_dir, struc[:-3] + "clean_{}".format(the_lig) + "." + arpeggio_suff)
                    if os.path.isfile(arpeggio_file_old_name_supp):
                        os.rename(arpeggio_file_old_name_supp, arpeggio_file_new_name_supp)
                    elif os.path.isfile(arpeggio_file_old_name_clean):
                        os.rename(arpeggio_file_old_name_clean, arpeggio_file_new_name_clean)
            else:
                log.error("Something went wrong when running Arpeggio for {} :(".format(struc))
                pass

    move_arpeggio_output(supp_pdbs_dir, clean_pdbs_dir, arpeggio_dir, pdb_clean_dir, strucs, struc2ligs)

    ligand_contact_list = []
    for struc in strucs:
        all_ligs = ligs_df.query('struc_name == @struc').label_comp_id.unique().tolist()
        arpeggio_out1 = os.path.join(arpeggio_dir, "arpeggio_all_cons_split_" + struc.replace(".pdb", ".csv")) # output file 1
        arpeggio_out2 = os.path.join(arpeggio_dir,  "arpeggio_lig_cons_" + struc.replace(".pdb", ".csv")) # output file 2
        if os.path.isfile(arpeggio_out1) and os.path.isfile(arpeggio_out2):
            lig_cons_split = pd.read_csv(arpeggio_out1)
            arpeggio_lig_cons = pd.read_csv(arpeggio_out2)
        else:
            if len(all_ligs) == 0:
                log.warning("No LOIs in {}, so skipping!".format(struc))
                continue
            else:
                lig_cons_split, arpeggio_lig_cons = process_arpeggio(struc, all_ligs, clean_pdbs_dir, arpeggio_dir, sifts_dir) ### NOW PROCESSES ALL LIGANDS ###
                #log.info("Arpeggio output processed for {}!".format(struc))
        ligand_contact = arpeggio_lig_cons["PDB_ResNum"].astype(str)
        ligand_contact_list.append(ligand_contact)

    log.info("ARPEGGIO section completed")

    ### BINDING SITE DEFINITION SECTION

    pdb_paths = [os.path.join(clean_pdbs_dir, file) for file in os.listdir(clean_pdbs_dir)]

    ligs = ligs_df.label_comp_id.unique().tolist()
    string_name = "{}_BS_def_OC_{}_{}_{}".format(uniprot_id, oc_method, oc_metric, oc_dist)
    bs_def_out = os.path.join(results_dir, "{}.csv".format(string_name))
    attr_out = os.path.join(results_dir, "{}.attr".format(string_name))
    chimera_script_out = os.path.join(results_dir, "{}.com".format(string_name))
    if os.path.isfile(bs_def_out) and os.path.isfile(attr_out) and os.path.isfile(chimera_script_out):
        log.info("Binding sites already defined!")
        pass
    else:
        def_bs_oc(
            results_dir, pdb_paths, uniprot_id,
            bs_def_out, attr_out, chimera_script_out,
            arpeggio_dir, 
            )
        
        #def_bs_oc(results_dir, pdb_files, prot, bs_def_out, attr_out, chimera_script_out, arpeggio_dir, metric = oc_metric, dist = oc_dist, method = oc_method, alt_fmt = False):
        print("Binding sites were sucessfully defined!")
    
    bs_definition = pd.read_csv(bs_def_out)



    



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Clusters ligands, define, and characterise binding sites.")
    parser.add_argument("input_dir", type = str, help = "Path to directory containing input structures")
    parser.add_argument("uniprot_id", type = str, help = "Uniprot ID of the protein")
    parser.add_argument("--override", help = "Override any previously generated files.", action = "store_true")
    parser.add_argument("--override_variants", help = "Override any previously generated files (ONLY VARIANTS SECTION).", action = "store_true")
    parser.add_argument("--variants", help = "Retrieves Human variants form MSA and generates tables.", action = "store_true")

    args = parser.parse_args()

    main(args)

