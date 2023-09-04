### IMPORTS ###

import os
import re
import Bio
import math
import scipy
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
swissprot = config["dbs"].get("swissprot")
stampdir = config["other"].get("stampdir")
oc_dist = float(config["clustering"].get("oc_dist"))
oc_metric = config["clustering"].get("oc_metric")
oc_method = config["clustering"].get("oc_method")
mes_sig_t = float(config["other"].get("mes_sig_t"))
msa_fmt = config["other"].get("msa_fmt")

### SETTING UP LOG

logging.basicConfig(filename = "fragsys_custom.log", format = '%(asctime)s %(name)s [%(levelname)-8s] - %(message)s', level = logging.INFO)

log = logging.getLogger("FRAGSYS_CUSTOM")

### FUNCTIONS

def setup_dirs(dirs):
    """
    Creates directories if they don't exist.
    """
    for dirr in dirs:
        if os.path.isdir(dirr):
            continue
        else:
            os.mkdir(dirr)

def generate_STAMP_domains(pdbs_dir, domains_out, roi = "ALL"):
    """
    Genereates domains file, needed to run STAMP.
    """
    with open(domains_out, "w+") as fh:
        for pdb in os.listdir(pdbs_dir):
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
    runs transform to obtain set of transformed coordinates
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
    moves STAMP output files to appropriate directory
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
        
def get_lig_data(supp_pdbs_dir, ligs_df_path):
    """
    From a directory containing a set of structurally superimposed pdbs,
    writes a csv file indicating the name, chain and residue number of the
    ligand(s) of interest in every pdb
    """
    ligs_df = pd.DataFrame([])
    for struc in os.listdir(supp_pdbs_dir):
        struc_path = os.path.join(supp_pdbs_dir, struc)
        df = PDBXreader(inputfile = struc_path).atoms(format_type = "pdb", excluded=())
        hetatm_df = df.query('group_PDB == "HETATM"')
        hetatm_df = df[df.group_PDB == "HETATM"]
        ligs = hetatm_df.label_comp_id.unique().tolist()
        lois = [lig for lig in ligs if lig not in non_relevant]
        for loi in lois:
            loi_df = hetatm_df[hetatm_df.label_comp_id == loi]
            lois_df_un = loi_df.drop_duplicates(["label_comp_id", "label_asym_id"])[["label_comp_id", "label_asym_id", "auth_seq_id"]]
            lois_df_un["struc_name"] = struc
            ligs_df = ligs_df.append(lois_df_un)
    ligs_df = ligs_df[["struc_name","label_comp_id", "label_asym_id", "auth_seq_id"]]
    ligs_df.to_csv(ligs_df_path, index = False)
    return ligs_df

### MAIN FUNCTION

def main(args):
    """
    Main function of the script. Calls all other functions.
    """
    log.info("Logging initiated")

    input_dir = args.input_dir
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

    ### STAMP SECTION

    domains_out = os.path.join(results_dir, "{}_stamp.domains".format(input_id))
    if os.path.isfile(domains_out):
        pass
    else:
        generate_STAMP_domains(input_dir, domains_out)

    prefix = "{}_stamp".format(input_id)
    n_strucs = len(os.listdir(input_dir))
    matrix_file = prefix + "." + str(n_strucs-1)

    if os.path.isfile(os.path.join(stamp_out_dir, matrix_file)):
        pass
    else:
        ec = stamp(
            domains_out,
            prefix, os.path.join(results_dir, prefix + ".out")
        )
        if ec == 0:
            pass
        else:
            log.error("Something went wrong with STAMP")

    c = 0 # counting the number of superseded pdbs
    structure_files = os.listdir(input_dir)
    for file in structure_files:
        if os.path.isfile(os.path.join(supp_pdbs_dir, file)): # only when they alaready have been transformed
            c += 1
    if c == len(structure_files):
        pass
    else:
        log.info("Proceeding to run TRANSFORM")
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

    



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Clusters ligands, define, and characterise binding sites.")
    parser.add_argument("input_dir", type = str, help = "Path to directory containing input structures")
    # parser.add_argument("output_dir", type = str, help = "Path to directory where output will be written")
    parser.add_argument("--override", help = "Override any previously generated files.", action = "store_true")
    parser.add_argument("--override_variants", help = "Override any previously generated files (ONLY VARIANTS SECTION).", action = "store_true")
    parser.add_argument("--variants", help = "Retrieves Human variants form MSA and generates tables.", action = "store_true")

    args = parser.parse_args()

    main(args)

