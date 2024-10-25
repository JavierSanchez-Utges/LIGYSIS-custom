# LIGYSIS pipeline for user job submission

This repository contains a customised version of our original ligand site analysis [**LIGYSIS**](https://github.com/JavierSanchez-Utges/ligysis_dev) pipeline, **LIGYSIS<sub>CUSTOM</sub>**, used for the analysis of protein-ligand complexes submitted to our **LIGYSIS** web server found [here](https://www.compbio.dundee.ac.uk/ligysis/). The code for the web server can be found [here](https://github.com/JavierSanchez-Utges/ligysis_flask).

## Pipeline methodology

The pipeline can be summarised in the following steps:
1. Entry, protein and gene names extraction from UniProt using the user-provided UniProt accession.
   
   **Note:** Empty strings are returned if the user provided <i>unknown</i> as the `uniprot_id` argument (only to be used if the structures are from protein not present in UniProt).
2. Structure <i>cleaning</i> using [`clean_pdb.py`](https://github.com/harryjubb/pdbtools) script.
3. Structural superimposition using [STAMP](http://www.compbio.dundee.ac.uk/downloads/stamp/).
4. <i>Simplification</i> of superposed files. This consists in keeping protein atoms only for one of the superposed structures, and heteroatoms for the rest. This is done to generate a lower-weight superposition (all ligands to a single protein scaffold).
5. Mapping of PDB residues to UniProt by means of a pairwise alignment: protein chain sequence to UniProt sequence associated to user-provided UniProt accession.

   **Note:** If `uniprot_id` is <i>unknown</i>, this will step will generate a pseudo-mapping where PDB residue numbers are mapped to themselves. For the programme to work correctly, if the protein is not in UniProt, structures should present the same numbering scheme.
7. [Relative solvent accessibility](https://en.wikipedia.org/wiki/Relative_accessible_surface_area) (RSA) and secondary structure elements calculation with [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) via [ProIntVar](https://github.com/bartongroup/prointvar).
8. Protein-ligand interactions calculation with [pdbe-rpeggio](https://github.com/PDBeurope/arpeggio).
9. Ligand clustering into binding sites using [SciPy](https://scipy.org/).
10. Generation of [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) visualisation scripts.
11. Multiple sequence alignment with [jackhmmer](http://hmmer.org/).
12. Shenkin amino acid divergence score calculation [[1](https://doi.org/10.1002/prot.340110408), [2](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009335)].
13. Missense enrichment score calculation with [VarAlign](https://github.com/bartongroup/SM_varalign) [[3](https://www.biorxiv.org/content/10.1101/127050v2), [4](https://onlinelibrary.wiley.com/doi/full/10.1002/pro.3783), [5](https://www.nature.com/articles/s42003-024-06117-5)].
14. RSA-based clustering label and functional score calculation [[6](https://www.nature.com/articles/s42003-024-05970-8)].

The final output of the pipeline consists of multiple tables collating the results from the different steps of the analysis for each residue, and for the defined ligand binding sites. These data include relative solvent accessibility (RSA), secondary structure, PDB/UniProt residue number, alignment column, divergence score, missense enrichment score, p-value, etc.

Refer to notebook [15](analysis/15_ML_predicting_rsa_labels.ipynb) to predict RSA cluster labels for your binding sites of interest.

## Dependencies
Third party dependencies for these notebooks include:
- [pdbe-arpeggio](https://github.com/PDBeurope/arpeggio) [(GNU GPL v3.0 License)](https://github.com/harryjubb/arpeggio/blob/master/LICENSE)
- [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) [(Boost Software License)](https://swift.cmbi.umcn.nl/gv/dssp/)
- [Hmmer](http://hmmer.org/) [(BSD-3 Clause License)](http://eddylab.org/software/hmmer/Userguide.pdf)
- [STAMP](http://www.compbio.dundee.ac.uk/downloads/stamp/) [(GNU GPL v3.0 License)](http://www.compbio.dundee.ac.uk/manuals/stamp.4.4/stamp.html)
- [ProIntVar](https://github.com/bartongroup/prointvar) [(MIT License)](https://github.com/bartongroup/ProIntVar/blob/master/LICENSE.md)
- [ProteoFAV](https://github.com/bartongroup/ProteoFAV) [(MIT License)](https://github.com/bartongroup/ProteoFAV/blob/master/LICENSE.md)
- [VarAlign](https://github.com/bartongroup/SM_varalign) [(MIT License)](https://github.com/bartongroup/SM_VarAlign/blob/master/LICENSE)

Other standard python libraries:
- [Biopython](https://biopython.org/) [(BSD 3-Clause License)](https://github.com/biopython/biopython/blob/master/LICENSE.rst)
- [Keras](https://keras.io/) [(Apache v2.0 License)](https://github.com/keras-team/keras/blob/master/LICENSE)
- [Numpy](https://numpy.org/) [(BSD 3-Clause License)](https://github.com/numpy/numpy/blob/main/LICENSE.txt)
- [Pandas](https://pandas.pydata.org/) [(BSD 3-Clause License)](https://github.com/pandas-dev/pandas/blob/main/LICENSE)
- [Scipy](https://scipy.org/) [(BSD 3-Clause License)](https://github.com/scipy/scipy/blob/main/LICENSE.txt)
- [Scikit-learn](https://scikit-learn.org/stable/) [(BSD 3-Clause License)](https://github.com/scikit-learn/scikit-learn/blob/main/COPYING)
- [Tensorflow](https://www.tensorflow.org/) [(Apache v2.0 License)](https://github.com/tensorflow/tensorflow/blob/master/LICENSE)

For more information on the dependencies, refere to the .yml files in the [`envs`](envs/) directory. To install all the dependencies, refer to the [installation manual](INSTALL.md).

## Running LIGYSIS<sub>CUSTOM</sub>

This is a command line programme, which can be executed like this:

```sh
python fragsys_custom.py IN/Q9UGL1_cif Q9UGL1 pdb
```

The programme has three mandatory arguments:
- `input_dir` is the input directory containing the set of structures to be analysed in either [PDB](https://www.wwpdb.org/documentation/file-format) (<i>.pdb</i>, <i>.ent</i>) or [mmCIF](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/beginner%E2%80%99s-guide-to-pdbx-mmcif) (<i>.cif</i>) formats. In this example:`IN/Q9UGL1_cif`.

- `uniprot_id` is the corresponding [UniProt](https://www.uniprot.org/) accession identifier of the structures within `input_dir`. In this example: [Q9UGL1](https://www.uniprot.org/uniprotkb/Q9UGL1/entry), which corresponds to human Lysine-specific demethylase 5B.

- Finally, `struc_fmt`, indicating whether the structures are in `pdb` or `mmcif` format. Only these two formats are selected, and the programme will not run properly unless a structure format is required.

## Help and manual

To get help or information about the programme, run:

```sh
python fragsys_custom.py -h
```

which will print the manual of the programme:

```
usage: fragsys_custom.py [-h] [--override] [--override_variants] [--variants]
                         [--clust_method CLUST_METHOD]
                         [--clust_dist CLUST_DIST] [--hmm_iters HMM_ITERS]
                         [--cons_thresh_high CONS_THRESH_HIGH]
                         [--cons_thresh_low CONS_THRESH_LOW]
                         [--mes_thresh MES_THRESH]
                         input_dir uniprot_id {pdb,mmcif}

Clusters ligands, defines, and characterises binding sites.

positional arguments:
  input_dir             Path to directory containing input structures
  uniprot_id            UniProt ID of the protein
  {pdb,mmcif}           Format of the input structures (must be 'pdb' or
                        'mmcif')

optional arguments:
  -h, --help            show this help message and exit
  --override            Override any previously generated files.
  --override_variants   Override any previously generated files (ONLY VARIANTS
                        SECTION).
  --variants            Retrieves Human variants from MSA and generates
                        tables.
  --clust_method CLUST_METHOD
                        Ligand clustering method (default: average)
  --clust_dist CLUST_DIST
                        Ligand clustering distance threshold (default: 0.50)
  --hmm_iters HMM_ITERS
                        Number of iterations for JACKHMMER (default: 3)
  --cons_thresh_high CONS_THRESH_HIGH
                        Conservation high threshold (default: 75)
  --cons_thresh_low CONS_THRESH_LOW
                        Conservation low threshold (default: 25)
  --mes_thresh MES_THRESH
                        MES threshold (default: 1.0)
```

## References
1. Shenkin PS, Erman B, Mastrandrea LD. Information-theoretical entropy as a measure of sequence variability.
Proteins. 1991; 11(4):297–313. Epub 1991/01/01. [https://doi.org/10.1002/prot.340110408](https://doi.org/10.1002/prot.340110408)
PMID: 1758884.

2. **Utgés JS**, Tsenkov MI, Dietrich NJM, MacGowan SA, Barton GJ. Ankyrin repeats in context with human population variation. PLoS Comput Biol. 2021 Aug 24;17(8):e1009335. doi: [10.1371/journal.pcbi.1009335](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009335). PMID: 34428215; PMCID: PMC8415598.
 
3. MacGowan, SA, Madeira, F, Britto-Borges, T, Schmittner, MS, Cole, C, & Barton, GJ (2017). Human missense variation is constrained by domain structure and highlights functional and pathogenic residues. bioRxiv, 127050. [https://doi.org/10.1101/127050](https://www.biorxiv.org/content/10.1101/127050v2).

4. MacGowan SA, Madeira F, Britto-Borges T, Warowny M, Drozdetskiy A, Procter JB, Barton GJ. The Dundee Resource for Sequence Analysis and Structure Prediction. Protein Sci. 2020 Jan;29(1):277-297. doi: [10.1002/pro.3783](https://onlinelibrary.wiley.com/doi/full/10.1002/pro.3783). Epub 2019 Nov 28. PMID: 31710725; PMCID: PMC6933851.

5. MacGowan SA, Madeira F, Britto-Borges T, Barton GJ. A unified analysis of evolutionary and population constraint in protein domains highlights structural features and pathogenic sites. Commun Biol. 2024 Apr 11;7(1):447. doi: [10.1038/s42003-024-06117-5](https://www.nature.com/articles/s42003-024-06117-5). PMID: 38605212; PMCID: PMC11009406.
   
6. **Utgés JS**, MacGowan SA, Ives CM, Barton GJ. Classification of likely functional class for ligand binding sites identified from fragment screening. Commun Biol. 2024 Mar 13;7(1):320. doi: [10.1038/s42003-024-05970-8](https://www.nature.com/articles/s42003-024-05970-8). PMID: 38480979; PMCID: PMC10937669.
