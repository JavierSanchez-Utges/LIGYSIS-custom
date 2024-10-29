# LIGYSIS pipeline for user job submission

This repository contains a customised version of our original ligand site analysis [**LIGYSIS**](https://github.com/JavierSanchez-Utges/LIGYSIS) pipeline, which was employed to analyse all biologically relevant protien-ligand interactions on the [PDBe](https://www.ebi.ac.uk/pdbe/) [[1](https://europepmc.org/article/MED/31691821)], which results are served in the **LIGYSIS** [web server](https://www.compbio.dundee.ac.uk/ligysis/). The code for the web server can be found [here](https://github.com/JavierSanchez-Utges/LIGYSIS-web).

This customised version of **LIGYSIS** does not rely on the [PDBe-KB](https://www.ebi.ac.uk/pdbe/pdbe-kb/) [[2](https://academic.oup.com/nar/article/48/D1/D344/5580911)] and PDBe [APIs](https://www.ebi.ac.uk/pdbe/pdbe-rest-api) [[3](https://academic.oup.com/bioinformatics/article/37/21/3950/6291664)], but instead processes from scratch any set of structures in <i>.pdb</i> or <i>.cif</i> formats, and is employed to process user jobs in the [LIGYSIS web server](https://www.compbio.dundee.ac.uk/ligysis/).

## Pipeline methodology

The pipeline can be summarised in the following steps:
1. Entry, protein and gene names extraction from UniProt using the user-provided [UniProt](https://www.uniprot.org/) [[4](https://academic.oup.com/nar/article/49/D1/D480/6006196?login=true)] accession.
   
   **Note:** Empty strings are returned if the user provided <i>unknown</i> as the `uniprot_id` argument (only to be used if the structures are from protein not present in UniProt).
2. Structure <i>cleaning</i> using [`clean_pdb.py`](https://github.com/harryjubb/pdbtools) script.
3. Structural superimposition using [STAMP](http://www.compbio.dundee.ac.uk/downloads/stamp/) [[5](https://onlinelibrary.wiley.com/doi/10.1002/prot.340140216)].
4. <i>Simplification</i> of superposed files. This consists in keeping protein atoms only for one of the superposed structures, and heteroatoms for the rest. This is done to generate a lower-weight superposition (all ligands to a single protein scaffold).
5. Mapping of PDB residues to UniProt by means of a pairwise alignment: protein chain sequence to UniProt sequence associated to user-provided UniProt accession.

   **Note:** If `uniprot_id` is <i>unknown</i>, this will step will generate a pseudo-mapping where PDB residue numbers are mapped to themselves. For the programme to work correctly, if the protein is not in UniProt, structures should present the same numbering scheme.
7. [Relative solvent accessibility](https://en.wikipedia.org/wiki/Relative_accessible_surface_area) (RSA) [[6](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0080635#:~:text=It%20is%20defined%20as%20a,1%5D%E2%80%93%5B5%5D.)] and secondary structure elements calculation with [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) [[7](https://onlinelibrary.wiley.com/doi/10.1002/bip.360221211)] via [ProIntVar](https://github.com/bartongroup/prointvar) [[8](https://onlinelibrary.wiley.com/doi/full/10.1002/pro.3783)].
8. Protein-ligand interactions calculation with [pdbe-rpeggio](https://github.com/PDBeurope/arpeggio) [[9](https://www.sciencedirect.com/science/article/pii/S0022283616305332?via%3Dihub)].
9. Ligand clustering into binding sites using [SciPy](https://scipy.org/) [[10](https://www.nature.com/articles/s41592-019-0686-2)].
10. Generation of [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) [[11](https://onlinelibrary.wiley.com/doi/10.1002/pro.3943)] visualisation scripts.
11. Multiple sequence alignment with [jackhmmer](http://hmmer.org/) [[12](https://academic.oup.com/bioinformatics/article/14/9/755/259550)].
12. Shenkin amino acid divergence score calculation [[13](https://doi.org/10.1002/prot.340110408), [14](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009335)].
13. Missense enrichment score calculation with [VarAlign](https://github.com/bartongroup/SM_varalign) [[15](https://www.biorxiv.org/content/10.1101/127050v2), [16](https://www.nature.com/articles/s42003-024-06117-5)].
14. RSA-based clustering label and functional score calculation [[17](https://www.nature.com/articles/s42003-024-05970-8)].

**Note:** this is done on a separate script: [predict_rsa_labels.py](predict_rsa_labels.py), which requires a different environment due to dependencies crashing. This programme must be executed on the `DEEP_LEAERNING` [environment](ENVS/DEEP_LEARNING.yml).

The final output of the pipeline consists of multiple tables collating the results from the different steps of the analysis for each residue, and for the defined ligand binding sites. These data include relative solvent accessibility (RSA), secondary structure, PDB/UniProt residue number, alignment column, divergence score, missense enrichment score, p-value, etc.

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
- [BioPython](https://biopython.org/) [(BSD 3-Clause License)](https://github.com/biopython/biopython/blob/master/LICENSE.rst) [[18](https://academic.oup.com/bioinformatics/article/25/11/1422/330687)]
- [Keras](https://keras.io/) [(Apache v2.0 License)](https://github.com/keras-team/keras/blob/master/LICENSE) [19]
- [NumPy](https://numpy.org/) [(BSD 3-Clause License)](https://github.com/numpy/numpy/blob/main/LICENSE.txt) [[20](https://www.nature.com/articles/s41586-020-2649-2)]
- [Pandas](https://pandas.pydata.org/) [(BSD 3-Clause License)](https://github.com/pandas-dev/pandas/blob/main/LICENSE)  [[21](https://doi.org/10.25080/Majora-92bf1922-00a)]
- [SciPy](https://scipy.org/) [(BSD 3-Clause License)](https://github.com/scipy/scipy/blob/main/LICENSE.txt)
- [Scikit-learn](https://scikit-learn.org/stable/) [(BSD 3-Clause License)](https://github.com/scikit-learn/scikit-learn/blob/main/COPYING) [[22](https://www.jmlr.org/papers/v12/pedregosa11a.html)]
- [TensorFlow](https://www.tensorflow.org/) [(Apache v2.0 License)](https://github.com/tensorflow/tensorflow/blob/master/LICENSE) [[23](https://static.googleusercontent.com/media/research.google.com/en//pubs/archive/45166.pdf)]

For more information on the dependencies, refer to the .yml files in the [`ENVS`](ENVS/) directory. To install all the dependencies, refer to the [installation manual](INSTALL.md).

## Environments

The `ENVS` folder contains three `.yml` files describing the necessary packages and dependencies for the different parts of the pipeline and analysis.
  -  [LIGYSIS](ENVS/LIGYSIS.yml) is needed to run **LIGYSIS**.
    
  -  [DEEP_LEARNING](ENVS/DEEP_LEARNING.yml) contains the packages necessary to do predict the RSA Cluster labels and functional scores with [predict_rsa_labels.py](predict_rsa_labels.py).

Note that there are no `.yml` files for the `ARPEGGIO`, `CLEAN_PDB`, `DSSP`, `HMMER` environments, as these are created from the command line without the need of a `.yml` file.

## Installation

For complete installation instructions refer [here](INSTALL.md).

### Downloading SwissProt

This is the database used for our analysis, but can be changed according to the user purposes, e.g. TrEMBL. What is important is to add the correct path in the [configuration file](ligysis_config.txt). To download SwissProt, follow the next steps.

```sh
# download SwissProt in fasta format (88MB)
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# decrompress the file
gzip -d uniprot_sprot.fasta.gz
```

### Downloading gnomAD v2.1

This is the database used for our analysis, but can be changed according to the user purposes, e.g. v > 2.1. What is important is to add the correct path in the [configuration file](ligysis_config.txt). To download gnomAD v2.1 [[24](https://www.nature.com/articles/s41586-020-2308-7)], follow the next steps.
```
# download gnomAD Exomves vcf (large file 58GB)
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz
```

For more information, refer to [gnomAD](https://gnomad.broadinstitute.org/).

After downloading gnomAD, it is required to run [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) [[25](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0974-4)] on it, as VarAlign uses its annotations. Information on how to do so [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html).

## Execution

**LIGYSIS** can be run like this:

```sh
python ligysis_custom.py IN/Q9UGL1_cif Q9UGL1 pdb
```

This needs to be done within the `LIGYSIS` [environment](ENVS/LIGYSIS.yml).

The programme uses relative paths, so it is recommended to run it in the repository directory, where it can directly read from `./IN` and write output to `./OUT` (you need to create this second directory).

The programme has three mandatory arguments:
- `input_dir` is the input directory containing the set of structures to be analysed in either [PDB](https://www.wwpdb.org/documentation/file-format) (<i>.pdb</i>, <i>.ent</i>) or [mmCIF](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/beginner%E2%80%99s-guide-to-pdbx-mmcif) (<i>.cif</i>) formats. In this example:`IN/Q9UGL1_cif`.

- `uniprot_id` is the corresponding [UniProt](https://www.uniprot.org/) accession identifier of the structures within `input_dir`. In this example: [Q9UGL1](https://www.uniprot.org/uniprotkb/Q9UGL1/entry), which corresponds to human Lysine-specific demethylase 5B.

- Finally, `struc_fmt`, indicating whether the structures are in `pdb` or `mmcif` format. Only these two formats are selected, and the programme will not run properly unless a structure format is required.

To carry out the last step and add the RSA-derived Cluster labels and functional scores, the [predict_rsa_labels.py](predict_rsa_labels.py) needs to be executed on the `DEEP_LEARNING` [environment](ENVS/DEEP_LEARNING.yml). This is how to run it:

```sh
python predict_rsa_labels.py Q9UGL1_cif
```

This script only requires a single mandatory argument, which is `input_dir` the name of the input directory. From this `input_dir`, the programme will find the relevant files in the `./OUT` directory. In this example, it is `Q9UGL1_cif`. This script will use a [multilayer perceptron model](RSA_pred_model.h5) [[6](https://www.nature.com/articles/s42003-024-05970-8)] to predict RSA-based Cluster labels and functional scores for each defined binding site. 

## Help and manual

To get help or information about the **LIGYSIS** pipeline, run:

```sh
python ligysis_custom.py -h
```

which will print the manual of the programme:

```
usage: ligysis_custom.py [-h] [--override] [--override_variants] [--variants]
                         [--clust_method CLUST_METHOD]
                         [--clust_dist CLUST_DIST] [--hmm_iters HMM_ITERS]
                         [--cons_thresh_high CONS_THRESH_HIGH]
                         [--cons_thresh_low CONS_THRESH_LOW]
                         [--mes_thresh MES_THRESH]
                         input_dir uniprot_id {pdb,mmcif}

LIGYSIS: a ligand binding site analysis pipeline that clusters ligands,
defines, and characterises binding sites.

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

To get help or information about the RSA-label and scores prediction script, run:

```sh
python predict_rsa_labels.py -h
```

which will print the manual of the programme:

```
usage: predict_rsa_labels.py [-h] input_id

This script predicts RSA cluster labels and calculates RSA-based functional
score (FS)

positional arguments:
  input_dir    This is the Input ID, Job ID or input directory, i.e..
              name of the directory where the binding site table resides.

options:
  -h, --help  show this help message and exit
```

### Optional command line arguments

- `--clust_method` is the clustering algorithm employed to cluster the ligands into binding sites. Average linkage clustering, `average`, is used by default. For other methods, check `scipy.cluster.hierarchy` documentation [here](https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html).

- `--clust_dist` is the distance thredhold at which the clustering tree or dendrogram is cut to obtain ligand clusters or ligand binding sites. The threshold value is `0.5` by default, i.e., on average the ligands within a cluster share half of their binding residues.

- `--hmm_iters` is the number of iterations for the MSA build using `jackHMMER`, default is `3`. More iterations will retrieve more remote homologues, sequences with an evolutionary link further away in time (less sequence similarity).

- `--cons_thresh_low` and `--cons_thresh_high` are the thresholds employed for the Normalised Shenkin divergence score [[2](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009335), [6](https://www.nature.com/articles/s42003-024-05970-8)], ranging 0-100, to determine the classes of <i>conserved</i> and <i>unconserved</i>. Default values are `25` and `75`, so only columns presenting ≤ 25% of the maximum divergence will be classified as <i>conserved</i> and columns with a divergence ≥75% as <i>unconserved</i> or <i>divergent</i>.

- `--mesh_thresh` is the threshold employed to classify columns into <i>enriched</i> or <i>depleted</i> in human missense variation, relative to the average of the other columns. The missense enrichment score, MES, is an odds ratio [[2](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009335), [3](https://www.biorxiv.org/content/10.1101/127050v2), [5](https://www.nature.com/articles/s42003-024-06117-5)], so values > 1 indicate enrichment and < 1 depletion. A values of 1 represents average variation or neutrality within a column. Default threshold is `1.0`.

## Citation

If you use **LIGYSIS** pipeline, please cite:

**Utgés JS**, MacGowan SA, Ives CM, Barton GJ. Classification of likely functional class for ligand binding sites identified from fragment screening. Commun Biol. 2024 Mar 13;7(1):320. doi: [10.1038/s42003-024-05970-8](https://www.nature.com/articles/s42003-024-05970-8). PMID: 38480979; PMCID: PMC10937669.

**Utgés JS** & Barton GJ. Comparative evaluation of methods for the prediction of protein-ligand binding sites, 08 August 2024, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-4849153/v1](https://doi.org/10.21203/rs.3.rs-4849153/v1)

## References

1. Armstrong DR, Berrisford JM, Conroy MJ, Gutmanas A, Anyango S, Choudhary P, Clark AR, Dana JM, Deshpande M, Dunlop R, Gane P, Gáborová R, Gupta D, Haslam P, Koča J, Mak L, Mir S, Mukhopadhyay A, Nadzirin N, Nair S, Paysan-Lafosse T, Pravda L, Sehnal D, Salih O, Smart O, Tolchard J, Varadi M, Svobodova-Vařeková R, Zaki H, Kleywegt GJ, Velankar S. PDBe: improved findability of macromolecular structure data in the PDB. Nucleic Acids Res. 2020 Jan 8;48(D1):D335-D343. doi: [10.1093/nar/gkz990](https://europepmc.org/article/MED/31691821). PMID: 31691821; PMCID: PMC7145656.

2. PDBe-KB consortium. PDBe-KB: a community-driven resource for structural and functional annotations. Nucleic Acids Res. 2020 Jan 8;48(D1):D344-D353. doi: [10.1093/nar/gkz853](https://academic.oup.com/nar/article/48/D1/D344/5580911). PMID: 31584092; PMCID: PMC6943075.

3. Nair S, Váradi M, Nadzirin N, Pravda L, Anyango S, Mir S, Berrisford J, Armstrong D, Gutmanas A, Velankar S. PDBe aggregated API: programmatic access to an integrative knowledge graph of molecular structure data. Bioinformatics. 2021 Nov 5;37(21):3950-3952. doi: [10.1093/bioinformatics/btab424](https://academic.oup.com/bioinformatics/article/37/21/3950/6291664). PMID: 34081107; PMCID: PMC8570819.

4. UniProt Consortium. UniProt: the universal protein knowledgebase in 2021. Nucleic Acids Res. 2021 Jan 8;49(D1):D480-D489. doi: [10.1093/nar/gkaa1100](https://academic.oup.com/nar/article/49/D1/D480/6006196?login=true). PMID: 33237286; PMCID: PMC7778908.

5. Russell RB, Barton GJ. Multiple protein sequence alignment from tertiary structure comparison: assignment of global and residue confidence levels. Proteins. 1992 Oct;14(2):309-23. doi: [10.1002/prot.340140216](https://onlinelibrary.wiley.com/doi/10.1002/prot.340140216). PMID: 1409577.

6. Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO. Maximum allowed solvent accessibilites of residues in proteins. PLoS One. 2013 Nov 21;8(11):e80635. doi: [10.1371/journal.pone.0080635](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0080635#:~:text=It%20is%20defined%20as%20a,1%5D%E2%80%93%5B5%5D.). PMID: 24278298; PMCID: PMC3836772.

7. Kabsch W, Sander C. Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers. 1983 Dec;22(12):2577-637. doi: [10.1002/bip.360221211](https://onlinelibrary.wiley.com/doi/10.1002/bip.360221211). PMID: 6667333.

8. MacGowan SA, Madeira F, Britto-Borges T, Warowny M, Drozdetskiy A, Procter JB, Barton GJ. The Dundee Resource for Sequence Analysis and Structure Prediction. Protein Sci. 2020 Jan;29(1):277-297. doi: [10.1002/pro.3783](https://onlinelibrary.wiley.com/doi/full/10.1002/pro.3783). Epub 2019 Nov 28. PMID: 31710725; PMCID: PMC6933851.

9. Jubb HC, Higueruelo AP, Ochoa-Montaño B, Pitt WR, Ascher DB, Blundell TL. Arpeggio: A Web Server for Calculating and Visualising Interatomic Interactions in Protein Structures. J Mol Biol. 2017 Feb 3;429(3):365-371. doi: [10.1016/j.jmb.2016.12.004](https://www.sciencedirect.com/science/article/pii/S0022283616305332?via%3Dihub). Epub 2016 Dec 10. PMID: 27964945; PMCID: PMC5282402.

10. Virtanen P, Gommers R, Oliphant TE, Haberland M, Reddy T, Cournapeau D, Burovski E, Peterson P, Weckesser W, Bright J, van der Walt SJ, Brett M, Wilson J, Millman KJ, Mayorov N, Nelson ARJ, Jones E, Kern R, Larson E, Carey CJ, Polat İ, Feng Y, Moore EW, VanderPlas J, Laxalde D, Perktold J, Cimrman R, Henriksen I, Quintero EA, Harris CR, Archibald AM, Ribeiro AH, Pedregosa F, van Mulbregt P; SciPy 1.0 Contributors. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nat Methods. 2020 Mar;17(3):261-272. doi: 10.1038/s41592-019-0686-2. Epub 2020 Feb 3. Erratum in: Nat Methods. 2020 Mar;17(3):352. doi: [10.1038/s41592-020-0772-5](https://www.nature.com/articles/s41592-019-0686-2). PMID: 32015543; PMCID: PMC7056644.

11. Pettersen EF, Goddard TD, Huang CC, Meng EC, Couch GS, Croll TI, Morris JH, Ferrin TE. UCSF ChimeraX: Structure visualization for researchers, educators, and developers. Protein Sci. 2021 Jan;30(1):70-82. doi: [10.1002/pro.3943](https://onlinelibrary.wiley.com/doi/10.1002/pro.3943). Epub 2020 Oct 22. PMID: 32881101; PMCID: PMC7737788.

12. Eddy SR. Profile hidden Markov models. Bioinformatics. 1998;14(9):755-63. doi: [10.1093/bioinformatics/14.9.755](https://academic.oup.com/bioinformatics/article/14/9/755/259550). PMID: 9918945.

13. Shenkin PS, Erman B, Mastrandrea LD. Information-theoretical entropy as a measure of sequence variability.
Proteins. 1991; 11(4):297–313. Epub 1991/01/01. [https://doi.org/10.1002/prot.340110408](https://doi.org/10.1002/prot.340110408)
PMID: 1758884.

14. **Utgés JS**, Tsenkov MI, Dietrich NJM, MacGowan SA, Barton GJ. Ankyrin repeats in context with human population variation. PLoS Comput Biol. 2021 Aug 24;17(8):e1009335. doi: [10.1371/journal.pcbi.1009335](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009335). PMID: 34428215; PMCID: PMC8415598.
 
15. MacGowan, SA, Madeira, F, Britto-Borges, T, Schmittner, MS, Cole, C, & Barton, GJ (2017). Human missense variation is constrained by domain structure and highlights functional and pathogenic residues. bioRxiv, 127050. [https://doi.org/10.1101/127050](https://www.biorxiv.org/content/10.1101/127050v2).

16. MacGowan SA, Madeira F, Britto-Borges T, Barton GJ. A unified analysis of evolutionary and population constraint in protein domains highlights structural features and pathogenic sites. Commun Biol. 2024 Apr 11;7(1):447. doi: [10.1038/s42003-024-06117-5](https://www.nature.com/articles/s42003-024-06117-5). PMID: 38605212; PMCID: PMC11009406.
   
17. **Utgés JS**, MacGowan SA, Ives CM, Barton GJ. Classification of likely functional class for ligand binding sites identified from fragment screening. Commun Biol. 2024 Mar 13;7(1):320. doi: [10.1038/s42003-024-05970-8](https://www.nature.com/articles/s42003-024-05970-8). PMID: 38480979; PMCID: PMC10937669.

18. Cock PJ, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJ. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009 Jun 1;25(11):1422-3. doi: [10.1093/bioinformatics/btp163](https://academic.oup.com/bioinformatics/article/25/11/1422/330687). Epub 2009 Mar 20. PMID: 19304878; PMCID: PMC2682512.
   
19. Chollet, F., & others. (2015). Keras. Retrieved from [https://keras.io](https://keras.io).

20. Harris CR, Millman KJ, van der Walt SJ, Gommers R, Virtanen P, Cournapeau D, Wieser E, Taylor J, Berg S, Smith NJ, Kern R, Picus M, Hoyer S, van Kerkwijk MH, Brett M, Haldane A, Del Río JF, Wiebe M, Peterson P, Gérard-Marchant P, Sheppard K, Reddy T, Weckesser W, Abbasi H, Gohlke C, Oliphant TE. Array programming with NumPy. Nature. 2020 Sep;585(7825):357-362. doi: [10.1038/s41586-020-2649-2](https://www.nature.com/articles/s41586-020-2649-2). Epub 2020 Sep 16. PMID: 32939066; PMCID: PMC7759461.

21. McKinney, W. (2010). Data Structures for Statistical Computing in Python. In S. van der Walt & J. Millman (Eds.), Proceedings of the 9th Python in Science Conference (pp. 56–61). doi: [10.25080/Majora-92bf1922-00a](https://doi.org/10.25080/Majora-92bf1922-00a)

22. Pedregosa, F., Varoquaux, G., Gramfort, A., Michel, V., Thirion, B., Grisel, O., Blondel, M., Prettenhofer, P., Weiss, R., Dubourg, V., Vanderplas, J., Passos, A., Cournapeau, D., Brucher, M., Perrot, M., & Duchesnay, É. (2011). Scikit-learn: Machine Learning in Python. Journal of Machine Learning Research, 12(85), 2825–2830. [https://www.jmlr.org/papers/v12/pedregosa11a.html](https://www.jmlr.org/papers/v12/pedregosa11a.html).

23. Martín Abadi, Ashish Agarwal, Paul Barham, Eugene Brevdo, Zhifeng Chen, Craig Citro, Greg S. Corrado, Andy Davis, Jeffrey Dean, Matthieu Devin, Sanjay Ghemawat, Ian Goodfellow, Andrew Harp, Geoffrey Irving, Michael Isard, Rafal Jozefowicz, Yangqing Jia, Lukasz Kaiser, Manjunath Kudlur, Josh Levenberg, Dan Mané, Mike Schuster, Rajat Monga, Sherry Moore, Derek Murray, Chris Olah, Jonathon Shlens, Benoit Steiner, Ilya Sutskever, Kunal Talwar, Paul Tucker, Vincent Vanhoucke, Vijay Vasudevan, Fernanda Viégas, Oriol Vinyals, Pete Warden, Martin Wattenberg, Martin Wicke, Yuan Yu, and Xiaoqiang Zheng. TensorFlow: Large-scale machine learning on heterogeneous systems, 2015. [Paper here](https://static.googleusercontent.com/media/research.google.com/en//pubs/archive/45166.pdf) Software available from [tensorflow.org](tensorflow.org).

24. Karczewski KJ, Francioli LC, Tiao G, Cummings BB, Alföldi J, Wang Q, Collins RL, Laricchia KM, Ganna A, Birnbaum DP, Gauthier LD, Brand H, Solomonson M, Watts NA, Rhodes D, Singer-Berk M, England EM, Seaby EG, Kosmicki JA, Walters RK, Tashman K, Farjoun Y, Banks E, Poterba T, Wang A, Seed C, Whiffin N, Chong JX, Samocha KE, Pierce-Hoffman E, Zappala Z, O'Donnell-Luria AH, Minikel EV, Weisburd B, Lek M, Ware JS, Vittal C, Armean IM, Bergelson L, Cibulskis K, Connolly KM, Covarrubias M, Donnelly S, Ferriera S, Gabriel S, Gentry J, Gupta N, Jeandet T, Kaplan D, Llanwarne C, Munshi R, Novod S, Petrillo N, Roazen D, Ruano-Rubio V, Saltzman A, Schleicher M, Soto J, Tibbetts K, Tolonen C, Wade G, Talkowski ME; Genome Aggregation Database Consortium; Neale BM, Daly MJ, MacArthur DG. The mutational constraint spectrum quantified from variation in 141,456 humans. Nature. 2020 May;581(7809):434-443. doi: 10.1038/s41586-020-2308-7. Epub 2020 May 27. Erratum in: Nature. 2021 Feb;590(7846):E53. doi: 10.1038/s41586-020-03174-8. Erratum in: Nature. 2021 Sep;597(7874):E3-E4. doi: [10.1038/s41586-021-03758-y](https://www.nature.com/articles/s41586-020-2308-7). PMID: 32461654; PMCID: PMC7334197.

25. McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P, Cunningham F. The Ensembl Variant Effect Predictor. Genome Biol. 2016 Jun 6;17(1):122. doi: [10.1186/s13059-016-0974-4](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0974-4). PMID: 27268795; PMCID: PMC4893825.
