# LIGYSIS pipeline for user job submission

This repository contains a customised version of our original ligand site analysis [**LIGYSIS**](https://github.com/JavierSanchez-Utges/ligysis_dev) pipeline, **LIGYSIS<sub>CUSTOM</sub>**, used for the analysis of protein-ligand complexes submitted to our **LIGYSIS** web server found [here](https://www.compbio.dundee.ac.uk/ligysis/). The code for the web server can be found [here](https://github.com/JavierSanchez-Utges/ligysis_flask).

## Running LIGYSIS<sub>CUSTOM</sub>

This is a command line programme, which can be executed like this:

```sh
python fragsys_custom.py IN/CMTR1 Q8N1G2 pdb
```

The programme has three mandatory arguments:
- `input_dir` is the input directory containing the set of structures to be analysed in either [PDB](https://www.wwpdb.org/documentation/file-format) (<i>.pdb</i>, <i>.ent</i>) or [mmCIF](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/beginner%E2%80%99s-guide-to-pdbx-mmcif) (<i>.cif</i>) formats. In this example:`IN/CMTR1`.

- `uniprot_id` is the corresponding [UniProt](https://www.uniprot.org/) accession identifier of the structures within `input_dir`. In this example: [Q8N1G2](https://www.uniprot.org/uniprotkb/Q8N1G2/entry), which corresponds to human Cap-specific mRNA (nucleoside-2'-O-)-methyltransferase 1.

- Finally, `struc_fmt`, indicating whether the structures are in `pdb` or `mmcif` format. Only these two formats are selected, and the programme will not run properly unless a structure format is required.

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
