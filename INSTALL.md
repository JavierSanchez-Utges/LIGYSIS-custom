# Installation

## Installation of DSSP

DSSP is incompatible with all other environments, and so must go on its environment of its own. You can install locally as well, all we will need is the location of the executable. The version of the libboost library must be specified to be this one, otherwise, dssp will not run.

```
conda create -n DSSP salilab::dssp=3.0.0 libboost=1.73.0
```

## Installation of HMMER

Installing HMMER on its own environment. Just easier...

```
conda create -n HMMER hmmer=3.4
```

## Installation of STAMP
The following instructions are to install STAMP. For more information refer to the [STAMP installation instructions](https://www.compbio.dundee.ac.uk/downloads/stamp/INSTALL).

```
# download STAMP
curl -O https://www.compbio.dundee.ac.uk/downloads/stamp/stamp.4.4.2.tar.gz

# decompress STAMP
tar -xzvf stamp.4.4.2.tar.gz

# change directory
cd stamp.4.4.2
```
To install STAMP, run the BUILD script in this directory using:
```
# building STAMP
./BUILD <system-type>
```
where \<system-type\> is one of:

- linux
- osx 
- dec
- sgi
- sun

The executables will be installed in bin/\<system-type\>/.

For more information refer to the [STAMP manual](https://www.compbio.dundee.ac.uk/manuals/stamp.4.4/stamp.html)

## Installation of LIGYSIS<sub>CUSTOM</sub>

The first step to install **LIGYSIS<sub>CUSTOM</sub>** is to Git Clone the repository.

```
# git clone LIGYSIS from repository
git clone -b revamped https://github.com/JavierSanchez-Utges/fragsys_custom.git
```

### Installation of pdbe-arpeggio

```
# change directory to environments directory
cd fragsys_custom/ENV

# install pdbe-arpeggio environment
conda create -n ARPEGGIO python=3.9 gemmi openbabel biopython -c conda-forge

# activating pdbe-arpeggio environment
conda activate ARPEGGIO

# install pdbe-arpeggio
pip install pdbe-arpeggio 

# test pdbe-arpeggio with help function
pdbe-arpeggio -h
```

The next step is to install the three Conda environments needed to run the pipeline and analyse the results. This can be done with Conda using the .yml files in the [`ENVS`](ENVS/) directory.

```
# install other environments

# install deep_learning environment
conda env create -f deep_learning_env.yml

# install arpeggio environment
conda create -n arpeggio-env python=3.9 gemmi openbabel biopython -c conda-forge

# activating arpeggio environment
conda activate arpeggio-env

# install pdbe-arpeggio
pip install pdbe-arpeggio 

# test pdbe-arpeggio with help function
pdbe-arpeggio -h
```

## Installation of VarAlign

The following instructions are to install VarAlign. Fore more information refer to the [VarAlign repository](https://github.com/bartongroup/SM_VarAlign/tree/JSU_branch).

```
# change directory to FRAGSYS envs directory
cd ../fragsys_custom/ENVS/

# install varalign environment
conda env create -f varalign_env.yml

# VarAlign installation (Grabbed from URL)

# change directory to main working directory
cd ../..

# git clone specific branch of VarAlign from repository
git clone -b JSU_branch https://github.com/bartongroup/SM_VarAlign.git

# change directory to VarAlign directory
cd SM_VarAlign

# activate varalign_env environment
conda activate varalign_env

# install VarAlign
pip install .
```

## Installation of ProIntVar

The following instructions are to install ProIntVar. Fore more information refer to the [ProIntVar repository](https://github.com/bartongroup/ProIntVar/tree/JSU_branch).

```
# ProIntVar installation (Grabbed from URL)

# change directory to main working directory
cd ..

# git clone specific branch of ProIntVar from repository
git clone -b JSU_branch https://github.com/bartongroup/ProIntVar.git

### change directory to ProIntVar directory
cd ProIntVar

# pip install ProIntVar dependencies
pip install -r requirements.txt

#then
python setup.py install
```

