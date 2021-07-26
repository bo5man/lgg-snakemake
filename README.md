<!--
<p align="center">
  <img width="100%" height="100%" src="https://github.com/tgac-vumc/lgg-snakemake/blob/master/DAG_all.svg">
</p>
-->
## Installation

For the installation of this pipeline any Python install compatable Conda is required.

The pipeline itself will run on Python 3.9.5, snakemake 6.4.1 and R 4.0.5. For exact dependencies view `lgg-snakemake.yaml`.

### Using Conda/Mamba

for easy installation you need (Mini)Conda.

Miniconda installation from folder where you want to install Miniconda:

```
cd </path/to/files/dir/>
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

follow the instructions of the installation process, give the location where you want Miniconda to be installed and answer YES to add Miniconda to your path.

go to the directory where the analysis needs to be performed

The pipeline `lgg-snakemake` makes use of a GitHub repository, `mnp-training`, of supplimentary material of the reference
```
Capper, D., Jones, D., Sill, M. et al. DNA methylation-based classification of central nervous system tumours. Nature 555, 469–474 (2018). https://doi.org/10.1038/nature26000
```
Please see `https://doi.org/10.1038/nature26000` for the full article and `https://github.com/mwsill/mnp_training/` for the Github repository.

The `conumee` Bioconductor package in the latter repository is applied to perform copy number variation analysis. To get the reference objects stored in `./CNV_data`, Git large file storage needs to be installed before cloning the repository `mnp_training`.

Follow the installation instructions at `https://github.com/git-lfs/git-lfs/wiki/Installation`

After installation of `Git Large File Storage`, clone the `mnp_training` and `lgg-snakemake` repositories in the analysis directory.

```
cd </path/to/analysis/dir>
git clone https://github.com/mwsill/mnp_training/
git clone https://github.com/tgac-vumc/lgg-snakemake/
cd lgg-snakemake
```

install Mamba as drop-in replacement for Conda with Mamba's improved installation-performance:

```
conda install -c conda-forge mamba
```

create the environment using Mamba:

```
mamba env create --name lgg-snakemake --file lgg-snakemake.yaml
```

activate the environment by:

```
conda activate lgg-snakemake
```

## Preparing analysis

### Prepare the data

go to analysis dir and prepare analysis by copy or create links to .idat files:

```
cd </path/to/analysis/dir>

mkdir data
cd data
mkdir allsamples
mkdir fullcohort
```

to link all files from a folder:

```
for file in <path/to/idat/files>/*.idat
do ln -s $file
done
```

The `fullcohort`-directory should contain all idat-files of the full cohort.
The `allsamples`-directory should contain all idat-files of the full cohort AND all idat-files of related patients, for this project this includes the `inhouse`-data.

Link the `fullcohort`-data to the `fullcohort`-directory and the `allsamples`-directory.
Link the `inhouse`-data to the `allsamples`-directory.

The directory tree of `data` should look something like this

```
data
├── allsamples
│   ├─── <sentrix1fullcohort>_Grn.idat
│   ├─── <sentrix1fullcohort>_Red.idat
│   ├─── <sentrix2fullcohort>_Grn.idat
│   ├─── <sentrix2fullcohort>_Red.idat
│   ├─── ...
│   ├─── <sentrix1inhouse>_Grn.idat
│   ├─── <sentrix1inhouse>_Red.idat
│   ├─── <sentrix2inhouse>_Grn.idat
│   ├─── <sentrix2inhouse>_Red.idat
│   └─── ...
└── fullcohort
    ├─── <sentrix1fullcohort>_Grn.idat
    ├─── <sentrix1fullcohort>_Red.idat
    ├─── <sentrix2fullcohort>_Grn.idat
    ├─── <sentrix2fullcohort>_Red.idat
    └─── ...
```
### Prepare the snakemake settings

Open the configuration file `config.yaml` to check the settings that snakemake will use and change according to your needs.

Fill all settings in accordingly. The meaning of each setting is explained in the `config.yaml` file.

All paths in the `config.yaml` must be documented relative to the directory that contains the `Snakefile`.

## Running analysis

Make sure that snakemake is able to find the excecutive file `Snakefile` by performing a dry-run:

```
cd ../lgg-snakemake
snakemake -n
```

Check the rules that are planned to be performed, conform the rule-graph.

A visualization of the order of rules to be performed can be viewed by running the following command and opening the DAG-file

```
snakemake --forceall --rulegraph | dot -Tsvg > DAG.svg
```

When ready, run the analysis

```
snakemake
```

Useful snakemake options

`-j , --cores, --jobs` : Use at most N cores in parallel (default: 1). If N is omitted, the limit is set to the number of available cores. Example `-j5` or `--cores=5`

`-n , --dryrun` : Do not execute anything. but show rules which are planned to be performed.

`-k , --keep-going` : Go on with independent jobs if a job fails.

`-f , --force` : Force the execution of the selected target or the first rule regardless of already created output. Example `-f analysis_probes`

`-R , --forcerun` : Force the re-execution or creation of the given rules or files. Use this option if you changed a rule and want to have all its output in your workflow updated. Example `-R create_clinicalDF`

`-U , --until` : Runs the pipeline until it reaches the specified rules or files. Only runs jobs that are dependencies of the specified rule or files, does not run sibling DAGs. Example `-U analysis_probes`

for all options go to https://snakemake.readthedocs.io/en/v6.4.1/executing/cli.html#all-options
