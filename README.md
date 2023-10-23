A pipeline for running [LDAK]([http://dougspeed.com/]) methods on a computing cluster.

Uses [gwf](https://gwf.app/)

How to use:
- Install conda requirements
```
conda create --name plink
conda activate plink
conda install -c bioconda plink
conda create --name ldak
conda activate ldak
conda install gwf
conda install -c conda-forge julia
```

- Get ProcessSumstats from [https://github.com/ymer/process_sumstats](github.com/ymer/process_sumstats)
```
git clone https://github.com/ymer/process_sumstats
```

- Install julia packages.
Go to the process_sumstats directory, and

```
julia julia_install.jl
```

- Make the ldak file executable. Or download the newest version from [dougspeed.com](http://dougspeed.com/)

- Download and untar a reference data set from LDAK. [https://genetics.ghpc.au.dk/doug/ref.tgz](https://genetics.ghpc.au.dk/doug/ref.tgz).  
This has SNP position as identifier, and genetic distances. It is called `ref_snp`.

- Download a reference data set from MAGMA. [https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip](https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip).  
This has Rsid as identifier. It is called `ref_rsid`.

- Go to the directory of the pipeline you are interested in running. Eg 'prediction'.

- Open workflow.py and change the settings.  
The paths for Julia and plink can be found with the conda command 'which'.  
Add the paths to the just downloaded `ref_snp` and `ref_rsid` files, as well as the study you want to analyze
(`individuals`) 

- Open studies.csv and add the studies you are interested in running ldak on.  
Delete the example on the second line.  
Note that some sumstat files may require special commands for ProcessSumstats. These can all be added in the process_sumstats_arguments column.

- Run the pipeline with `gwf run`

- Look at the progress of the pipeline with `gwf status`

For great justice.

