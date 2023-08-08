A pipeline for running [LDAK]([http://dougspeed.com/]) methods on a computing cluster.

Uses [gwf](https://gwf.app/)

How to use:
- Install conda requirements
´´´
conda create --name plink
conda activate plink
conda install -c bioconda plink
conda create --name ldak
conda activate ldak
conda install gwf
conda install -c conda-forge julia
´´´

- Get ProcessSumstats from [https://github.com/ymer/process_sumstats](github.com/ymer/process_sumstats)
´´´
git clone https://github.com/ymer/process_sumstats
´´´

- Install julia packages.
Go to the process_sumstats directory, and
´´´
julia julia_install.jl
´´´

- Make the ldak file executable.
Or download the newest version from [dougspeed.com](http://dougspeed.com/)

- Go to the directory of the pipeline you are interested in running. Eg 'prediction'.

- Open workflow.py and change the settings.
The paths for Julia and plink can be found with the conda command `which`

- Run the pipeline with `gwf run`


