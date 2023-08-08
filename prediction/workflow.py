### SETTINGS ###

# input files
studypth = ""
studies = ["study1", "study2"]
individuals = ""

# program paths
ldak = ".././ldak5.2.linux"
plink = "~/miniforge3/envs/plink/bin/plink"
julia = "~/miniforge3/envs/ldak/bin/julia"
process_sumstats = "../../process_sumstats/ProcessSumstats.jl"

# settings
nthreads = 12


### CODE ###
from os import makedirs
from gwf import Workflow
gwf = Workflow(defaults = {"account": "perjektet"})

def bim(s): return f'{s}.bim'
def pl(s): return [f'{s}.{ext}' for ext in ['bim', 'bed', 'fam']]

def newdir(s):
    makedirs(s, exist_ok = True)
    return(s)


### Download files that are used for all summary stats pipelines

d = "files"

# Download and extract the first 64 SNP annotations for the BLD-LDAK Model. (See [BLD-LDAK Annotations](http://dougspeed.com/bldldak/)).

d_bld = newdir(f'{d}/bld'); s_bld = f'{d_bld}/bld'
bld_files = [f'{s_bld}{i}' for i in range(0, 65)]
gwf.target(
    f'download_snp_annotations_{d}',
    inputs = [],
    outputs = bld_files
    ) << f"""
        wget -O {s_bld}.zip https://genetics.ghpc.au.dk/doug/bld.zip
        unzip -d {d_bld} {s_bld}.zip
    """


# Download and extract a non-Finnish European reference panel constructed from 1000 Genomes Project data.

d_ref = newdir(f'{d}/ref'); 
s_ref_orig = f'{d_ref}/ref'
gwf.target(
    f'download_reference_{d}',
    inputs = [],
    outputs = pl(s_ref_orig)
    ) << f"""
        wget -O {s_ref_orig}.tgz https://genetics.ghpc.au.dk/doug/ref.tgz
        tar -xzvf {s_ref_orig}.tgz -C {d_ref}
    """



### Pipeline for each individual sumstat file

def workflows(studyname, d):

    # redefine folders for individual sumstats
    d_ref = newdir(f'{d}/ref') 

# Format and clean the summary statistics. (Using [ProcessSumstats](https://github.com/ymer/process_sumstats))

    d_sumstats = newdir(f'{d}/sumstats')
    sumstats = f'{d_sumstats}/{studyname}'
    snplist = f'{sumstats}.snplist'
    gwf.target(
        f'process_sumstats_{d}',
        memory = '64g',
        inputs = pl(s_ref_orig),
        outputs = [snplist, sumstats]
        ) << f"""
            {julia} {process_sumstats} --in {studypth}/{studyname} --outdir {d_sumstats} --write-snplist --filter-reference {s_ref_orig}.bim
        """


# Filter the reference genome to only include the SNPs remaining after the above sumstats cleaning.

    s_ref = f'{d_ref}/ref_filter'
    gwf.target(
        f'filter_ref_{d}',
        memory = '64g',
        inputs = pl(s_ref_orig) + [snplist],
        outputs = pl(s_ref)
        ) << f"""
            {plink} --bfile {s_ref_orig} --extract {snplist} --make-bed --out {s_ref}
        """


# Download details of long-range linkage disquilibrium regions and identify which reference panel SNPs they contain.(See [High-LD Regions](http://dougspeed.com/high-ld-regions/)).

    d_highld = newdir(f'{d}/highld')
    highld = f'{d_highld}/highld.txt'
    highld_predictors = f'{d_highld}/genes.predictors.used'
    gwf.target(
        f'get_high_ld_regions_{d}',
        memory = '64g',
        inputs = pl(s_ref), 
        outputs = [highld_predictors]
        ) << f"""
            wget -O {highld} http://dougspeed.com/wp-content/uploads/highld.txt
            {ldak} --cut-genes {d_highld} --bfile {s_ref} --genefile {highld} 
        """


# Compute LDAK weightings. (See [LDAK Weightings](http://dougspeed.com/calculate-weightings/))

    d_sections = newdir(f'{d}/sections')
    weights_ch = []
    for ch in range(21, 23):
        d_section = newdir(f'{d_sections}/s{ch}')
        weights_short = f'{d_section}/weights.short'
        weights_ch.append(weights_short)
        gwf.target(
            f'compute_weightings_{ch}_{d}',
            inputs = pl(s_ref) + [snplist],
            outputs = [weights_short],
            memory = '64g'
        ) << f"""
            {ldak} --cut-weights {d_section} --bfile {s_ref} --extract {snplist} --chr {ch}
            {ldak} --calc-weights-all {d_section} --bfile {s_ref} --extract {snplist} --chr {ch}
            """

    weights = f'{d}/weights'
    cat_command = 'cat ' + ' '.join(weights_ch) + f' > {weights}'
    gwf.target(
        f'join_weights_{d}',
        memory = '64g',
        inputs = weights_ch,
        outputs = [weights],
        ) << f"""
        cat {' '.join(weights_ch)} > {weights}
        """


# Calculate the tagging file and heritability matrix assuming the BLD-LDAK Model. (See [Calculate Taggings](http://dougspeed.com/calculate-taggings/)).

    d_bldldak = newdir(f'{d}/bld.ldak/')
    s_bldldak = f'{d_bldldak}/bld.ldak'
    matrix = f'{s_bldldak}.matrix'
    tagging = f'{s_bldldak}.tagging'
    gwf.target(
        f'calc_tagging_{d}',
        memory = '64g', cores = nthreads, walltime = '24:00:00',
        inputs = pl(s_ref) + [snplist],
        outputs = [matrix, tagging]
        ) << f"""
            {ldak} --calc-tagging {s_bldldak} --bfile {s_ref} --extract {snplist} --ignore-weights YES --power -.25 --annotation-number 64 --annotation-prefix {s_bld} --window-cm 1 --save-matrix YES --max-threads {nthreads}
        """


# Estimate the heritability contributed by each SNP. (See [SNP Heritability](http://dougspeed.com/snp-heritability/)).

    ind_hers = f'{s_bldldak}.ind.hers'
    gwf.target(
        f'estimate_heritability_{d}',
        memory = '64g',
        inputs = [matrix, tagging],
        outputs = [ind_hers]
        ) << f"""
            {ldak} --sum-hers {s_bldldak} --tagfile {tagging} --summary {sumstats} --matrix {matrix}
        """


# Calculate predictor-predictor correlations. (See [MegaPRS](http://dougspeed.com/megaprs/))

    d_cors = newdir(f'{d}/cors/')
    s_cors = f'{d_cors}/cors'
    cors_files = [f'{s_cors}.cors.{ext}' for ext in ['bim', 'bin', 'noise', 'root']]
    gwf.target(
        f'calc_cors_{d}',
        memory = '64g', cores = nthreads,
        inputs = pl(s_ref) + [snplist],
        outputs = cors_files
        ) << f"""
            {ldak} --calc-cors {s_cors} --bfile {s_ref} --window-cm 3 --extract {snplist} --max-threads {nthreads}
        """


# Construct the prediction model (see [MegaPRS](http://dougspeed.com/megaprs/)).

    d_bayesr = newdir(f'{d}/bayesr')
    s_bayesr = f'{d_bayesr}/bayesr'
    effects = f'{s_bayesr}.effects'
    gwf.target(
        f'mega_prs_{d}',
        memory = '64g', cores = nthreads,
        inputs = [snplist, ind_hers, sumstats, highld_predictors] + cors_files,
        outputs = [effects]
        ) << f"""
            {ldak} --mega-prs {s_bayesr} --model bayesr --ind-hers {ind_hers} --summary {sumstats} --cors {s_cors} --cv-proportion .1 --high-LD {highld_predictors} --window-cm 1 --extract {snplist} --max-threads {nthreads}
        """


# Calculate scores

    s_scores = f'{d}/scores'
    scores = f'{s_scores}.profile'
    gwf.target(
        f'calc_prs_{d}',
        memory = '64g', cores = nthreads,
        inputs = [effects],
        outputs = [scores]
        ) << f"""
            {ldak} --calc-scores {s_scores} --scorefile {effects} --bfile {individuals} --power 0    
    """

for studyname in studies:
    studyname = newdir(studyname)
    workflows(studyname, studyname)

