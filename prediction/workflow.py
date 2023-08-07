### SETTINGS
ldak = ".././ldak5.2.linux"
studyname = "ACC_DIURNAL_INACT_RAW.csv"
studypth = "../../sumstats/"
julia = "~/miniforge3/envs/julia/bin/julia"
d = "files"
process_sumstats = "../../process_sumstats/ProcessSumstats.jl"
nthreads = 8

### CODE
from os import makedirs
from gwf import Workflow
gwf = Workflow(defaults = {"account": "perjektet"})

def bim(s): return f'{s}.bim'
def pl(s): return [f'{s}.{ext}' for ext in ['bim', 'bed', 'fam']]

def newdir(s):
    makedirs(s, exist_ok = True)
    return(s)

d_bld = newdir(f'{d}/bld'); s_bld = f'{d_bld}/bld'
bld_files = [f'{s_bld}{i}' for i in range(0, 65)]
gwf.target(
    f'download_snp_annotations',
    inputs = [],
    outputs = bld_files
    ) << f"""
        wget -O {s_bld}.zip https://genetics.ghpc.au.dk/doug/bld.zip
        unzip -d {d_bld} {s_bld}.zip
    """

d_ref = newdir(f'{d}/ref'); s_ref = f'{d_ref}/ref'
gwf.target(
    f'download_reference',
    inputs = [],
    outputs = pl(s_ref)
    ) << f"""
        wget -O {s_ref}.tgz https://genetics.ghpc.au.dk/doug/ref.tgz
        tar -xzvf {s_ref}.tgz -C {d_ref}
    """

d_highld = newdir(f'{d}/highld')
highld = f'{d_highld}/highld.txt'
genes_predictors = f'{d_highld}/genes.predictors.used'
gwf.target(
    f'get_high_ld_regions',
    inputs = pl(s_ref), 
    outputs = [genes_predictors]
    ) << f"""
        wget -O {highld} http://dougspeed.com/wp-content/uploads/highld.txt
        {ldak} --cut-genes {d_highld} --bfile {s_ref} --genefile {highld} 
    """

d_sumstats = newdir(f'{d}/sumstats')
sumstats = f'{d_sumstats}/{studyname}'
snplist = f'{sumstats}.snplist'
gwf.target(
    f'process_sumstats',
    memory = '64g',
    inputs = pl(s_ref),
    outputs = [snplist, sumstats]
    ) << f"""
        {julia} {process_sumstats} --in {studypth}/{studyname} --outdir {d_sumstats} --write-snplist --filter-reference {s_ref}.bim
    """

d_sections = newdir(f'{d}/sections')
weights_ch = []
for ch in range(21, 23):
    d_section = newdir(f'{d_sections}/s{ch}')
    weights_short = f'{d_section}/weights.short'
    weights_ch.append(weights_short)
    gwf.target(
        f'compute_weightings_{ch}',
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
    f'join_weights',
    memory = '64g',
    inputs = weights_ch,
    outputs = [weights],
    ) << f"""
    cat {' '.join(weights_ch)} > {weights}
    """


d_bldldak = newdir(f'{d}/bld.ldak/')
s_bldldak = f'{d_bldldak}/bld.ldak'
matrix = f'{s_bldldak}.matrix'
tagging = f'{s_bldldak}.tagging'
gwf.target(
    f'calc_tagging',
    memory = '64g', cores = nthreads, walltime = '24:00:00',
    inputs = pl(s_ref) + [snplist],
    outputs = [matrix, tagging]
    ) << f"""
        {ldak} --calc-tagging {s_bldldak} --bfile {s_ref} --extract {snplist} --ignore-weights YES --power -.25 --annotation-number 64 --annotation-prefix {s_bld} --window-cm 1 --save-matrix YES --max-threads {nthreads}
    """

ind_hers = f'{s_bldldak}.ind.hers'
gwf.target(
    f'estimate_heritability',
    memory = '64g',
    inputs = [matrix, tagging],
    outputs = [ind_hers]
    ) << f"""
        {ldak} --sum-hers {s_bldldak} --tagfile {tagging} --summary {sumstats} --matrix {matrix} --check-sums NO
    """


d_cors = newdir(f'{d}/cors/')
s_cors = f'{d_cors}/cors'
cors_files = [f'{s_cors}.cors.{ext}' for ext in ['bim', 'bin', 'noise', 'root']]
gwf.target(
    f'calc_cors',
    memory = '64g', cores = nthreads,
    inputs = pl(s_ref) + [snplist],
    outputs = cors_files
    ) << f"""
        {ldak} --calc-cors {s_cors} --bfile {s_ref} --window-cm 3 --extract {snplist} --max-threads {nthreads}
    """


gwf.target(
    f'mega_prs',
    memory = '64g', cores = 1,
    inputs = [snplist, ind_hers, sumstats] + cors_files,
    outputs = []
    ) << f"""
        {ldak} --mega-prs bayesr --model bayesr --ind-hers {ind_hers} --summary {sumstats} --cors {s_cors} --cv-proportion .1 --high-LD {highld} --window-cm 1 --extract {snplist} --max-threads {nthreads}
    """



