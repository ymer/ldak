### SETTINGS
ldak = ".././ldak5.2.linux"
studyname = "ACC_DIURNAL_INACT_RAW.csv"
studypth = "../../sumstats/"
julia = "~/miniforge3/envs/julia/bin/julia"
d = "files"
process_sumstats = "../../process_sumstats/ProcessSumstats.jl"


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
gwf.target(
    f'get_high_ld_regions',
    inputs = [], 
    outputs = [highld]
    ) << f"""
        wget -O {highld} http://dougspeed.com/wp-content/uploads/highld.txt
    """

genes_predictors = f'{d_highld}/genes.predictors.used'
gwf.target(
    f'process_high_ld_regions',
    memory = '64g',
    inputs = [highld] + pl(s_ref),
    outputs = [genes_predictors]
    ) << f"""
        {ldak} --cut-genes {d_highld} --bfile {s_ref} --genefile {highld} 
    """

d_sumstats = f'{d}/sumstats'
makedirs(d_sumstats, exist_ok = True)
sumstats = f'{d_sumstats}/{studyname}'
snplist = f'{d_sumstats}/{studyname}.snplist'
gwf.target(
    f'process_sumstats',
    memory = '64g',
    inputs = pl(s_ref),
    outputs = [snplist, sumstats],
    cores = 8
    ) << f"""
        {julia} {process_sumstats} --in {studypth}/{studyname} --outdir {d_sumstats} --write-snplist --filter-reference {s_ref}.bim
    """

d_sections = f'{d}/sections'
weights_chr = []
makedirs(d_sections, exist_ok=True)
for ch in range(21, 23):
    section = f'{d_sections}/s{ch}'
    weights_short = f'{section}/weights.short'
    weights_chr.append(weights_short)
    makedirs(section, exist_ok=True)
    gwf.target(
        f'compute_weightings_{ch}',
        inputs = pl(s_ref) + [snplist],
        outputs = [weights_short],
        memory = '64g'
    ) << f"""
        {ldak} --cut-weights {section} --bfile {s_ref} --extract {snplist} --chr {ch}
        {ldak} --calc-weights-all {section} --bfile {s_ref} --extract {snplist} --chr {ch}
        """

weights = f'{d}/weights'
gwf.target(
    f'join_weights',
    memory = '64g',
    inputs = weights_chr,
    outputs = [weights],
    ) << f"""
        cat {d_sections}/s{{21..22}}/weights.short > {weights}
    """


d_bldldak = f'{d}/bld.ldak/'
makedirs(d_bldldak, exist_ok = True)
stub_bldldak = f'{d_bldldak}/bld.ldak'
matrix = f'{stub_bldldak}.matrix'
tagging = f'{stub_bldldak}.tagging'

gwf.target(
    f'calc_tagging',
    memory = '64g', cores = 8, walltime = '24:00:00',
    inputs = pl(s_ref) + [snplist],
    outputs = [matrix, tagging]
    ) << f"""
        {ldak} --calc-tagging {stub_bldldak} --bfile {s_ref} --extract {snplist} --ignore-weights YES --power -.25 --annotation-number 64 --annotation-prefix {s_bld} --window-cm 1 --save-matrix YES
    """

ind_hers = f'{stub_bldldak}.ind.hers'
gwf.target(
    f'estimate_heritability',
    memory = '64g',
    inputs = [matrix, tagging],
    outputs = [ind_hers]
    ) << f"""
        {ldak} --sum-hers {stub_bldldak} --tagfile {tagging} --summary {sumstats} --matrix {matrix}
    """


d_cors = f'{d}/cors/'
makedirs(d_cors, exist_ok=True)
stub_cors = f'{d_cors}/cors'
cors_bim = f'{stub_cors}.cors.bim'
gwf.target(
    f'calc_cors',
    memory = '64g',
    inputs = pl(s_ref) + [snplist],
    outputs = [cors_bim]
    ) << f"""
        {ldak} --calc-cors {stub_cors} --bfile {s_ref} --window-cm 3 --extract {snplist}
    """


gwf.target(
    f'mega_prs',
    memory = '64g',
    inputs = [snplist, ind_hers, sumstats],
    outputs = []
    ) << f"""
        {ldak} --mega-prs bayesr --model bayesr --ind-hers {ind_hers} --summary {sumstats} --cors cors --cv-proportion .1 --pseudos height --high-LD {highld} --window-cm 1 --extract {snplist}
    """



