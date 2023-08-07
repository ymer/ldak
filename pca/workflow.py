### SETTINGS
ldak = ".././ldak5.2.linux"
genome = "../../study_genome/genome"
d = "files"


### CODE
from os import makedirs
from gwf import Workflow
gwf = Workflow(defaults = {"account": "perjektet"})

def bim(s): 
    return f'{s}.bim'

def newdir(s):
    makedirs(s, exist_ok = True)
    return(s)

d_thin = newdir(f'{d}/thin')
weights_thin = f'{d_thin}/weights.thin'
gwf.target(
    f'create_weights',
    inputs = [bim(genome)],
    outputs = {'weights': weights_thin}
    ) << f"""
    {ldak} --thin {d_thin}/thin --bfile {genome} --window-prune .98 --window-kb 100
    awk < {d_thin}/thin.in '{{print $1, 1}}' > {weights_thin}
    """

d_kinship = newdir(f'{d}/kinship')
s_kinship = f'{d_kinship}/LDAK-Thin'
kinship_matrix = [f'{s_kinship}.{ext}' for ext in 
        ["grm.bin", "grm.id", "grm.details", "grm.adjust"]]
kinship = gwf.target(
    f'compute_kinship',
    inputs = [weights_thin, bim(genome)],
    outputs = kinship_matrix
    ) << f"""
    {ldak} --calc-kins-direct {s_kinship} --bfile {genome} --weights {weights_thin} --power -.25
    """

pca_vector = f'{s_kinship}.vect'
pca_values = f'{s_kinship}.values'
pca = gwf.target(
    f'compute_pca',
    inputs = kinship_matrix,
    outputs = [pca_vector, pca_values]
        ) << f"""
        {ldak} --pca {s_kinship} --grm {s_kinship} --axes 20
    """

pca_loadings = f'{s_kinship}.load'
loadings = gwf.target(
    f'obtain_loadings',
    inputs = [pca_vector, pca_values],
    outputs = [pca_loadings]
    ) << f"""
        {ldak} --calc-pca-loads {s_kinship} --pcastem {s_kinship} --grm {s_kinship} --bfile {genome}
    """

