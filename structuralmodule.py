from rhofold import Rhofolder
from hdock import hdocking, getdockingscore
from strscorers import RASPscore, get3Dstructurescore, Rgscore, isFolded, getMMGBSAbindingscore
from minimization import pdbeditor, proteinprep, aptamerprep, openmmrelax, delwaters, removehatms
from Bio.PDB import PDBIO, Atom, Residue, Chain, Model, Structure
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO

import os
from tqdm import tqdm
import pandas as pd
import random
import subprocess
import pickle

def mkdir(folderpath):
    if not os.path.exists(folderpath):
        os.makedirs(folderpath)

def writefasta(id, sequence, outputfile):
    with open(outputfile, 'w') as f:
        f.write(f"> Individual {id}\n{sequence}\n")

def readfasta(folderpath):
    with open(folderpath+'aptamer.fasta', "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            sequence = str(record.seq)

    return sequence

def seperatepdb(filesfolder, file):
    with open(filesfolder+file, 'r') as file:
        protein_lines = []
        aptamer_lines = []

        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain_id = line[21]
                if chain_id == 'A':
                    protein_lines.append(line)
                elif chain_id == 'B':
                    aptamer_lines.append(line)

    with open(filesfolder+'i_protein.pdb', 'w') as protein_file:
        protein_file.writelines(protein_lines)

    with open(filesfolder+'i_aptamer.pdb', 'w') as aptamer_file:
        aptamer_file.writelines(aptamer_lines)

def combinepdb(filesfolder):
    
    os.remove(filesfolder + "complex.pdb")
    complex_lines = []

    with open(filesfolder+'protein.pdb', 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                complex_lines.append(line)

    with open(filesfolder+'aptamer.pdb', 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                complex_lines.append(line)
    
    with open(filesfolder+'complex.pdb', 'w') as aptamer_file:
        aptamer_file.writelines(complex_lines)

def pdbconv(filesfolder, mode):

    atoms_data = []

    if mode == 'A':
        with open(filesfolder + 'i_aptamer.pdb', 'r') as pdb_file:
            pdb_content = pdb_file.readlines()
    
        for line in pdb_content:
            cnts = line.strip().split()
            atom = (
                int(cnts[1]),        # atom ID
                cnts[2],             # atom name
                cnts[3][1:],         # residue name
                cnts[4],             # chain ID
                int(cnts[5]),        # residue ID
                float(cnts[6]),      # x coordinate
                float(cnts[7]),      # y coordinate
                float(cnts[8]),      # z coordinate
                float(0.00),         # occupancy
                float(1.00),         # temperature factor (B-factor)
                cnts[2][0]           # element
            )
            atoms_data.append(atom)

    elif mode == 'P':
        with open(filesfolder + 'i_protein.pdb', 'r') as pdb_file:
            pdb_content = pdb_file.readlines()

        for line in pdb_content:
            cnts = line.strip().split()
            atom = (
                int(cnts[1]),        # atom ID
                cnts[2],             # atom name
                cnts[3],             # residue name
                cnts[4],             # chain ID
                int(cnts[5]),        # residue ID
                float(cnts[6]),      # x coordinate
                float(cnts[7]),      # y coordinate
                float(cnts[8]),      # z coordinate
                float(0.00),         # occupancy
                float(1.00),         # temperature factor (B-factor)
                cnts[2][0]           # element
            )
            atoms_data.append(atom)

    structure = Structure.Structure("aptamer")
    model = Model.Model(0)

    if mode == 'A':
        chain = Chain.Chain("B")
    
    elif mode == 'P':
        chain = Chain.Chain("A")

    current_residue = None

    for atom_data in atoms_data:
        atom_id, atom_name, residue_name, chain_id, residue_id, x, y, z, occupancy, b_factor, element = atom_data

        if current_residue is None or current_residue.get_id()[1] != residue_id:
            current_residue = Residue.Residue((' ', residue_id, ' '), residue_name, ' ')
            chain.add(current_residue)

        atom = Atom.Atom(
            atom_name,  # Atom name
            [x, y, z],  # Coordinates as a list
            occupancy,  # Occupancy
            b_factor,   # B-factor
            ' ',        # Alternate location indicator
            atom_name,  # Full atom name (for display in PDB)
            atom_id,    # Serial number (atom ID)
            element     # Element symbol
        )
        current_residue.add(atom)

    model.add(chain)
    structure.add(model)

    io = PDBIO()
    io.set_structure(structure)

    if mode == 'A':
        io.save(filesfolder + "aptamer.pdb")
        os.remove(filesfolder + "i_aptamer.pdb")
    
    elif mode == 'P':
        io.save(filesfolder + "protein.pdb")
        os.remove(filesfolder + "i_protein.pdb")

def characterizeproteinchain(input, output):

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", input)

    for model in structure:
        for chain in model:
            chain.id = 'P'

    io = PDBIO()
    io.set_structure(structure)
    io.save(output)


### ======================================================================== ###
### ============================== BIRTHPLACE ============================== ###
### ======================================================================== ###
"""                                    1                                     """
"""               Generating a random set of aptamer sequences               """
 
def sequencepoolbirth(samplesize, minlength, maxlength, peoplebefore):

    all_sequences = set()  # set to store unique sequences

    while len(all_sequences) < samplesize:
        length = random.randint(minlength, maxlength)
        sequence = ''.join(random.choices('AUGC', k=length))
        all_sequences.add(sequence)  # automatically handles duplicates

    #all_sequences = list(all_sequences)
    all_sequences = {peoplebefore + i + 1: seq for i, seq in enumerate(all_sequences)}

    return all_sequences
 
def mkpopulationdf(selectedseqdf, peoplebefore):

    df = selectedseqdf
    
    df['id'] = [f'{i+1}' for i in range(peoplebefore, peoplebefore+len(df))]

    return df

def seqprepare(sequences, peoplebefore):
    
    sequences = {peoplebefore + i: seq for i, seq in enumerate(sequences)}
    
    return sequences

### ======================================================================== ###
### ================================= END ================================== ###
### ======================================================================== ###


### ======================================================================== ###
### =========================== STRUCTURE MODULE =========================== ### 
### ======================================================================== ###
"""                                    2                                     """
"""Performing RhoFold and  seperate folded/ unfolded structures with Rg score"""

def evaluateaptamerstructures(args):
    aptamerseq, idfolder, nu = args
    aptamerid = int(idfolder.split('/')[-2])

    #seperatepdb(idfolder,file='complex.pdb')
    #pdbconv(idfolder, mode = 'A')
    #pdbconv(idfolder, mode = 'P')
    #combinepdb(idfolder)

    try:
        score_Rg, normscore_Rg = Rgscore(idfolder, nt_num = len(aptamerseq), nu = nu)
        #normscore_Rg = random.uniform(-10, 10)
    except:
        score_Rg = float('+inf')
        normscore_Rg = float('+inf')
    
    with open(f'{idfolder}Rg_score.txt', 'w') as file:
        file.write(str(normscore_Rg))

    return {'Rg_norm': normscore_Rg, 'id': aptamerid}


def structurespredictor(individuals, memory, generation, nu, allgenerationsfolder, numworkers):

    ######## -------------- DATA PREPARATION -------------- ####
    genfolder = allgenerationsfolder + f'gen{generation}'
    mkdir(genfolder)
    searchpool = pd.DataFrame(list(individuals.items()), columns = ['ID','Sequence'])
    searchpool.to_csv(f'{genfolder}/searchpool.csv', index=False) # this is rhofold's input

    for id, sequence in individuals.items(): 
        idfolder = genfolder + f'/{id}' 
        mkdir(idfolder)
        writefasta(id, sequence, idfolder+'/aptamer.fasta')
    ######## -------------- DATA PREPARATION -------------- ####

    ######## ----------- STRUCTURE MODULE (GPU) ----------- ####
    Rhofolder(inputdf=f'{genfolder}/searchpool.csv', outputfolder=genfolder)
    ######## ----------- STRUCTURE MODULE (GPU) ----------- ####

    ######## ---------- Rg STRUCTURE SCORER (CPU) --------- ####  
    with mp.Pool(numworkers) as pool:
        args = [(seq, genfolder + f'/{id}/', nu) for id, seq in individuals.items()]
        results = list(pool.imap_unordered(evaluateaptamerstructures, args))

    validated = pd.DataFrame(results, columns=['id', 'Rg_norm'])
    population = pd.DataFrame(columns=['id', 'seq', 'Rg_norm', 'RASP_norm', 'hdock_norm', 'MMGBSA_norm'])
    population = pd.concat([population, validated], ignore_index=True)

    population = population.sort_values(by='Rg_norm', ascending=True)#.head(selection)
    population.reset_index(drop=True, inplace=True)
    
    for index, individual in population.iterrows():
        idx = int(individual['id'])
        individualfolder = allgenerationsfolder + f'gen{generation}/{idx}/'
        population.at[index, 'seq'] = readfasta(individualfolder)             # bring sequence to population df

    foldedidlist = validated['id'].tolist() # list of current gen individuals
    memory.append(validated['id'].tolist()) # append new seqs to previous memory

    os.remove(f'{genfolder}/searchpool.csv')
    ######## ---------- Rg STRUCTURE SCORER (CPU) --------- ####

    return population, memory

 
### ======================================================================== ###
### ================================= END ================================== ### 
### ======================================================================== ###


### ======================================================================== ###
### ============================== MINIMIZER =============================== ### 
### ======================================================================== ###
"""                                    3                                     """
"""                     Minimizing with OpenMM (ff14/OL3)                    """

import subprocess
from multiprocessing import Pool
from tqdm import tqdm

# minimizer for aptamer

def aptaminimizer(args):
    workingdir, ministeps = args

    folderin = workingdir+'/temp/'
    os.makedirs(folderin, exist_ok=True)
    shutil.copy(f"{workingdir}/aptamer.pdb", folderin+"aptamer_unprep.pdb")

    aptamerprep(folderin)
    removehatms(folderin, mode = 'A')

    try:
        openmmrelax(folderin, ministeps=ministeps, mode = 'A')
        delwaters(folderin, mode = 'A')

        shutil.copy(folderin+'minimized_prediction.pdb', f'{workingdir}/minimized_aptamer.pdb')
        shutil.rmtree(folderin)
        os.remove(workingdir + "/aptamer.pdb")

    except Exception as e:
        print(f'> > > Minimization @ {workingdir} failed. Error was {e}.\n Passing.. < < <')
        shutil.copy(f"{workingdir}/aptamer.pdb", f"{workingdir}/minimized_aptamer.pdb") # act like it happened so it doesn't crush the rest of the code (will be eliminated by energy criteria)
      
def comminimizer(args):
    workingdir, ministeps = args

    folderin = workingdir+'/temp/'
    os.makedirs(folderin, exist_ok=True)
    shutil.copy(f"{workingdir}/complex.pdb", folderin+"complex_prep.pdb")
    
    try:    # If the doesnt have H atoms, it will fail. So, we need to copy the original file and pass it to the next step
        removehatms(folderin, mode = 'C')
    except:
        shutil.copy(f"{workingdir}/complex.pdb", f"{folderin}complex_prep_noH.pdb")
    
    try:
        openmmrelax(folderin, ministeps=ministeps, mode = 'C')
        delwaters(folderin, mode = 'C')

        shutil.copy(folderin+'minimized_prediction.pdb', f'{workingdir}/minimized_complex.pdb')
        shutil.rmtree(folderin)
        os.remove(workingdir + "/complex.pdb")

    except Exception as e:
        print(f'> > > Minimization @ {workingdir} failed. Error was {e}.\n Passing.. < < <')
        shutil.copy(f"{workingdir}/complex.pdb", f"{workingdir}/minimized_complex.pdb") # act like it happened so it doesn't crush the rest of the code (will be eliminated by energy criteria)
        
def accminimizer(foldedseqs, mode, minimizationsteps, generation, allgenerationsfolder, numworkers): # multi-CPU acceleration

    tasks = [
        (f"{allgenerationsfolder}gen{generation}/{aptamerid}", str(minimizationsteps))
        for aptamerid in foldedseqs['id']
    ]

    if mode == 'A':
        with Pool(numworkers) as pool:
            list(tqdm(pool.imap(aptaminimizer, tasks), total=len(tasks), desc='Minimizing 3D structures'))

    elif mode == 'C':
        with Pool(numworkers) as pool:
            list(tqdm(pool.imap(comminimizer, tasks), total=len(tasks), desc='Minimizing 3D structures'))



### ======================================================================== ###
### ================================= END ================================== ###
### ======================================================================== ###


### ======================================================================== ###
### ==================== MINIMIZED STRUCTURE EVALUATION ==================== ### 
### ======================================================================== ###
"""                                    4                                     """
"""    Evaluating sructures: RASP for RNA - DRPS for protein-RNA complex     """

import multiprocessing as mp  
import shutil 

def evaluatewithRASP(args):
    aptameridx, genfolder, foldedthresshold = args
    aptamerstructurefolder = genfolder + f'{aptameridx}/'

    try:
        foldedflag = isFolded(aptamerstructurefolder, theta=foldedthresshold)
        RASPscore(aptamerstructurefolder)
        
        score_3D, normscore_3D = get3Dstructurescore(aptamerstructurefolder) 
        
        score_3D = score_3D * foldedflag
        normscore_3D = normscore_3D * foldedflag # if not folded, score will be 0 (big RASP value >> will be discarded)
    
    except:
        print(f'{aptameridx} : excepted from batch')
        score_3D = float('+inf')
        normscore_3D = float('+inf')
    
    return {'RASP_norm': normscore_3D, 'id': aptameridx}


def structurescorer(foldedseqs, foldingthresshold, generation, allgenerationsfolder, numworkers):

    genfolder = allgenerationsfolder + f'gen{generation}/'

    ######## --------- RASP STRUCTURE SCORER (CPU) -------- #### 
    with mp.Pool(numworkers) as pool:
        args = [(idx, genfolder, foldingthresshold) for idx in foldedseqs['id']]
        results = list(tqdm(pool.imap_unordered(evaluatewithRASP, args), total=len(foldedseqs), leave=True, desc='Evaluating aptamer structure'))

    validated = pd.DataFrame(results, columns=['id', 'RASP_norm'])
    foldedseqs = foldedseqs.merge(validated, on='id', how='left')
    foldedseqs = foldedseqs[['id', 'seq', 'Rg_norm', 'RASP_norm_y', 'hdock_norm', 'MMGBSA_norm']]
    foldedseqs = foldedseqs.rename(columns={'RASP_norm_y': 'RASP_norm'})
    population = foldedseqs.sort_values(by='RASP_norm', ascending=True)#.head(populationsize)
    population.reset_index(drop=True, inplace=True)
    ######## --------- RASP STRUCTURE SCORER (CPU) -------- ####
    
    return population


### ======================================================================== ###
### ================================= END ================================== ### 
### ======================================================================== ###


### ======================================================================== ###
### =============================== DOCKING ================================ ### 
### ======================================================================== ###
"""                                    5                                     """
"""                           Docking with HDock                             """

def docking(params):
    idx, proteintarget, predictionoutputfolder, population, spacingstep, anglestep = params
    
    sequence = population.loc[population['id'] == idx, 'seq'].values[0]
    individualfolder = predictionoutputfolder + f'{idx}/'
    
    try: 
        hdocking(protein=proteintarget, workingdir=individualfolder, internalid=idx, spacingstep=spacingstep, anglestep=anglestep)
        fitnessscore, normfitnessscore = getdockingscore(sequence, individualfolder + 'complex.pdb')
    except Exception as e:
            print(f'{idx} : excepted from batch because: {e}')
            fitnessscore = float('+inf')
            normfitnessscore = float('+inf')
    
    with open(f'{individualfolder}Hdock_score.txt', 'w') as file:
        file.write(str(normfitnessscore))

    return (idx, normfitnessscore)

def dockaptamers(population, generation, proteintarget, spacingstep, anglestep, allgenerationsfolder, nworkers):
    
    ######## ----------- DOCKING MODULE (MULTICORE) ----------- ####
    predictionoutputfolder = allgenerationsfolder + f'gen{generation}/'
    
    #population = population.head(capacity)
    individuals = population['id'].tolist()
    params = [(idx, proteintarget, predictionoutputfolder, population, spacingstep, anglestep) for idx in individuals]

    with mp.Pool(nworkers) as pool:
        results = list(tqdm(pool.imap(docking, params), total=len(params), desc=f'Performing molecular docking'))

    for id, fitnessscore in results:            
        population.loc[population['id'] == id, 'hdock_norm'] = fitnessscore
    ######## ----------- DOCKING MODULE (MULTICORE) ----------- ####
    
    return population


### ======================================================================== ###
### ================================= END ================================== ###
### ======================================================================== ###


### ======================================================================== ###
### ============================= MMGBSA MODULE ============================ ###
### ======================================================================== ###
"""                                    6                                     """
"""              Estimate binding affinity with Amber MMPBSA.py              """

def ambermmgbsa(params):
    idx, predictionoutputfolder = params
    
    individualfolder = predictionoutputfolder + f'{idx}/'
    mmgbsas = '/media/arma/DATA.ext4/giannis_phd/geneticalgo2/mmgbsas'
    workingdir = individualfolder + 'mmgbsas'
    shutil.copytree(mmgbsas, workingdir)

    shutil.copy(individualfolder+'minimized_complex.pdb',workingdir)

    try:
        cmd = f'python3 energyest.py'
        subprocess.run(cmd, shell=True, check=True, cwd=workingdir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        shutil.copy(workingdir+'/bindingenergy.dat', individualfolder+'MMGBSA.txt')
        shutil.copy(workingdir+'/apt_decomp.csv', individualfolder+'MMGBSA_decomp.csv')
        shutil.rmtree(workingdir)
        mmgbsascore = getMMGBSAbindingscore(individualfolder)
    except:
        print(f'> > > MMGBSA @ {idx} failed < < <')
        mmgbsascore = float('+inf')
        decomp = {'Residue': ['-   -'], 'TOTAL': [0.0]}
        pd.DataFrame(decomp).to_csv(individualfolder+'MMGBSA_decomp.csv', index=False)
        #pass
    
    return {'MMGBSA': mmgbsascore, 'id': idx} 

def bindingestimator(population, generation, allgenerationsfolder, numworkers):
    
    predictionoutputfolder = allgenerationsfolder + f'gen{generation}/'
    
    individuals = population['id'].tolist()
    params = [(aptameridx, predictionoutputfolder) for aptameridx in individuals]

    with mp.Pool(numworkers) as pool:
        results = list(tqdm(pool.imap(ambermmgbsa, params), total=len(params), desc=f'Calculating frame MM/GBSA'))

    for result in results:
        mmgbsascore = result['MMGBSA']
        id = result['id']
        seq_length = int(population.loc[population['id'] == id, 'seq'].str.len().values[0])
        
        population.loc[population['id'] == id, 'MMGBSA_norm'] = mmgbsascore / seq_length
    
    return population 


### ======================================================================== ###
### ================================= END ================================== ###
### ======================================================================== ###


### ======================================================================== ###
### =========================== MUTATION MODULE ============================ ###
### ======================================================================== ###
"""                                    7                                     """
"""                        Currate population data                           """


def get_decomp_scores(genfolder):
    
    results = [] 
    
    for root, dirs, files in os.walk(genfolder):
        if 'MMGBSA_decomp.csv' in files:
            csv_file = os.path.join(root, 'MMGBSA_decomp.csv')

            # extract idx
            parts = os.path.normpath(csv_file).split(os.sep)
            idx = int(next((part for part in parts if part.isdigit()), None))

            # extract decomp scores (tupples)
            df = pd.read_csv(csv_file)
            df['Residue'] = df['Residue'].str.strip()
            base_scores = tuple(df['TOTAL'].round(2))
            
            results.append({'id': idx, 'Decomp': base_scores})

    return pd.DataFrame(results)

def summarizepopulation(population, mode, populationsize, parentalmemory, generation, allgenerationsfolder):
    
    curfolder = allgenerationsfolder + f'gen{generation}/'

    # - - - - - currate current sequences - - - - - #
    selectedpopulation = population.sort_values(by='MMGBSA_norm', ascending=True).head(populationsize) # current selection with MMGBSA alone
    for index, individual in selectedpopulation.iterrows():
        idx = int(individual['id'])
    # - - - - - currate current sequences - - - - - #

    #- - - - - bring memory individuals - - - - - - #
    if generation > 0:
        selectedpopulation = pd.concat([selectedpopulation, parentalmemory], ignore_index=True)
        for _, individual in selectedpopulation.iterrows():
            idx = str(individual['id']).split('.')[0]
            try:
                alreadycalculated = allgenerationsfolder + f'gen{generation-1}/{idx}/'
                individualfolder = allgenerationsfolder + f'gen{generation}/{idx}/'
                shutil.copytree(alreadycalculated, individualfolder)
            except:
                pass  
    #- - - - - bring memory individuals - - - - - - #
    
    if mode == 'Normal':
        parentalmemory.to_csv(f'{curfolder}parentalmem.csv')
        selectedpopulation.sort_values(by='MMGBSA_norm').to_csv(f'{curfolder}summary.csv', index=False)
    elif mode == 'MinimaJump':
        selectedpopulation.sort_values(by='MMGBSA_norm').to_csv(f'{curfolder}enhancersummary.csv', index=False)
        
    if generation == 0:
        decomp = get_decomp_scores(curfolder)
        selectedpopulation['generation'] = generation
        combined = pd.merge(selectedpopulation, decomp, on='id', how='inner')
        combined.sort_values(by='MMGBSA_norm').to_csv(f'{allgenerationsfolder}history.csv', index=False)
    else:
        decomp = get_decomp_scores(curfolder)
        selectedpopulation['generation'] = generation
        combined = pd.merge(selectedpopulation, decomp, on='id', how='inner')
        
        history = pd.read_csv(f'{allgenerationsfolder}history.csv')
        
        # Append only entries not already in the history file (based on the 'sequence' column)
        new_entries = combined.sort_values(by='MMGBSA_norm')[
            ~combined.sort_values(by='MMGBSA_norm')['seq'].isin(history['seq'])
        ]

        updated_history = pd.concat([history, new_entries], ignore_index=True)
        updated_history.to_csv(f'{allgenerationsfolder}history.csv', index=False)

    return selectedpopulation

 
### ======================================================================== ###
### ================================= END ================================== ###
### ======================================================================== ###


### ======================================================================== ###
### ================================= END ================================== ###
### ======================================================================== ###
"""                                    ~                                     """
"""                                 ensemble                                 """

from adaptivemutation import exportmutates

def structuremodule(inpsequences, generation, peoplebefore, 
                    initpopulation, nu, ssthreshold, 
                    nworkers, minimizationsteps, hdockspacingstep, 
                    hdockanglestep, proteintarget, allgenerationsfolder):
    
    # 0. Fix nsga genrated input to be processed by the structure module
    birthpool = seqprepare(inpsequences, peoplebefore)
    
    # 1. Structure Module (RhoFold - Rg selection)
    foldedaptamersdf, memory = structurespredictor(individuals=birthpool, 
                                                    memory = [], # empty list (no memory)
                                                    generation = generation,
                                                    nu = nu, 
                                                    allgenerationsfolder = allgenerationsfolder,
                                                    numworkers = nworkers) 

    # 2. Fisrt Structure Refinement Module (OpenMM)        
    accminimizer(foldedseqs = foldedaptamersdf, 
                    mode = 'A',
                    minimizationsteps = minimizationsteps,
                    generation = generation, 
                    allgenerationsfolder = allgenerationsfolder, 
                    numworkers = nworkers)  

    # 3. Initial Structure Evaluation Module (RASP)
    foldedaptamersdf = structurescorer(foldedseqs = foldedaptamersdf, 
                                        foldingthresshold = ssthreshold,
                                        generation = generation, 
                                        allgenerationsfolder = allgenerationsfolder,
                                        numworkers = nworkers)

    # 4. Docking module (HDock)
    foldedaptamersdf = dockaptamers(population = foldedaptamersdf, 
                                    generation = generation, 
                                    proteintarget = proteintarget, 
                                    spacingstep=hdockspacingstep, 
                                    anglestep=hdockanglestep,
                                    allgenerationsfolder = allgenerationsfolder, 
                                    nworkers = nworkers)

    # 5. Second Structure Refinement Module (OpenMM)        
    accminimizer(foldedseqs = foldedaptamersdf, 
                    mode = 'C',
                    minimizationsteps = minimizationsteps,
                    generation = generation, 
                    allgenerationsfolder = allgenerationsfolder, 
                    numworkers = nworkers)  

    # 6. Bindding Energy Estimation Module (MMGBSA.py)
    foldedaptamersdf = bindingestimator(population = foldedaptamersdf, 
                                        generation = generation, 
                                        allgenerationsfolder = allgenerationsfolder, 
                                        numworkers = nworkers)

    # 7. Sumarizing Population
    population = summarizepopulation(population = foldedaptamersdf, 
                                        mode = 'Normal',
                                        populationsize = initpopulation,
                                        parentalmemory = pd.DataFrame(columns=['id', 'seq', 'Rg_norm', 'RASP_norm', 'hdock_norm', 'MMGBSA_norm']), # empty df (no parental memory)
                                        generation = generation, 
                                        allgenerationsfolder = allgenerationsfolder)

    mutprobs = exportmutates(history = f'{allgenerationsfolder}history.csv', population = population, savedir=allgenerationsfolder, generation=generation)
    
    # Update peoplebefore
    peoplebefore = peoplebefore + initpopulation
    
    return population, peoplebefore, mutprobs


### ======================================================================== ###
### ================================= END ================================== ###
### ======================================================================== ###
