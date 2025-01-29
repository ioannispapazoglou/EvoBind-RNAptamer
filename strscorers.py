from Bio.PDB import PDBParser
import numpy as np
import subprocess
import shutil
import os


#=========================================================#
#                         RASP                            #
#=========================================================#


def RASPscore(workingdir):

    shutil.copy('./3pbins/rasp_fd', workingdir)

    command = f"./rasp_fd -e all -p minimized_aptamer.pdb > RASP.txt"

    subprocess.run(command, shell=True, check=True, cwd=workingdir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    os.remove(os.path.join(workingdir, 'rasp_fd'))

#./rasp_fd -e all -p MOLECULE.pdb > outputname.txt 	

def get3Dstructurescore(workingdir):

    score_3D = None

    file_3D = workingdir + 'RASP.txt'

    with open(file_3D, 'r') as file:
        
        for line in file:
            values = line.split()
            score_3D = float(values[0])
            normscore_3D = float(values[2])
            break

    return score_3D, normscore_3D


#=========================================================#
#                  SECONDARY STRUCTURE                    #
#=========================================================#


def parsect(file_path):
    pairs = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines[1:]:
            fields = line.split()
            resN = int(fields[0])
            resM = int(fields[4])
            pairs.append((resN, resM))
    return pairs

def pairs2bin(pairs):
    max_res = max(res for res, _ in pairs)
    matrix = [[0] * max_res for _ in range(max_res)]

    for resN, resM in pairs:
        if resM != 0:
            matrix[resN - 1][resM - 1] = 1
            matrix[resM - 1][resN - 1] = 1

    return matrix

def bin2db(matrix):
    n = len(matrix)
    dot_bracket = ["."] * n
    memory1 = []  # Stack for primary level nesting ()
    memory2 = []  # Stack for secondary level nesting []

    for i in range(n):
        for j in range(i + 1, n):
            if matrix[i][j] == 1:
                if not memory1 or (memory1 and memory1[-1] < i):
                    memory1.append(i)
                    dot_bracket[i] = "i"
                    dot_bracket[j] = "i"
                else:
                    memory2.append(i)
                    dot_bracket[i] = "["
                    dot_bracket[j] = "]"
    
    return "".join(dot_bracket)

def isFolded(predfolder, theta=0.7):

    ctfile = f'{predfolder}ss.ct'
    
    pairs = parsect(ctfile)
    matrix = pairs2bin(pairs)
    dotbracket = bin2db(matrix)

    if (list(dotbracket).count('i') / len(dotbracket)) <= theta: 
        return float(0)
    else:
        return float(1)


#=========================================================#
#                        MMGBSA                           #
#=========================================================#


def getMMGBSAbindingscore(workingdir):
    delta = None
    with open(workingdir+'MMGBSA.txt', 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith("DELTA TOTAL"):
                delta = float(line.split()[2])
                break

    return delta


#=========================================================#
#                       GYRATION                          #
#=========================================================#


atomic_masses = {
    'C': 12.01,  # Carbon
    'N': 14.01,  # Nitrogen
    'O': 16.00,  # Oxygen
    'P': 30.97,  # Phosphorus
    'H': 1.01    # Hydrogen (though excluded)
}

def compute_center_of_mass(atom_list):
    """Compute the center of mass for a list of atoms """
    total_mass = 0
    center_of_mass = np.zeros(3)
    for atom in atom_list:
        mass = atom.mass if hasattr(atom, 'mass') else 1.0  # Simplification: assume all atoms have a mass of 1.0
        center_of_mass += mass * atom.coord
        total_mass += mass
    center_of_mass /= total_mass
    return center_of_mass

def compute_radius_of_gyration(pdb_file):
    """Compute the radius of gyration for an RNA structure """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('RNA', pdb_file)

    atom_list = [atom for atom in structure.get_atoms() if atom.element != 'H']

    center_of_mass = compute_center_of_mass(atom_list)

    rg_square_sum = 0.0
    for atom in atom_list:
        distance = np.linalg.norm(atom.coord - center_of_mass)
        rg_square_sum += distance ** 2

    radius_of_gyration = np.sqrt(rg_square_sum / len(atom_list))
    return radius_of_gyration


def Rgscore(workingdir, nt_num, nu):

    pdb_file = workingdir + '/aptamer.pdb'
    rg_score = compute_radius_of_gyration(pdb_file)

    nu = 0.5
    rg_score_norm = rg_score / (nt_num ** nu)

    return rg_score, rg_score_norm


def Rgscore_imm(pdb_file):

    rg = compute_radius_of_gyration(pdb_file)

    return rg

def Rgscore_imm_norm(pdb_file, nt_num, nu):

    rg_score = compute_radius_of_gyration(pdb_file)

    rg_score_norm = rg_score / (nt_num ** nu)

    return rg_score, rg_score_norm