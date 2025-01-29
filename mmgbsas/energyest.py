from Bio.PDB import PDBParser, PDBIO
import pandas as pd
import re

def pdbeditor(mode):
    # ~~~~~~~~~~~ #
    # Combination / Seperation of system components
    # ~~~~~~~~~~~ #
    if mode == 'SEP':
        with open('minimized_complex.pdb', 'r') as file:
            protein_lines = []
            aptamer_lines = []

            for line in file:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    chain_id = line[21]
                    if chain_id == 'A':
                        protein_lines.append(line)
                    elif chain_id == 'B':
                        aptamer_lines.append(line)

        with open('protein_unprep.pdb', 'w') as protein_file:
            protein_file.writelines(protein_lines)

        with open('aptamer_prep.pdb', 'w') as aptamer_file:
            aptamer_file.writelines(aptamer_lines)
    
    elif mode == 'COM':
        complex_lines = []

        with open('protein_prep.pdb', 'r') as file:
            for line in file:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    complex_lines.append(line)

        with open('aptamer_prep.pdb', 'r') as file:
            for line in file:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    complex_lines.append(line)

        with open('complex_prep.pdb', 'w') as aptamer_file:
            aptamer_file.writelines(complex_lines)
        

def proteinprep(folderin):
    # ~~~~~~~~~~~ #
    # Labeling issues
    # ~~~~~~~~~~~ #

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", folderin+"protein_unprep.pdb")

    model = structure[0]
    chain = model['A']

    res_num = 0

    for model in structure:
        for chain in model:
            for residue in chain:
                res_num += 1

                # delete initial H (messes up force field)
                if res_num == 1:
                    atoms_to_delete = []
                    for atom in residue:
                        if atom.serial_number == 2:
                            atoms_to_delete.append(atom)
                    
                    for atom in atoms_to_delete:
                        residue.detach_child(atom.id)

                # specify HIS protonation states
                if residue.get_resname() == "HIS":
                    has_hd1 = any(atom.get_name() == "HD1" for atom in residue)
                    residue_id = residue.get_id()

                    if has_hd1:
                        # Rename HIS to HID if it has an HD1 atom
                        residue.resname = "HID"
                    else:
                        # Rename HIS to HIE if it doesn't have an HD1 atom
                        residue.resname = "HIE"

    io = PDBIO()
    io.set_structure(structure)
    io.save(folderin + "protein_prep.pdb")


def extractrnascores(file_path):

    with open(file_path, "r") as file:
        data = file.readlines()

    start_line = None
    for i, line in enumerate(data):
        if "Total Energy Decomposition:" in line:
            start_line = i
            break

    decomposition_data = data[start_line + 1:]

    rna_data = []
    for line in decomposition_data:
        if line.strip() == "":
            break
        if "Backbone" in line or "Sidechain" in line:
            continue
        if re.match(r"\s*(?:[AUCG]|[AUCG][35])\s+\d+,L", line):
            rna_data.append(line.strip())

    columns = [
        "Residue", "Location", "Internal", "Internal_Std_Dev", "Internal_Std_Err",
        "van_der_Waals", "VDW_Std_Dev", "VDW_Std_Err", "Electrostatic", "Electrostatic_Std_Dev",
        "Electrostatic_Std_Err", "Polar_Solvation", "Polar_Solv_Std_Dev", "Polar_Solv_Std_Err",
        "Non_Polar_Solvation", "NonPolar_Solv_Std_Dev", "NonPolar_Solv_Std_Err", "TOTAL", "Total_Std_Dev", "Total_Std_Err"
    ]

    rna_df = pd.DataFrame([line.split(",") for line in rna_data], columns=columns)
    selected_columns = ["Residue", "TOTAL"]
    rna_df = rna_df[selected_columns]

    rna_df.to_csv("apt_decomp.csv", index=False)

# EXCECUTION :

import os
import shutil

folderin = './'

pdbeditor(mode='SEP')
proteinprep(folderin)

#os.remove(folderin+'minimized_complex.pdb')
os.remove('./protein_unprep.pdb')

import subprocess

cmd = 'tleap -f leap.in'
subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

cmd = 'MMPBSA.py -O -i mmgbsa.in -cp complex.prmtop -rp protein.prmtop -lp rna.prmtop -y complex.inpcrd -o bindingenergy.dat'
subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

extractrnascores("FINAL_DECOMP_MMPBSA.dat")

#cmd = 'MMPBSA.py --clean'
#subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
