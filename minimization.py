from simtk.unit import *
from openmm.app import *
from openmm import *
import openmm as mm
from sys import stdout
from Bio.PDB import PDBParser, PDBIO, Atom



def pdbeditor(folderin, mode):
    # ~~~~~~~~~~~ #
    # Combination / Seperation of system components
    # ~~~~~~~~~~~ #
    if mode == 'SEP':
        with open(folderin+'complex_unprep.pdb', 'r') as file:
            protein_lines = []
            aptamer_lines = []

            for line in file:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    chain_id = line[21]
                    if chain_id == 'P':
                        protein_lines.append(line)
                    elif chain_id == 'A':
                        aptamer_lines.append(line)

        with open(folderin+'protein.pdb', 'w') as protein_file:
            protein_file.writelines(protein_lines)

        with open(folderin+'aptamer.pdb', 'w') as aptamer_file:
            aptamer_file.writelines(aptamer_lines)
    
    elif mode == 'COM':
        complex_lines = []

        with open(folderin+'protein_prep.pdb', 'r') as file:
            for line in file:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    complex_lines.append(line)

        with open(folderin+'aptamer_prep.pdb', 'r') as file:
            for line in file:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    complex_lines.append(line)

        with open(folderin+'complex_prep.pdb', 'w') as aptamer_file:
            aptamer_file.writelines(complex_lines)
        


def aptamerprep(folderin):
    # ~~~~~~~~~~~ #
    # Deletion of initial PO2 group 
    # ~~~~~~~~~~~ #

    parser = PDBParser()
    structure = parser.get_structure("aptamer", folderin+"aptamer_unprep.pdb")

    model = structure[0]
    chain = model['0'] 

    atom_names_to_remove = {"OP1", "OP2", "P"}

    first_residue = next(chain.get_residues())

    atoms_to_remove = [atom for atom in first_residue if atom.get_name() in atom_names_to_remove]

    for atom in atoms_to_remove:
        first_residue.detach_child(atom.get_id())

    io = PDBIO()
    io.set_structure(structure)
    io.save(folderin+"aptamer_prep.pdb")

    """
    serial_numbers_to_remove = {1, 2, 3}  # O - P - O of first residue

    for residue in chain:
        atoms_to_remove = [atom for atom in residue if atom.get_serial_number() in serial_numbers_to_remove]
        for atom in atoms_to_remove:
            residue.detach_child(atom.get_id())
    """

def proteinprep(folderin):
    # ~~~~~~~~~~~ #
    # Addition of terminal OXT
    # ~~~~~~~~~~~ #

    parser = PDBParser()
    structure = parser.get_structure("protein", folderin+"protein_unprep.pdb")

    model = structure[0]
    chain = model['A'] 

    res_num = 0
    atm_num = 0

    for model in structure:
        for chain in model:
            for residue in chain:
                res_num += 1
                atm_num += len(residue)

    for residue in chain:
        residx = residue.get_id()[1]
        for atom in residue:
            if (atom.get_name() == 'O') and (residx == res_num):
                    serial_number=atm_num+1
                    coordinates=atom.get_coord()
                    bfactor=atom.get_bfactor()
                    occupancy=atom.get_occupancy()
                    altloc=atom.get_altloc()
                    element=atom.element

    newatm = Atom.Atom( name="OXT",
                        coord=[round(coord + 0.1, 3) for coord in coordinates],
                        bfactor=bfactor,
                        occupancy=occupancy,
                        altloc=altloc,
                        fullname=" OXT",
                        serial_number=atm_num+1,
                        element=element )
    residue.add(newatm)

    io = PDBIO()
    io.set_structure(structure)
    io.save(folderin+"protein_prep.pdb")



def removehatms(folderin, mode):
    # ~~~~~~~~~~~ #
    # Deletion of initial hydrogens
    # ~~~~~~~~~~~ #

    parser = PDBParser(QUIET=True)
    if mode == 'C':
        structure = parser.get_structure("system", folderin+'complex_prep.pdb')
    elif mode == 'A':
        structure = parser.get_structure("system", folderin+'aptamer_prep.pdb')

    for model in structure:
        for chain in model:
            for residue in chain:
                atoms_to_remove = [atom for atom in residue if atom.element == 'H']
                for atom in atoms_to_remove:
                    residue.detach_child(atom.get_id())

    io = PDBIO()
    io.set_structure(structure)
    if mode == 'C':
        io.save(folderin+'complex_prep_noH.pdb')
    elif mode == 'A':
        io.save(folderin+'aptamer_prep_noH.pdb')



def calculatedpredictedenergy(simulationobj):

    initialenergy = simulationobj.context.getState(getEnergy=True).getPotentialEnergy()
    
    return initialenergy._value * KcalPerKJ



def openmmrelax(folderin, ministeps, mode):
    
    # ~~~~~~~~~~~ #
    # Designed by following openmm guide @ http://docs.openmm.org/latest/userguide/application/03_model_building_editing.html
    # For detailed explanation of each step, read here ^ ^ ^ ^ ^
    # ~~~~~~~~~~~ #
    
    #print('Minimizing.. ')
    if mode == 'C':
        pdb = PDBFile(folderin+'complex_prep_noH.pdb')
    elif mode == 'A':
        pdb = PDBFile(folderin+'aptamer_prep_noH.pdb')

    modeller = Modeller(pdb.topology, pdb.positions)

    # system build (protonation - simulation box)
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    modeller.addHydrogens(forcefield, pH=7.4)
    modeller.addSolvent(forcefield, padding=1 * nanometer) # tip3p, NaCl is used here by default
  

    system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, nonbondedCutoff=1 * nanometer, constraints=HBonds)
    integrator = LangevinIntegrator(310 * kelvin, 1 / picosecond, 0.002 * picoseconds)

    platform = mm.Platform.getPlatformByName('CUDA')
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)
    simulation.reporters.append(StateDataReporter(stdout, ministeps, step=True, potentialEnergy=True, temperature=True))
    simulation.minimizeEnergy(maxIterations=int(ministeps))

    position = simulation.context.getState(getPositions=True).getPositions()
    
    if mode == 'C':
        app.PDBFile.writeFile(simulation.topology, position, open(folderin+'minimized_system.pdb', 'w'))
    elif mode == 'A':
        app.PDBFile.writeFile(simulation.topology, position, open(folderin+'minimized_aptamer.pdb', 'w'))
        
    #energy = simulation.context.getState(getEnergy=True).getP
    #print(f'Energy at Minima is', energy._value * KcalPerKJ, 'kcal/mol')



def delwaters(folderin, mode):
    # ~~~~~~~~~~~ #
    # Deletion of waters-ions from minimized structure
    # ~~~~~~~~~~~ #
    if mode == 'C':
        with open(folderin+'minimized_system.pdb', 'r') as infile, open(folderin+'minimized_prediction.pdb', 'w') as outfile:
            for line in infile:
                if line.startswith("ATOM"):
                    outfile.write(line)
    elif mode == 'A':
        with open(folderin+'minimized_aptamer.pdb', 'r') as infile, open(folderin+'minimized_prediction.pdb', 'w') as outfile:
            for line in infile:
                if line.startswith("ATOM"):
                    outfile.write(line)


# EXCECUTION :
"""
import os
import shutil

folderin = './temp/'
os.makedirs(folderin, exist_ok=True)
shutil.copy("./complex.pdb", folderin+"complex_unprep.pdb")

pdbeditor(folderin, mode='SEP')
proteinprep(folderin)
aptamerprep(folderin)
pdbeditor(folderin, mode='COM')
removehatms(folderin)
openmmrelax(folderin)
delwaters(folderin)

shutil.copy(folderin+'minimized_prediction.pdb', './minimized_prediction.pdb')
shutil.rmtree(folderin)
"""
# - - - - - - 