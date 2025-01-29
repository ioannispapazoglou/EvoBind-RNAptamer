import subprocess
import shutil
import os

def hdocking(protein, workingdir, internalid, spacingstep, anglestep):

    # copy execution files
    shutil.copy(protein, workingdir)
    shutil.copy('./3pbins/hdock', workingdir)
    shutil.copy('./3pbins/createpl', workingdir)

    hdock = ['./hdock', 
            protein, 
            'minimized_aptamer.pdb', 
            '-out', f'Hdock{internalid}.out',
            '-spacing', f'{spacingstep}',
            '-angle', f'{anglestep}'
            ]


    createpl = ['./createpl', 
                f'Hdock{internalid}.out',
                'complex.pdb', 
                '-nmax', '1', 
                '-complex', 
                '-models'
                ]

    subprocess.run(hdock, check=True, cwd=workingdir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    #print('-------------------------------------------')
    subprocess.run(createpl, check=True, cwd=workingdir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    #print("===========================================")   

    # rename results
    shutil.move(workingdir+'model_1.pdb', workingdir+'complex.pdb')

    # remove execution files
    os.remove(os.path.join(workingdir, os.path.basename(protein)))
    os.remove(os.path.join(workingdir, 'hdock'))
    os.remove(os.path.join(workingdir, 'createpl'))
    os.remove(os.path.join(workingdir, f'Hdock{internalid}.out'))


#./hdock 7F6G_prep_nolig.pdb ligand_R08-0060942_S_000009.pdb 
#           -out Hdock.out 
#           -rsite "txt file with binding site residues"

#./createpl Hdock.out top10.pdb 
#           -nmax 10 
#           -complex 
#           -models

def getdockingscore(seq, cmplexpdbpath): 
    
    score = None

    with open(cmplexpdbpath, 'r') as file:
        for line in file:
            if line.startswith("REMARK Score:"):
                score = float(line.split(":")[1].strip())
                break

    normscore = score / len(seq) 

    return score, normscore
