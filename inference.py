import time
import argparse
import numpy as np

from structuralmodule import structuremodule
from adaptivemutation import AdaptiveMutation

from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.problem import Problem
from pymoo.optimize import minimize

#################################################################################################
#################################################################################################

BASES = ["A", "C", "G", "U"]

class Evolution(Problem):
    def __init__(self, minlength, maxlength, population, ssthreshold, nworkers, minimizationsteps, hdockspacingstep, hdockanglestep, proteintarget, allgenerationsfolder):
        super().__init__(
            n_var=1 + maxlength,                            # 1 for length + maxlength for bases
            n_obj=2,                                        # Two objectives
            n_constr=0,                                     # No constraints
            xl=[minlength] + [0] * maxlength,               # Length min + Base min
            xu=[maxlength] + [len(BASES) - 1] * maxlength   # Length max + Base max
        )
        self.min_seq_length = minlength
        self.max_seq_length = maxlength
        self.gener = 0
        self.peoplebefore = 0
        self.population = population
        self.ssthreshold = ssthreshold
        self.nworkers = nworkers
        self.minimizationsteps = minimizationsteps
        self.hdockspacingstep = hdockspacingstep
        self.hdockanglestep = hdockanglestep
        self.proteintarget = proteintarget
        self.allgenerationsfolder = allgenerationsfolder

    def _evaluate(self, X, out, *args, **kwargs):
        population_size = X.shape[0]
        
        # Extract lengths and sequences
        lengths = X[:, 0].astype(int)
        raw_sequences = X[:, 1:].astype(int)

        decoded_sequences = [
            ''.join([BASES[int(base)] for base in raw_sequences[i, :lengths[i]]])
            for i in range(population_size)
        ]

        # Remove empty or invalid sequences
        decoded_sequences = [seq for seq in decoded_sequences if seq]

        # Structural module call
        populationdf, peoplebefore, mutprobs = structuremodule(
                                                inpsequences=decoded_sequences,
                                                generation=self.gener,
                                                peoplebefore=self.peoplebefore,
                                                initpopulation=self.population,
                                                nu=0.41,
                                                ssthreshold=self.ssthreshold,
                                                nworkers=self.nworkers,
                                                minimizationsteps=self.minimizationsteps,
                                                hdockspacingstep=self.hdockspacingstep,
                                                hdockanglestep=self.hdockanglestep,
                                                proteintarget=self.proteintarget,
                                                allgenerationsfolder=self.allgenerationsfolder
        )

        # Objective evaluation
        f1_values = np.array([
            populationdf.loc[populationdf['seq'] == seq, 'RASP_norm'].iloc[0]
            for seq in decoded_sequences
        ])

        f2_values = np.array([
            populationdf.loc[populationdf['seq'] == seq, 'MMGBSA_norm'].iloc[0]
            for seq in decoded_sequences
        ])

        # Assign the objectives
        out["F"] = np.column_stack((f1_values, f2_values))

        self.gener += 1
        self.peoplebefore = 0#peoplebefore # if update -- messes up with the verbose


#################################################################################################
#################################################################################################

start = time.time()

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Search for RNA aptamers againist a protein structure. Might need to change parameters:")

    # arguments
    parser.add_argument("--kmer", type=int, default=4, help="K-mer size used in adaptive mutation (default: 4).")
    parser.add_argument("--minlength", type=int, default=10, help="Minimum sequence length (default: 10).")
    parser.add_argument("--maxlength", type=int, default=30, help="Maximum sequence length (default: 30).")
    parser.add_argument("--population", type=int, default=10, help="Population size (default: 10).")
    parser.add_argument("--generations", type=int, default=10, help="Number of generations (default: 10).")
    parser.add_argument("--ssthreshold", type=float, default=0.7, help="Secondary structure threshold (default: 0.7).")
    parser.add_argument("--nworkers", type=int, default=30, help="Number of CPU's (default: 30).")
    parser.add_argument("--minimizationsteps", type=int, default=500, help="Number of minimization steps (default: 500).")
    parser.add_argument("--hdockspacingstep", type=float, default=2.0, help="HDock spacing step size (default: 2.0).")
    parser.add_argument("--hdockanglestep", type=float, default=30.0, help="HDock angle step size (default: 30.0).")
    parser.add_argument("--proteintarget", type=str, required=True, help="Path to the protein-target (.pdb) file.")
    parser.add_argument("--outputfolder", type=str, default="./nsga2autopilot/", help="Path to the output folder (default: ./nsga2autopilot/).")

    args = parser.parse_args()

    problem = Evolution(args.minlength, 
                        args.maxlength, 
                        args.population, 
                        args.ssthreshold,
                        args.nworkers, 
                        args.minimizationsteps, 
                        args.hdockspacingstep,
                        args.hdockanglestep, 
                        args.proteintarget, 
                        args.outputfolder)

    algorithm = NSGA2(pop_size=args.population,
                    mutation=AdaptiveMutation(directory=args.outputfolder, 
                                                k=args.kmer, 
                                                minlength=args.minlength, 
                                                maxlength=args.maxlength)
    )

    res = minimize(
        problem,
        algorithm,
        termination=('n_gen', args.generations),
        seed=1,
        save_history=True,
        verbose=True
    )

if __name__ == "__main__":
    main()
    
end = time.time()

print(f"Time passed: {(end - start)/3600:.2f} hours")
print('Done! Bye.')

#################################################################################################
#################################################################################################
