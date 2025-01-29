from scipy.special import softmax
from collections import Counter
import pandas as pd
import numpy as np
import pickle
import random



def kmerentropy(sequences, k=4):
    kmer_data = {}
    
    for seq in sequences:
        # Extract k-mers from the sequence
        kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]
        kmer_counts = Counter(kmers)
        total = sum(kmer_counts.values())
        freqs = {kmer: count / total for kmer, count in kmer_counts.items()} if total > 0 else {}

        # Calculate Shannon entropy for each k-mer
        for kmer, p in freqs.items():
            H = -p * np.log2(p) if p > 0 else 0
            if kmer in kmer_data:
                kmer_data[kmer] += H 
            else:
                kmer_data[kmer] = H
    
    return kmer_data



def kmerscoring(decomp, sequence, k, kmer_scores):
    
    base_scores = decomp.split('(')[1].split(')')[0].split(',')
    base_scores = [float(score.strip()) for score in base_scores]

    kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
    
    for i, kmer in enumerate(kmers):
        scores = []
        base_indices = range(i, i+k)  # Positions of the bases in the k-mer
        
        for j, base in zip(base_indices, kmer):
            scores.append(base_scores[j])
        
        average_score = round(sum(scores) / k, 2)
        
        if kmer not in kmer_scores:
            kmer_scores[kmer] = []
            
        kmer_scores[kmer].append(average_score) 
        
    return kmer_scores

import numpy as np



def boltzmann_weighting(kmer_scores, k_constant=1, T=1):
    
    representative_values = {}
    
    for kmer, values in kmer_scores.items():
        boltzmann_factors = np.exp(-np.array(values) / (k_constant * T))
        probabilities = boltzmann_factors / boltzmann_factors.sum()
        representative_value = np.sum(probabilities * values)
        representative_values[kmer] = representative_value
    
    return representative_values



def mutfunction(Hmer_scores, Bmer_scores):
    mutprobs = {}
    combined_values = {}

    for key in Hmer_scores.keys():
        if key in Bmer_scores:
            H_value = Hmer_scores[key]
            B_value = Bmer_scores[key]

            # Compute H * B term
            combined_values[key] = H_value * B_value

    # Normalize using softmax globally across all kmers
    keys = list(combined_values.keys())
    values = list(combined_values.values())
    softmax_values = softmax(values)

    # Assign normalized probabilities back to each key
    for i, key in enumerate(keys):
        mutprobs[key] = softmax_values[i]

    return mutprobs
    
    

def exportmutates(history, population, savedir, generation):
    
    data = pd.read_csv(history)#, nrows=1300)
    kmer_scores = {}
    
    evolutionterms = kmerentropy(data['seq'], k=4)

    for parentsequence in population['seq']:

        decomp = str(data.loc[data['seq'] == parentsequence, 'Decomp'].values)

        if decomp == "(0.0,)":  # In case MMGBSA decomp was skipped out of error, don't try to find any values
            continue
        
        else:
            kmer_scores = kmerscoring(decomp, sequence=parentsequence, k=4, kmer_scores=kmer_scores)  
            bindingterms = boltzmann_weighting(kmer_scores)
            mutprobs = mutfunction(evolutionterms, bindingterms)
    
    with open(f'{savedir}gen{generation}/mutprobs.pkl', 'wb') as file:
        pickle.dump(mutprobs, file)
    
    with open(f'{savedir}summutprobs.pkl', 'wb') as file:
        pickle.dump(mutprobs, file)

    return mutprobs


#########################################################################################################
#########################################################################################################


from pymoo.core.mutation import Mutation

def mutate_kmer(sequence, selpos, k, minlength, maxlength, seed=None):

    if seed is not None:
        random.seed(seed)

    kmer_start = selpos
    kmer_end = selpos + k
    kmer = sequence[kmer_start:kmer_end]

    mutation_type = random.choice(["substitution", "deletion", "insertion"])

    if mutation_type == "substitution":
        mutate_pos = random.randint(0, k - 1)
        bases = ["A", "U", "G", "C"]
        bases.remove(kmer[mutate_pos])  # Avoid replacing with the same base
        mutated_base = random.choice(bases)
        kmer = kmer[:mutate_pos] + mutated_base + kmer[mutate_pos + 1:]

    elif mutation_type == "deletion" and len(sequence) > minlength:
        mutate_pos = random.randint(0, k - 1)
        sequence = sequence[:kmer_start + mutate_pos] + sequence[kmer_start + mutate_pos + 1:]

    elif mutation_type == "insertion" and len(sequence) < maxlength:
        mutate_pos = random.randint(0, k)
        mutated_base = random.choice(["A", "U", "G", "C"])
        sequence = sequence[:kmer_start + mutate_pos] + mutated_base + sequence[kmer_start + mutate_pos:]

    else:
        # If deletion or insertion cannot be performed, fallback to substitution
        mutate_pos = random.randint(0, k - 1)
        bases = ["A", "U", "G", "C"]
        bases.remove(kmer[mutate_pos])  # Avoid replacing with the same base
        mutated_base = random.choice(bases)
        kmer = kmer[:mutate_pos] + mutated_base + kmer[mutate_pos + 1:]

    if mutation_type == "substitution":
        sequence = sequence[:kmer_start] + kmer + sequence[kmer_end:]

    return sequence    



BASES = ["A", "C", "G", "U"]



class AdaptiveMutation(Mutation):
    def __init__(self, directory, k=4, minlength=10, maxlength=30, leakage_rate=0.1):
        super().__init__()
        self.directory = directory
        self.k = k  
        self.minlength = minlength
        self.maxlength = maxlength
        self.leakage_rate = leakage_rate

    def _do(self, problem, X, **kwargs):
        
        # laod current mutprobs
        with open(f'{self.directory}summutprobs.pkl', 'rb') as f:
            mutprobs = pickle.load(f)

        for i in range(X.shape[0]):

            sequence_length = int(X[i, 0])
            raw_sequence = X[i, 1:sequence_length + 1].astype(int)
            sequence = ''.join([BASES[int(b)] for b in raw_sequence]) # Decode sequence

            # Step 1: Filter mutprobs to include only k-mers that exist in the sequence
            kmers_in_sequence = [
                kmer for kmer in mutprobs.keys() if kmer in sequence
            ]
            filtered_mutprobs = {kmer: mutprobs[kmer] for kmer in kmers_in_sequence}

            # Step 2: Select a k-mer and return its position in the original sequence
            if random.random() < self.leakage_rate:
                # Leakage: Do no follow mutprobs
                selected_kmer = random.choice(list(filtered_mutprobs.keys()))
            else:
                highest_prob_kmers = [
                    kmer for kmer, prob in filtered_mutprobs.items() if prob > 0.5
                ]
                if highest_prob_kmers:
                    selected_kmer = random.choice(highest_prob_kmers)  # Random among highest
                else:
                    selected_kmer = random.choice(list(filtered_mutprobs.keys()))

            # Step 3: Find the position(s) of the selected k-mer in the sequence
            positions = [pos for pos in range(len(sequence) - len(selected_kmer) + 1)
                        if sequence[pos:pos + len(selected_kmer)] == selected_kmer]

  
            selpos = random.choice(positions) # Randomly choose one
            
            # Mutate
            mutated_sequence = mutate_kmer(
                sequence, 
                selpos, 
                k=self.k,
                minlength=self.minlength, 
                maxlength=self.maxlength
            )
            
            encoded_sequence = [BASES.index(base) for base in mutated_sequence]
            encoded_sequence = np.pad(encoded_sequence, (0, max(0, self.maxlength - len(encoded_sequence))), constant_values=0)

            X[i, 0] = min(len(mutated_sequence), self.maxlength)
            X[i, 1:1 + self.maxlength] = encoded_sequence[:self.maxlength]
            
            #X[i, :] = [BASES.index(base) for base in mutated_sequence]  # Encode mutant

        return X