#!/usr/bin/env python3

from random import randint
from RandomizedMotifSearch import Score_Motifs, Most_Likely_Kmer, Pseudocount_Profile

def Gibbs_Sampler(DNA, KMER, NUM_KMER, NUM_SAMPLES):
    # Randomly generate list of integers with length KMER
    Random_Ints = [randint(0, len(DNA[0])-KMER) for a in range(NUM_KMER)]
    # Make a list of motifs of length KMER beginning at a randomly generated index
    Motifs = [DNA[i][r:r+KMER] for i,r in enumerate(Random_Ints)]
    # Initialize the best score as a score higher than the highest possible score
    Best_Score = [Score_Motifs(Motifs, KMER), Motifs]
    # Iterate through motifs
    for n in range(NUM_SAMPLES):
        r = randint(0, NUM_KMER-1)
        # Make a pseudocount profile for each motif except for one, randomly picked
        Current_Profile = Pseudocount_Profile([motif for index, motif in enumerate(Motifs) if index != r], KMER, NUM_KMER-1)
        Motifs = [Most_Likely_Kmer(DNA[index], KMER, Current_Profile) if index == r else Motif for index, Motif in enumerate(Motifs)]
        Current_Score = Score_Motifs(Motifs, KMER)
        if Current_Score < Best_Score[0]:
            Best_Score = [Current_Score, Motifs]
    return Best_Score

if __name__ == '__main__':
       
    with open('rosalind_ba2g.txt') as input_data:
        Kmer, Num_Kmer, Num_Samples = map(int, input_data.readline().split())
        DNA_Lst = [line.strip() for line in input_data.readlines()]
    
    # Initialize the best scoring motif as a score higher than the highest possible score
    Best_Motifs = [Kmer * Num_Kmer, None]
    
    # Repeat randomized motif seach 20 times
    for repeat in range(20):
        Current_Motifs = Gibbs_Sampler(DNA_Lst, Kmer, Num_Kmer, Num_Samples)
        if Current_Motifs[0] < Best_Motifs[0]:
            Best_Motifs = Current_Motifs

    # Print and save the answer
    with open('3b_Ans.txt', 'w') as Output_Data:
        Output_Data.write('\n'.join(Best_Motifs[1]))
