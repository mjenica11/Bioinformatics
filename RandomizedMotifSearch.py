#!/usr/bin/env python3

from random import randint

def HammingDistance(Seq1, Seq2):
    if len(Seq1) != len(Seq2):
        raise ValueError('Sequences are not equal length.')
    # make tuples at matching indices and return the number of tuples with diff values
    return sum(c1 != c2 for c1, c2 in zip(Seq1, Seq2))

def Score_Motifs(Motifs):
    score = 0
    for i in range(KMER):
        # make a col at each index by joining the characters at index i from each Motif 
        Motif  = ''.join([Motifs[j][i] for j in range(len(Motifs))])
        # make a string of each nt, each as long as Motif  
        # count the number of diffs between each nt string and each Motif; store the min value
        score += min([HammingDistance(Motif, x*len(Motif)) for x in 'ATCG'])
    return score

def Profile(Motifs):
    Prof = []
    for i in range(KMER):
        col = ''.join([Motifs[j][i] for j in range(NUM_KMER)])
        # count the freq of each nt in each col and divide by the length of the col
        Prof.append([col.count(nt)/len(col) for nt in 'ATCG'])
    return Prof

def Pseudocount_Profile(Motifs):
    Prof = []
    for i in range(KMER):
        # make a col at each index by joining the characters at index i from each Motif 
        col = ''.join([Motifs[j][i] for j in range(NUM_KMER)])
        # count the freq of each nt in each col and divide by the length of the col; 
        # add 1 to each col + 4 to the denominator
        Prof.append([col.count(nt)+1/len(col)+4 for nt in 'ATCG'])
    return(Prof)

def Most_Likely_Kmer(DNA, KMER, Prof):
    # make a dict enumerating each nt
    NT_Dct = {nt:index for index, nt in enumerate('ATCG')}
    # initialize the maximum probability
    Max_Prob = [-1, None]
    # Compute the probability of each KMER-mer, store store it if it's currently a maximum
    for i in range(len(DNA)-KMER + 1):
        Current_Prob = 1
        for j, nt in enumerate(DNA[i:i+KMER]):
            Current_Prob *= Prof[j][NT_Dct[nt]]
        if Current_Prob > Max_Prob[0]:
            Max_Prob = [Current_Prob, DNA[i:i+KMER]]
    return Max_Prob[1]

def Profile_To_Motif(Profile, DNA, KMER):
    return [Most_Likely_Kmer(seq, KMER, Profile) for seq in DNA]

def Random_Motif_Search(DNA, KMER, NUM_KMER):
    # Randomly generate list of integers with length KMER  
    Random_Ints = [randint(0, len(DNA[0])-KMER) for a in range(NUM_KMER)]
    # For each seq given (first item in tuples generates with enumerate),
    # make a list of Motifs of length KMER begining at a randomly generated index 
    Motifs = [DNA_Lst[i][r:r+KMER] for i,r in enumerate(Random_Ints)]

    # Initialize the best score as a score higher than the highest possible score
    Best_Score = [Score_Motifs(Motifs), Motifs]

    # Iterate Motifs
    while True:
        Current_Profile = Pseudocount_Profile(Motifs)
        Motifs = Profile_To_Motif(Current_Profile, DNA_Lst, KMER)
        Current_Score = Score_Motifs(Motifs)
        if Current_Score < Best_Score[0]:
            Best_Score = [Current_Score, Motifs]
        else:
            return Best_Score

if __name__ == '__main__':

    with open('rosalind_ba2f.txt') as input_data:
        # split the first line and read in as integers
        KMER, NUM_KMER = map(int, input_data.readline().split())
        # read in sequences line by line and store as list
        DNA_Lst = [line.strip() for line in input_data.readlines()]

    # Initialize the best scoring Motif as a score higher than the highest possible score
    Best_Motifs = [KMER*NUM_KMER, None]

    # Repeat the randomized Motif search 1000 times
    for repeat in range(1000):
        Current_Motifs = Random_Motif_Search(DNA_Lst, KMER, NUM_KMER)
        if Current_Motifs[0] < Best_Motifs[0]:
            Best_Motifs = Current_Motifs

    # Print and save the answer
    with open('3Ans.txt', 'w') as output_data:
        output_data.write('\n'.join(Best_Motifs[1]))
