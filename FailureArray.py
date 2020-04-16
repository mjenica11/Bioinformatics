#!/usr/bin/env/ python3

from Bio import SeqIO

fasta = SeqIO.read('rosalind_kmp.txt','fasta')
Seq = str(fasta.seq)

def LPSArray(pat):
    M = len(pat)
    lps = [0]*M
    i = 0
    j = 1

    while j < M:
        if i == 0:
            if pat[i] == pat[j]:
                lps[j] = i + 1
                i += 1
                j += 1
            else:
                lps[j] = 0
                j += 1
        elif i != 0: 
            if pat[i] == pat[j]:
                lps[j] = i + 1
                i += 1
                j += 1
            else: 
                i = lps[i - 1]
                lps[j] = 0
    return lps

res = LPSArray(pat=Seq)
print(*res)
