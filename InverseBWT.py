#!/usr/env/python3

def Inverse_BWT(string):
    l = len(string)

    # store in dict the number of times each item in the string occurs
    def Rank_Items(i):
        # set default value to 0
        d[string[i]] = d.get(string[i],0) + 1
        return d[string[i]]

    d = {}
   
    # make a list of tuples; char in string, idx + 1, idx
    Ranked_Characters =  [(string[i], Rank_Items(i), i) for i in range(l)]

    # sort by first item in tuple (character in string)
    Sorted_Ranking = sorted(Ranked_Characters)

    # print both lists of tuples
    for i in range(l):
        r = str(Sorted_Ranking[i]) + (' ' * 5) + str(Ranked_Characters[i])
        print(r)

    # decode string
    i = 0
    Decoded = ''
    for j in range(l):
        # get the last item in each tuple (original index)
        i = Sorted_Ranking[i][2]
        # get the first item in each tuple (character) and add to string
        Decoded += Sorted_Ranking[i][0]
    return Decoded

if __name__ == "__main__":
    # read in string to be Decoded
    with open('rosalind_ba9j.txt', 'r') as input_data:
        STRING = input_data.read().replace('\n', '')

    # do inverse BWT
    IBWT = Inverse_BWT(STRING)

    # print and save result
    with open("hw5_result.txt", "wt", encoding='utf-8') as output_data:
        output_data.write(IBWT)
