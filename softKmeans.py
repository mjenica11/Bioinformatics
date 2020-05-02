#!/usr/bin/env python3

import sys
import numpy as np

def decide_centers(data, num_k):
    # nrow, ncol
    N, D = data.shape
    # initialize empty array; ncol = num k, same dim as data
    centers = np.zeros((num_k, D))
    old_idx = []
    # randomly choose k number of indices from the data
    # then save them as centers in a random order
    for k in range(num_k):
        idx = np.random.choice(N)
        while idx in old_idx:
            idx = np.random.choice(N)
        old_idx.append(idx)
        centers[k] = data[idx]
    return centers

def centers_to_clusters(centers, data, beta):
    # initialize matrix w/ ncol = num k and nrow = num data points
    N,_ = data.shape
    num_k,_ = centers.shape
    resp_mat = np.zeros((N, num_k))
    # generate responsibility matrix
    for n in range(N):
        resp_mat[n] = np.exp(-beta * np.linalg.norm(centers-data[n], ord=2, axis=1))
        resp_mat /= resp_mat.sum(axis=1, keepdims=True)
    return resp_mat

def clusters_to_centers(data, resp_mat, num_k):
    N, dim = data.shape
    centers = np.zeros((num_k, dim))
    # update center to the weighted center of gravity
    for k in range(num_k):
        centers[k] = resp_mat[:,k].dot(data) / resp_mat[:,k].sum()
    return centers

def soft_kMeans(data, num_k, beta, max_iters=100):
    centers = decide_centers(data, num_k)
    for _ in range(max_iters):
        resp = centers_to_clusters(centers, data, beta)
        centers = clusters_to_centers(data, resp, num_k)
    return centers

def main():
    with open('rosalind_ba8d.txt') as f:
        k, m = map(int, f.readline().split())
        b = float(f.readline().strip())
        mat = [[float(num) for num in line.split(' ')] for line in f]
        mat = np.array(mat)

    Centers = soft_kMeans(data=mat, num_k=k, beta=b)
    Centers = np.flipud(Centers)
    txt = np.savetxt('res.txt', Centers, fmt="%.3f")
    out = open('res.txt', 'r')
    print(out.read())

if __name__ == "__main__":
    main()

