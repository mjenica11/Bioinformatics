#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def initialize_centers(x, num_k):
    # nrow, ncol
    N, D = x.shape
    # initialize empty array; ncol = num k, same dim as data
    centers = np.zeros((num_k, D))
    used_idx = []
    # randomly choose k number of indices from the data
    # then save them as centers in a random order
    for k in range(num_k):
        idx = np.random.choice(N)
        while idx in used_idx:
            idx = np.random.choice(N)
        used_idx.append(idx)
        centers[k] = x[idx]
    return centers

def cluster_responsibilities(centers, x, beta):
    # initialize matrix w/ ncol = num k and nrow = num data points
    N,_ = x.shape
    K, D = centers.shape
    R = np.zeros((N, K))

    for n in range(N):
        R[n] = np.exp(-beta * np.linalg.norm(centers-x[n], ord=2, axis=1))
        R /= R.sum(axis=1, keepdims=True)
    return R

def update_centers(x, r, num_k):
    N, D = x.shape
    centers = np.zeros((num_k, D))
    for k in range(num_k):
        centers[k] = r[:,k].dot(x) / r[:,k].sum()
    return centers

def square_dist(a, b):
    return (a-b) ** 2

def cost_func(x, r, centers, num_k):
    cost = 0
    for k in range(num_k):
        norm = np.linalg.norm(x - centers[k], 2)
        cost += (norm * np.expand_dims(r[:, k], axis=1)).sum()
    return cost

def soft_k_means(x, num_k, max_iters=20, beta=.1):
    centers = initialize_centers(x, num_k)
    prev_cost = 0
    for _ in range(max_iters):
        r = cluster_responsibilities(centers, x, beta)
        centers = update_centers(x, r, num_k)
        cost = cost_func(x, r, centers, num_k)
        if np.abs(cost - prev_cost) < 1e-5:
            break
        prev_cost = cost
    return centers


def generate_samples(std=1, dim=2, dist=4):
    mu0 = np.array([0,0])
    mu1 = np.array([dist, dist])
    mu2 = np.array([0, dist])
    # Number of samples per class
    Nc = 300
    # Initialize array ncol Nc nrow dim 
    x0 = np.random.randn(Nc, dim) * std + mu0
    x1 = np.random.randn(Nc, dim) * std + mu1
    x2 = np.random.randn(Nc, dim) * std + mu2
    # stack vertically
    x = np.concatenate((x0, x1, x2), axis=0)
    return x

def main():
    X = generate_samples()
    soft_k_means(x=X, num_k=3)

if __name__ == "__main__"
    main()

