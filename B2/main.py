#!/usr/bin/env python

from numpy import sqrt, insert, cumsum, linspace
from numpy.random import seed, normal
from matplotlib import pyplot as plt

def wiener_process(T, N, random_seed=1108):
	seed(random_seed)
	dW = normal(0, sqrt(T/N), size=N)
	W = insert(cumsum(dW), 0, 0)
	return linspace(0, T, N+1), W

T = 10
N = 100
v0 = 2
D = 0.1
gamma = 1

dt = T/N
t, W = wiener_process(10, 100)
W2 = W**2
eta = (W[1:] - W[:-1]) / dt
v = [v0]
for i in range(1, N):
	dvdt = -gamma*v[-1] + sqrt(2*D)*eta[i]
	v.append(v[-1] + dvdt*dt)
plt.plot(t[:-1], v)
plt.show()
