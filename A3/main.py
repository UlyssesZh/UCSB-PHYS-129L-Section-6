#!/usr/bin/env python

from numpy import exp, sin, pi, tile, eye, zeros, ones, arange
from numpy.linalg import eigh, norm, matrix_power
from scipy.linalg import null_space, logm
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt

# 1
# copied from section 4
def is_nn_hopping(n, state1, state2):
	if state1 == state2:
		return 0
	d = state1 ^ state2 # find sites where the spin differs
	if state1&d==0 or state2&d==0: # one of the states has spin up on those sites
		return 0
	return int(d%3==0 and d//3&d//3-1==0) + int(d==(1<<n-1)+1) # condition 1: adjacent; condition 2: periodic

def get_spin(n, state, i):
	return 1/2 - (state >> i%n & 1)

def xxx_element(n, state1, state2):
	result = 0
	if state1 == state2:
		for i in range(n):
			result += 1/4 - get_spin(n, state1, i) * get_spin(n, state1, i+1)
	result -= is_nn_hopping(n, state1, state2)/2
	return result

def xxx_hamiltonian(n):
	return [[xxx_element(n, i, j) for j in range(1<<n)] for i in range(1<<n)]

def boltzman(energies, temp):
	numerator = exp(-energies / temp)
	return numerator / sum(numerator)

def transition_prob_matrix(n, temp=1):
	hamil = xxx_hamiltonian(n)
	energies, basis = eigh(hamil)
	prob = boltzman(energies, temp)
	basis_amp = (basis.conj() * basis).real
	return basis_amp.T @ tile(prob, (1<<n, 1)) @ basis_amp

N = 3
P = transition_prob_matrix(N)
print(P)

# 2
stationary = null_space(P.T - eye(1<<N)).T
stationary /= stationary.sum(axis=1, keepdims=True)
if len(stationary) != 1:
	raise ValueError('Stationary distribution is not unique.')
stationary = stationary[0]
print(stationary)

# 3
def iteration(transition, initial, max_iter=1000, tol=1e-6):
	current = initial
	for i in range(max_iter):
		next = current @ transition
		if norm(next - current) < tol:
			return next
		current = next
	raise ValueError('Did not converge in {} iterations.'.format(max_iter))

initial = zeros(1<<N)
initial[0b000] = 1
print(iteration(P, initial))

initial = zeros(1<<N)
initial[0b000] = 0.5
initial[0b101] = 0.5
print(iteration(P, initial))

initial = ones(1<<N) / (1<<N)
print(iteration(P, initial))

# 4
def transition_prob_matrix_magnon(n, temp=1):
	k = arange(n)
	energies = 2*sin(pi*k/n)**2
	prob = boltzman(energies, temp)
	return tile(prob, (n, 1))

P = transition_prob_matrix_magnon(N)
print(P)

# The difference between the transition matrix in site basis and in magnon basis
# is that they have different dimensions.
# The probability in the magnon transition matrix
# is the conditional probility for a magnon to transit to another magnon,
# neglecting the possibility of transiting to non-magnon states.

# 5
stationary = null_space(P.T - eye(N)).T
stationary /= stationary.sum(axis=1, keepdims=True)
if len(stationary) != 1:
	raise ValueError('Stationary distribution is not unique.')
stationary = stationary[0]
print(stationary)

# 6
initial = zeros(N)
initial[1] = 1
print(iteration(P, initial))

initial = zeros(N)
initial[1] = 0.5
initial[2] = 0.5 # According to the worksheet, it should be initial[4], but N=3 forbids such large k.
print(iteration(P, initial))

initial = ones(N) / N
print(iteration(P, initial))

# 7
steps = 100
Q = logm(matrix_power(P, steps))/steps
initial = zeros(N)
initial[1] = 1
solution = solve_ivp(lambda t, y: Q.T @ y, (0, 20), initial)
for k in range(N):
	plt.plot(solution.t, solution.y[k], label=r'$\pi_{}$'.format(k))
plt.xlabel('$t$')
plt.legend()
plt.show()
