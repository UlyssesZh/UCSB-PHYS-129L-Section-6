#!/usr/bin/env python

from numpy import insert, cumsum, sqrt, linspace, exp, full_like, geomspace, rint, mean, std, array
from numpy.random import seed, normal
from matplotlib import pyplot as plt

# a
def wiener_process(T, N, random_seed=1108):
	seed(random_seed)
	dW = normal(0, sqrt(T/N), size=N)
	W = insert(cumsum(dW), 0, 0)
	return linspace(0, T, N+1), W

def strato_int(f, t, x):
	integrand = f(t, x)
	dx = x[1:] - x[:-1]
	return insert(cumsum((integrand[1:] + integrand[:-1]) / 2 * dx), 0, 0)

def ito_int(f, t, x):
	integrand = f(t, x)
	dx = x[1:] - x[:-1]
	return insert(cumsum(integrand[:-1] * dx), 0, 0)

def ito_to_strato(f, t, x, d=1e-6):
	fx = ((f(t, x+d) - f(t, x-d)) / (2*d))[:-1]
	dt = t[1:] - t[:-1]
	return ito_int(f, t, x) + insert(cumsum(fx * dt), 0, 0) / 2

t, x = wiener_process(10, 100)
args = (lambda t, x: x, t, x)
plt.plot(t, ito_int(*args), label='Ito')
plt.plot(t, strato_int(*args), label='Strato')
plt.plot(t, ito_to_strato(*args), label='Ito to Strato')
plt.legend()
plt.show()

# b
def geometric_brownian(mu, sigma, T, N, seed=1108):
	t, W = wiener_process(T, N, seed)
	return t, exp(mu*t + sigma*W)
plt.plot(*geometric_brownian(0.1, 0.2, 10, 100))
plt.show()

# c
t, x = geometric_brownian(0.1, 0.2, 10, 100)
args = (lambda t, x: full_like(t, 1), t, x)
plt.plot(t, ito_int(*args))
plt.show()

# d
plt.plot(t, strato_int(*args))
plt.show()

# e
ito_mean = []
strato_mean = []
ito_std = []
strato_std = []
Ns = rint(geomspace(10, 1e4, 6)).astype(int)
for N in Ns:
	ito_final = []
	strato_final = []
	for i in range(500):
		t, x = geometric_brownian(0.1, 0.2, 10, N, seed=i)
		args = (lambda t, x: full_like(t, 1), t, x)
		ito_final.append(ito_int(*args)[-1])
		strato_final.append(strato_int(*args)[-1])
	ito_mean.append(mean(ito_final))
	strato_mean.append(mean(strato_final))
	ito_std.append(std(ito_final))
	strato_std.append(std(strato_final))
fig, axs = plt.subplots(2, 2)
axs[0,0].plot(Ns, ito_mean)
axs[0,0].set_xscale('log')
axs[0,0].set_title('Ito Mean')
axs[0,1].plot(Ns, strato_mean)
axs[0,1].set_xscale('log')
axs[0,1].set_title('Strato Mean')
axs[1,0].plot(Ns, ito_std)
axs[1,0].set_xscale('log')
axs[1,0].set_title('Ito Std')
axs[1,1].plot(Ns, strato_std)
axs[1,1].set_xscale('log')
axs[1,1].set_title('Strato Std')
plt.show()

# fc
t, x = geometric_brownian(0.1, 0.2, 10, 100)
args = (lambda t, x: x**2, t, x)
plt.plot(t, ito_int(*args))
plt.show()

# fd
plt.plot(t, strato_int(*args))
plt.show()

# fe
ito_mean = []
strato_mean = []
ito_std = []
strato_std = []
Ns = rint(geomspace(10, 1e4, 6)).astype(int)
for N in Ns:
	ito_final = []
	strato_final = []
	for i in range(500):
		t, x = geometric_brownian(0.1, 0.2, 10, N, seed=i)
		args = (lambda t, x: x**2, t, x)
		ito_final.append(ito_int(*args)[-1])
		strato_final.append(strato_int(*args)[-1])
	ito_mean.append(mean(ito_final))
	strato_mean.append(mean(strato_final))
	ito_std.append(std(ito_final))
	strato_std.append(std(strato_final))
fig, axs = plt.subplots(2, 2)
axs[0,0].plot(Ns, ito_mean)
axs[0,0].set_xscale('log')
axs[0,0].set_title('Ito Mean')
axs[0,1].plot(Ns, strato_mean)
axs[0,1].set_xscale('log')
axs[0,1].set_title('Strato Mean')
axs[1,0].plot(Ns, ito_std)
axs[1,0].set_xscale('log')
axs[1,0].set_title('Ito Std')
axs[1,1].plot(Ns, strato_std)
axs[1,1].set_xscale('log')
axs[1,1].set_title('Strato Std')
plt.show()

# The bare trajectory is the same for different integraters,
# but different for the functional dynamics.
# Also, the statistics of the bare trajectory do not seem to depend on N,
# but those of the functional dynamics do.

# g
# Actually C(t, tau) = <F(t) F(t+tau)> depends on t in addition to tau
# (because <F(t)> ~ exp(3 mu t) depends on t, also notice that F(0) = 0 always),
# but the worksheet seems to assume that it does not.
# I will adopt C(tau) = <F(tau)> instead (i.e. C(tau) ~ C(0, tau) / F(0)).
dt = 0.1
N = 300
t = array([5, 10, 20, 30])

t_indices = rint(t / dt).astype(int)
corr = []
for i in range(500):
	tt, x = geometric_brownian(0.1, 0.2, dt*N, N, seed=i)
	f = strato_int(lambda t, x: x**2, tt, x)
	corr.append([f[j] for j in t_indices])
corr = mean(corr, axis=0)
fig, ax = plt.subplots()
ax.plot(t, corr)
ax.set_yscale('log')
plt.show()
