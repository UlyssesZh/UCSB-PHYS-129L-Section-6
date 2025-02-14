# Section 6

## A

### Task 1

Zero, because a random variable taking a particular value in a continuous distribution is zero.

Assume the question is asking for the probability that the nearest star is **within** distance $R$.
The probability that $k$ stars are in a volume $V$ is given by the Poisson distribution $P(k,V)=(nV)^ke^{-nV}/k!$.
The probability that the nearest star is within distance $R$ is
one minus the probability that there are no stars in the sphere with radius $R$,
so it is

```math
1 - P(0,4\pi R^3/3) = 1 - e^{-4\pi nR^3/3}.
```

### Task 2

```math
\tilde x(\omega) = \frac{F\delta(\omega-\omega_{\mathrm f})}{-\omega^2+i\gamma\omega+\omega_0^2}.
```

```math
x(t) = \frac{Fe^{i\omega_{\mathrm f}t}}{-\omega_{\mathrm f}^2+i\gamma\omega_{\mathrm f}+\omega_0^2}.
```

```math
\int_{t=0}^{2\pi/\omega_{\mathrm f}}\mathrm{Re}\,F\,\mathrm{Re}\,dx
= \int_0^{2\pi/\omega_{\mathrm f}}dt
\,\mathrm{Re}(Fe^{i\omega_{\mathrm f}t})
\,\mathrm{Re}\,\frac{i\omega_{\mathrm f}Fe^{i\omega_{\mathrm f}t}}{-\omega_{\mathrm f}^2+i\gamma\omega_{\mathrm f}+\omega_0^2}
= \frac{\pi\omega_{\mathrm f}\gamma F}{(\omega_0-\omega_{\mathrm f})^2+\gamma^2\omega_{\mathrm f}^2}.
```

### Task 3

*On the worksheet, this task is titled task 2.*

See `A3/main.py`.

## B

### Task 1

See `B1/main.py`.

### Task 2

See `B2/main.py`.
