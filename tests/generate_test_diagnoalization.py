import numpy as np
from numpy import linalg

h0 = np.array([[0, -2.8284271247461903], [-2.82842712474619, 4.0]])
h2 = np.array([[0, -2.8284271247461903], [-2.82842712474619, 3.6085624472639157]])

w0, v0 = linalg.eigh(h0)
w2, v2 = linalg.eigh(h2)

def ei(x):
    from cmath import exp
    return exp(1j * x)

x1 = v0.conj().T.dot(v2)
x2 = x1.conj().T
unitary1 = np.diag(list(map(ei, w2)))

right = x1.dot(unitary1).dot(x2)
left = right.conj().T

energy_diag = np.diag(w0)

change = (left.dot(energy_diag).dot(right) - energy_diag).diagonal()

print(h2, w2, v2, unitary1, x1, x2, right, left, change, sep="\n\n")
