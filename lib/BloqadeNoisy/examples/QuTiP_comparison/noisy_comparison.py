import qutip as qt
import numpy as np
from matplotlib import pyplot as plt

def n():
    return (1-qt.sigmaz())/2

def put(nqubits, qubit, op):
    ops = [qt.identity(2) for q in range(nqubits)]
    ops[qubit] = op
    return qt.tensor(*reversed(ops))

def state(*states):
    return qt.tensor(*[qt.basis(2, s) for s in reversed(states)])

def rectlattice(h, w, scale = 1):
    return [(scale * j, scale*i) for i in range(w) for j in range(h)]

C = np.pi * 2 * 862690

def dist(c1, c2):
    return np.sqrt((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2)

def rydberg_h(atoms, Omega, delta):
    N = len(atoms)
    sumOfX = sum([put(N, i, qt.sigmax()) for i in range(N)])
    sumOfn = sum([put(N, i, n()) for i in range(N)])
    ryd_op = []
    for i in range(N):
        for j in range(i+1, N):
            ryd_op.append(
                C/dist(atoms[i], atoms[j])**6 * put(N, i, n()) * put(N, j, n())
            )
    return Omega/2 * sumOfX - delta * sumOfn + sum(ryd_op)

Omega = 2*np.pi
R = (C/Omega)**(1/6)
N = 3
atoms = [(0, 0),(8, 0),(18,0)]
h = rydberg_h(atoms, Omega, 0)

rate = 1/10

destroy = qt.basis(2,0) * qt.basis(2,1).dag()
c_ops = [np.sqrt(rate) * put(N, i, destroy) for i in range(N)]

times = np.arange(0, 4.0, 1e-2)

sol_noisy = qt.mesolve(
    h, 
    state(*np.zeros(N, dtype = int)),
    times,
    c_ops = c_ops,
    e_ops = [put(N, i, n()) for i in range(3)]
)

sol = qt.sesolve(
    h,
    state(*np.zeros(N, dtype = int)),
    times,
    e_ops = [put(N, i, n()) for i in range(3)]
)

np.savetxt("./3q_expect_nonoise.csv", sol.expect,  header = "", delimiter=",")
np.savetxt("./3q_expect.csv", sol_noisy.expect,  header = "", delimiter=",")