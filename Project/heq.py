
import qutip as qt
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
import string
import time

Ec = 1
level = 4
Nc = 50
Ej = np.array([0.1, 1.0, 5.0, 10.0])
ng = np.arange(-2, 2, 0.01)

# Operators
N = qt.charge(Nc)  # charge operator
Ic = qt.qeye(2*Nc+1)  # Identity operator
T = qt.tunneling(2*Nc+1)  # tunneling operator

def Hamil(ng, Ej):
    '''
    Function that calculates eigenvalues of 
    Hamiltonian in charge basis.
    Input: gate charge ng, Josephson energy Ej
    Output: lowest three eigenenergies
    '''
    H = 4*Ec*(N - ng*Ic)**2 + 0.5*Ej*T
    return H.eigenenergies(eigvals=level)