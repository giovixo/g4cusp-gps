import numpy as np
from scipy.integrate import solve_ivp
import sympy as sp

"""
radiationSolver
Library for the numerical solution of the Bateman equation for radioactive decay-
with applications to Geant4 simulations

Version 1.0 - 2023-12-04

(c) Riccardo Campana - INAF/OAS
    riccardo.campana@inaf.it
"""

def bateman_equation(t, y, lambdas):
    """
    Bateman equation for a system of radioactive decay modes.
    
    Input:
        t (float): Time
        y (array): Array of concentrations at time t (can be arbitrarily scaled)
        lambdas (array): Array of decay constants for each decay mode
        
    Output:
        dydt (array): Array of derivatives of concentrations at time t
    """
    n = len(y)  # Number of decay modes
    dydt = np.zeros(n)
    
    # Compute derivatives for each decay mode
    for i in range(n):
        dydt[i] = -lambdas[i] * y[i]
        
        # Add contributions from other decay modes
        for j in range(n):
            if j == i-1:
                dydt[i] += lambdas[j] * y[j]
    
    return dydt



def solve_bateman_equation(lambdas, y0, t_span, t_eval=None):
    """
    Solve the Bateman equation numerically.
    
    Input:
        lambdas (array): Array of decay constants for each decay mode
        y0 (array): Initial concentrations
        t_span (tuple): Time span (start, end)
        t_eval (array, optional): Times in which to evaluate the solution
    Output:
        sol (OdeSolution): Solution() object containing the result, from scipy.integrate.solve_ivp
    """
    if t_eval is None:
        t_eval = np.linspace(t_span[0], t_span[1], 100)
    sol = solve_ivp(lambda t, y: bateman_equation(t, y, lambdas), t_span, y0, t_eval=t_eval, method='DOP853')
    return sol
  

def solve_bateman_equation_Onega(decay_constants, initial_concentrations, times):
    """
    Solve the Bateman equations using the Onega matrix approach
    Ref:    R.J. Onega, Am. J. Phys. 37 (1969) 1019.
    See also:   Amaku et al. Computer Physics Communications 181 (2010) 21–23
                Hauf et al. arXiv:1307.0996 (2013)
    """

    #Construct the lambda matrix from the decay constants
    n= len(decay_constants)
    L = np.zeros((n,n))
    
    for i,l in enumerate(decay_constants):
        col = np.zeros(n)
        col[i] = -l
        if i < n-1:
            col[i+1] = l
        L[:,i] = col
    # Get eigenvalues and eigenvectors of the lambda matrix
    # We have to use NOT-NORMALISED eigenvectors, so we use sympy instead of numpy
    mat_L = sp.Matrix(L)
    eigvals = np.real(np.array([complex(x) for x in mat_L.eigenvals().keys()], dtype=complex))
    eigvec = np.real(np.array([list(x[2][0]) for x in mat_L.eigenvects()], dtype=complex))
    # Derive C and C^-1 matrices
    C = eigvec.T
    Cm = np.linalg.inv(C)

    concentrations = np.zeros((n, len(times)))
    
    for j,t in enumerate(times):
        expLd = np.diag(np.exp(eigvals*t))
        a = np.dot(C,np.dot(expLd,Cm))
        concentrations[:,j] = np.dot(a, initial_concentrations)
    return concentrations
