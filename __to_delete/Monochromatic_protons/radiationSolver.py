import numpy as np
from scipy.integrate import solve_ivp
import sympy as sp

"""
radiationSolver
Library for the numerical solution of the Bateman equation for radioactive decay
with applications to Geant4 simulations.

Version 1.1 - 2025-06-10

(c) Riccardo Campana - INAF/OAS
    riccardo.campana@inaf.it

Changelog v1.1:
  - FIX (critical): solve_bateman_equation_Onega fails silently for degenerate
    decay constants (two isotopes with the same lambda). Sympy introduces tiny
    imaginary parts in the roots of the characteristic polynomial, making C
    near-singular (det ~ 1e-33). The inversion then produces values of order
    1e17 instead of the correct O(1) concentrations. Added a degenerate-
    eigenvalue guard that raises ValueError with a clear message rather than
    returning garbage.

  - FIX: eigenvalues and eigenvectors were extracted via two independent sympy
    calls (eigenvals() and eigenvects()) with no pairing guarantee. Both are
    now extracted from a single eigenvects() call, ensuring eigenvectors are
    always in the same order as the eigenvalues used to build expLd.

  - FIX: the inner loop in bateman_equation is O(n^2) but only the j==i-1
    condition ever fires. Replaced with a direct conditional assignment, making
    the function O(n).

  - NOTE (design): both solvers model a strictly linear decay chain where 100%
    of the decays from isotope i feed isotope i+1. Branching ratios are not
    supported; to handle branching the caller must split each chain into
    independent linear chains weighted by their cumulative branching ratio and
    sum the results.

  - NOTE (performance): the time loop in solve_bateman_equation_Onega has been
    vectorised: C and C^-1 are computed once and applied to all time steps in
    a single matrix operation, replacing the Python-level for loop.
"""


def bateman_equation(t, y, lambdas):
    """
    Right-hand side of the Bateman ODE for a linear decay chain.

    dN_0/dt = -lambda_0 * N_0
    dN_i/dt = -lambda_i * N_i  +  lambda_{i-1} * N_{i-1}   (i > 0)

    Input:
        t       (float) : time (not used explicitly; required by solve_ivp)
        y       (array) : concentrations at time t
        lambdas (array) : decay constants [s^-1]

    Output:
        dydt (array) : derivatives of concentrations at time t
    """
    n = len(y)
    dydt = np.empty(n)

    # FIX: replaced O(n^2) double loop with O(n) direct assignments.
    # The original inner loop had 'if j == i-1' which only fired for j = i-1,
    # making the second loop completely unnecessary.
    dydt[0] = -lambdas[0] * y[0]
    for i in range(1, n):
        dydt[i] = -lambdas[i] * y[i] + lambdas[i-1] * y[i-1]

    return dydt


def solve_bateman_equation(lambdas, y0, t_span, t_eval=None):
    """
    Solve the Bateman equation numerically using scipy's DOP853 integrator.

    Input:
        lambdas (array)       : decay constants [s^-1] for each isotope in the chain
        y0      (array)       : initial concentrations
        t_span  (tuple)       : (t_start, t_end)
        t_eval  (array, opt.) : evaluation times; defaults to 100 equally-spaced points

    Output:
        sol (OdeSolution) : result from scipy.integrate.solve_ivp
    """
    if t_eval is None:
        t_eval = np.linspace(t_span[0], t_span[1], 100)
    sol = solve_ivp(
        lambda t, y: bateman_equation(t, y, lambdas),
        t_span, y0, t_eval=t_eval, method='DOP853'
    )
    return sol


def solve_bateman_equation_Onega(decay_constants, initial_concentrations, times):
    """
    Solve the Bateman equations using the Onega matrix-diagonalisation approach.

    Ref: R.J. Onega, Am. J. Phys. 37 (1969) 1019.
    See also: Amaku et al., Comput. Phys. Commun. 181 (2010) 21-23;
              Hauf et al., arXiv:1307.0996 (2013).

    N(t) = C * diag(exp(lambda_k * t)) * C^{-1} * N(0)

    where C is the matrix of (un-normalised) eigenvectors of the decay matrix L,
    and the lambda_k are its eigenvalues (all real and negative for physical chains).

    Input:
        decay_constants       (array, length n) : lambda_i = ln2 / T_{1/2,i} [s^-1]
        initial_concentrations (array, length n) : N_i(0)
        times                  (array, length m) : evaluation times [s]

    Output:
        concentrations (array, shape n x m) : N_i(t_j)

    Raises:
        ValueError : if any two decay constants are equal (degenerate eigenvalues),
                     because the diagonalisation becomes numerically ill-conditioned.
                     Use solve_bateman_equation() instead for degenerate cases.
    """
    n = len(decay_constants)

    # --- Guard: degenerate eigenvalues make C singular ---
    # FIX: check before invoking sympy to give a clear error rather than
    # silently returning concentrations of order 1e17.
    dc = np.asarray(decay_constants, dtype=float)
    tol = 1e-12 * max(np.max(np.abs(dc)), 1.0)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(dc[i] - dc[j]) < tol:
                raise ValueError(
                    f"Degenerate decay constants: lambda[{i}]={dc[i]:.6g} and "
                    f"lambda[{j}]={dc[j]:.6g} are equal to within {tol:.1e}. "
                    "The Onega method requires distinct eigenvalues. "
                    "Use solve_bateman_equation() for degenerate cases."
                )

    # --- Build the lower-triangular decay matrix L ---
    # L[i,i]   = -lambda_i
    # L[i+1,i] =  lambda_i
    L = np.zeros((n, n))
    for i, l in enumerate(dc):
        L[i, i] = -l
        if i < n - 1:
            L[i + 1, i] = l

    # --- Diagonalise with sympy to obtain exact (un-normalised) eigenvectors ---
    # FIX: extract both eigenvalues and eigenvectors from a single eigenvects()
    # call to guarantee consistent ordering. The original code used eigenvals()
    # and eigenvects() independently, with no pairing guarantee.
    mat_L = sp.Matrix(L)
    eigdata = mat_L.eigenvects()   # list of (eigenvalue, multiplicity, [vectors])

    eigvals = np.array(
        [complex(ev) for ev, _mult, _vecs in eigdata], dtype=complex
    )
    eigvec = np.array(
        [list(vecs[0]) for _ev, _mult, vecs in eigdata], dtype=complex
    )

    # Eigenvalues of a real lower-triangular matrix with real entries are real;
    # any imaginary part is numerical noise from sympy's polynomial root-finding.
    if np.max(np.abs(np.imag(eigvals))) > 1e-10:
        raise ValueError(
            "Eigenvalues have non-negligible imaginary parts "
            f"(max |Im| = {np.max(np.abs(np.imag(eigvals))):.2e}). "
            "This suggests near-degenerate decay constants. "
            "Use solve_bateman_equation() instead."
        )

    eigvals = np.real(eigvals)
    eigvec  = np.real(eigvec)

    # C: columns are eigenvectors; C^{-1}: its inverse
    C  = eigvec.T
    Cm = np.linalg.inv(C)

    # --- Vectorised time evolution ---
    # FIX: replaced Python for-loop over times with a fully vectorised operation.
    # Shape:  eigvals (n,), times (m,) -> exp_mat (n, m)
    #         C (n,n), Cm (n,n), N0 (n,) -> concentrations (n, m)
    times_arr = np.asarray(times)
    exp_mat   = np.exp(eigvals[:, np.newaxis] * times_arr[np.newaxis, :])  # (n, m)
    # concentrations = C @ (exp_mat * (Cm @ N0)[:, np.newaxis])
    Cm_N0 = Cm @ np.asarray(initial_concentrations, dtype=float)           # (n,)
    concentrations = C @ (exp_mat * Cm_N0[:, np.newaxis])                  # (n, m)

    return concentrations
