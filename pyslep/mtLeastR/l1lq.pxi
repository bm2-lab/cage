import numpy.matlib as np
cimport numpy as np
from libc.math cimport floor

cdef extern from "eppMatrix.h":
    void _eppMatrix(double *X, double * V, int k, int n, double rho, double p)

cdef np.ndarray eppMatrix(double[:, :] v, int k, int n, double rho, double p):
    cdef double[:, :] x = np.zeros((n, k), dtype=np.double)
    _eppMatrix(&x[0,0], &v[0,0], k, n, rho, p)
    return np.asmatrix(x).T

cdef extern from "ep21R.h":
    void _ep21R(double * x, double *t, double * u, double * v, int n, int k)

cdef tuple ep21R(double[:, :] u, double[:, :] v, int n, int k):
    cdef double[:, :] x = np.zeros((u.shape[0], u.shape[1]), dtype=np.double)
    cdef double[:, :] t = np.zeros((v.shape[0], v.shape[1]), dtype=np.double)
    _ep21R(&x[0,0], &t[0,0], &u[0,0], &v[0,0], n, k)
    return (np.asmatrix(x), np.asmatrix(t))

cdef double initFactor(double x_norm, double[:, :] _Ax , double[:, :] _y, double z):
    cdef np.ndarray Ax = np.asmatrix(_Ax)
    cdef np.ndarray y = np.asmatrix(_y)
    cdef double ratio = (Ax.T * y - z * x_norm) / (Ax.T * Ax)
    return ratio
    
cdef void sll_opts(dict opts):
    # Starting point
    if 'init' in opts:
        if (opts['init'] != 0) and (opts['init'] != 1) and (opts['init'] != 2):
            opts['init'] = 0
        if ('x0' not in opts) and (opts['init'] == 1):
            opts['init'] = 0
    else:
        opts['init'] = 0
    # Termination
    if 'maxIter' in opts:
        if opts['maxIter'] < 1:
            opts['maxIter'] = 10000
    else:
        opts['maxIter'] = 10000

    if 'tol' not in opts:
        opts['tol'] = 1e-3

    if 'tFlag' in opts:
        if opts['tFlag'] < 0:
            opts['tFlag'] = 0
        elif opts['tFlag'] > 5:
            opts['tFlag'] = 5
        else:
            opts['tFlag'] = int(floor(opts['tFlag']))
    else:
        opts['tFlag'] = 0
    # Normalization
    if 'nFlag' in opts:
        if (opts['nFlag']) != 1 and (opts['nFlag'] != 2):
            opts['nFlag'] = 0
    else:
        opts['nFlag'] = 0
    # Regularization
    if 'rFlag' in opts:
        if opts['rFlag'] != 1:
            opts['rFlag'] = 0
    else:
        opts['rFlag'] = 0
    # Method (Line Search)
    if 'lFlag' in opts:
        if opts['lFlag'] != 1:
            opts['lFlag'] = 0
    else:
        opts['lFlag'] = 0
        
    if 'mFlag' in opts:
        if opts['mFlag'] != 1:
            opts['mFlag'] = 0
    else:
        opts['mFlag'] = 0


        
        

        
