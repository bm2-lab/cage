from __future__ import division
import numpy.matlib as np
import numpy.linalg as npl
cimport numpy as np
include "l1lq.pxi"

cpdef tuple mtLeastR(double[:,:] _A, double[:,:] _y, double z, dict opts):
    # Verify and initialize the parameters
    cdef np.ndarray A = np.asmatrix(_A)
    cdef np.ndarray y = np.asmatrix(_y)
    cdef np.ndarray x, funVal, ValueL
    cdef np.ndarray ind
    cdef np.ndarray mu, nu
    cdef long m, n, k
    cdef double q
    cdef double lamda
    cdef long i, j

    m = A.shape[0]
    n = A.shape[1]
    if y.size != m:
        raise ValueError('Check the length of y!')

    if z <= 0:
        raise ValueError('z should be positive!')
    
    sll_opts(opts)
    
    # Detailed initialization
    if 'ind' not in opts:
        raise IndexError('In mtLeastR, .ind should be specified')
    else:
        ind = opts['ind']
        k = ind.size - 1
        if ind[k] != m:
            raise ValueError('Check opts.ind')

    if 'q' not in opts:
        q = 2
        opts['q'] = 2
    else:
        q = opts['q']
        if q < 1:
            raise ValueError('q should be larger than 1')

    # Normalization
    cdef np.ndarray ind_zero
    if opts['nFlag'] != 0:
        if 'mu' in opts:
            mu = opts['mu']
            if mu.shape[1] != n:
                raise ValueError('Check the input .mu')
        else:
            mu = A.mean(axis=0)

        if opts['nFlag'] == 1:
            if 'nu' in opts:
                nu = opts['nu']
                if nu.shape[0] != n:
                    raise ValueError('Check the input .nu')
            else:
                nu = np.power(np.power(A, 2).sum(axis=0)/m, 0.5)
                nu = nu.T
        elif opts['nFlag'] == 2:
            if 'nu' in opts:
                nu = opts['nu']
                if nu.shape[0] != m:
                    raise ValueError('Check the input .nu')
            else:
                nu = np.power(np.power(A, 2).sum(axis=1)/n, 0.5)
        ind_zero = np.absolute(nu) <= 1e-10
        nu[ind_zero] = 1

    # Starting point initialization
    cdef np.ndarray ATy = np.zeros((n, k))
    cdef np.ndarray ind_i_m
    cdef np.ndarray ind_i
    cdef long m_i
    cdef np.ndarray tt
    cdef np.ndarray invNu
    cdef np.ndarray mu_invNu
    cdef double lamda_max = 0
    cdef np.ndarray Ax = np.zeros((m, 1))
    cdef double x_norm = 0
    cdef double ratio
    
    for i in xrange(k):
        ind_i_m = np.arange(ind[i]+1, ind[i+1]+1)
        ind_i = ind_i_m - 1
        if opts['nFlag'] == 0:
            tt = A[ind_i, :].T * y[ind_i, 0]
        elif opts['nFlag'] == 1:
            tt = A[ind_i, :].T * y[ind_i, 0] - y[ind_i, 0].sum() * mu.T
            tt = tt / nu[:, 0]
        else:
            invNu = y[ind_i, 0] / nu[ind_i, 0]
            tt = A[ind_i, :].T * invNu - invNu.sum() * mu.T
        ATy[:, i] = tt

    if opts['rFlag'] == 0:
        lamda = z
    else:
        if (z < 0) or (z > 1):
            raise ValueError('opts.rFlag=1, and z should be in [0,1]')
        q_bar = 0
        if q == 1:
            q_bar = np.inf
        elif q >= 1e6:
            q_bar = 1
        else:
            q_bar = q / (q - 1)
            
        for i in xrange(n):
            lamda_max = np.maximum(lamda_max, npl.norm(np.squeeze(np.asarray(ATy[i, :])), q_bar))
        lamda = z * lamda_max

    if opts['init'] == 2:
        x = np.zeros((n, k))
    else:
        if 'x0' in opts:
            x = opts['x0']
            if (x.shape[0] != n) or (x.shape[1] != k):
                raise ValueError('Check the input .x0')
        else:
            x = ATy

    
    for i in xrange(k):
        ind_i_m = np.arange(ind[i]+1, ind[i+1]+1)
        ind_i = ind_i_m - 1
        m_i = ind[i+1] - ind[i]

        if opts['nFlag'] == 0:
            Ax[ind_i, 0:1] = A[ind_i, :] * x[:, i]
            
        elif opts['nFlag'] == 1:
            invNu = x[:, i] / nu
            mu_invNu = mu * invNu
            Ax[ind_i, 0:1] = A[ind_i, :] * invNu - np.repmat(mu_invNu, m_i, 1)
        else:
            Ax[ind_i, 0:1] = A[ind_i, :] * x[:, i] - np.repmat(mu * x[:, i], m_i, 1)
            Ax[ind_i, 0:1] = Ax[ind_i, 0] / nu[ind_i, 0]
    if opts['init'] == 0:
        for i in xrange(n):
            x_norm = x_norm + npl.norm(np.squeeze(np.asarray(x[i, :])), q)
        if x_norm >= 1e-6:
            ratio = initFactor(x_norm, Ax, y, lamda)
            x = ratio * x
            Ax = ratio * Ax

    # The main program
    # The Armijo Goldstein line search schemes + accelerated gradient descent
    cdef long bFlag
    cdef double L, beta
    cdef np.ndarray xp
    cdef np.ndarray Axp
    cdef np.ndarray xxp
    cdef long iterStep
    cdef double alphap = 0, alpha = 1
    cdef np.ndarray s
    cdef np.ndarray As
    cdef np.ndarray ATAs
    cdef np.ndarray g
    cdef np.ndarray v
    cdef np.ndarray Av
    cdef double r_sum, l_sum
    cdef np.ndarray Axy
    cdef double norm_xxp, norm_xp
    cdef double gamma
    cdef np.ndarray t
    cdef np.ndarray tp
    cdef np.ndarray ATAx
    cdef np.ndarray ATAxp
    cdef np.ndarray s_t
    cdef np.ndarray xnew, tnew
    cdef np.ndarray Axnew
    cdef np.ndarray u, v_t
    cdef double tao
 
    if (opts['mFlag'] == 0) and (opts['lFlag'] == 0):
        bFlag = 0
        L = 1
        xp = x.copy()
        Axp = Ax.copy()
        xxp = np.zeros((n, k))
        

        funVal = np.squeeze(np.asarray(np.zeros(opts['maxIter'])))
        ValueL = np.squeeze(np.asarray(np.zeros(opts['maxIter'])))
        for iterStep in xrange(opts['maxIter']):
            # step 1
            beta = (alphap - 1) / alpha
            s = x + beta * xxp

            # step 2
            As = Ax + beta * (Ax - Axp)
            ATAs = np.zeros((n, k))
            
            for i in xrange(k):
                ind_i_m = np.arange(ind[i]+1, ind[i+1]+1)
                ind_i = ind_i_m - 1
                if opts['nFlag'] == 0:
                    tt = A[ind_i, :].T * As[ind_i, 0]
                elif opts['nFlag'] == 1:
                    tt = A[ind_i, :].T * As[ind_i, 0] - As[ind_i, 0].sum() * mu.T
                    tt = tt / nu[:, 0]
                else:
                    invNu = As[ind_i, 0] / nu[ind_i, 0]
                    tt = A[ind_i, :].T * invNu - invNu.sum() * mu.T
                ATAs[:, i] = tt
                
            g = ATAs - ATy
            xp = x.copy()
            Axp = Ax.copy()
            while 1:
                v = s - g / L
                x = eppMatrix(np.asarray(v, order='F'), n, k, lamda / L, q)
                v = x - s
                for i in xrange(k):
                    ind_i_m = np.arange(ind[i]+1, ind[i+1]+1)
                    ind_i = ind_i_m - 1
                    m_i = ind[i+1] - ind[i]
                    
                    if opts['nFlag'] == 0:
                        Ax[ind_i, 0:1] = A[ind_i, :] * x[:, i]
                    elif opts['nFlag'] == 1:
                        invNu = x[:, i] / nu
                        mu_invNu = mu * invNu
                        Ax[ind_i, 0:1] = A[ind_i, :] * invNu - np.repmat(mu_invNu, m_i, 1)
                    else:
                        Ax[ind_i, 0:1] = A[ind_i, :] * x[:, i] - np.repmat(mu * x[:, i], m_i, 1)
                        Ax[ind_i, 0:1] = Ax[ind_i, 0] / nu[ind_i, 0]
                Av = Ax - As
                r_sum = npl.norm(v, 'fro') ** 2
                l_sum = np.double(Av.T * Av)
            
                if r_sum <= 1e-20:
                    bFlag = 1
                    break
                if l_sum <= (r_sum * L):
                    break
                else:
                    L = np.maximum(2 * L, l_sum / r_sum)
            # step 3
            alphap = alpha
            alpha = (1 + np.sqrt(4 * alpha * alpha + 1)) / 2
            ValueL[iterStep] = L

            xxp = x - xp
            Axy = Ax - y
            funVal[iterStep] = np.double(Axy.T * Axy) / 2

            for i in xrange(n):
                funVal[iterStep] = funVal[iterStep] + lamda * npl.norm(np.squeeze(np.asarray(x[i, :])), q)

            if bFlag != 0:
                break

            if opts['tFlag'] == 0:
                if iterStep >= 1:
                    if np.absolute(funVal[iterStep] - funVal[iterStep-1]) <= opts['tol']:
                        break
            elif opts['tFlag'] == 1:
                if iterStep >= 1:
                    if np.absolute(funVal[iterStep] - funVal[iterStep-1]) <= (opts['tol'] * funVal[iterStep-1]):
                        break
            elif opts['tFlag'] == 2:
                if funVal[iterStep] <= opts['tol']:
                    break
            elif opts['tFlag'] == 3:
                norm_xxp = npl.norm(xxp, 'fro')
                if norm_xxp <= opts['tol']:
                    break
            elif opts['tFlag'] == 4:
                norm_xp = npl.norm(xp, 'fro')
                norm_xxp = npl.norm(xxp, 'fro')
                if norm_xxp <= (opts['tol'] * np.maximum(norm_xp, 1)):
                    break
            elif opts['tFlag'] == 5:
                if iterStep >= opts['maxIter']:
                    break

    # Reformulated problem + adaptive line search
    if (opts['mFlag'] == 1) and (opts['lFlag'] == 1) and (opts['q'] == 2):
        L = 1
        bFlag = 0
        gamma = 1
        xp = x.copy()
        Axp = Ax.copy()
        xxp = np.zeros((n, k))
        t = np.sqrt(np.power(x, 2).sum(axis=1))
        tp = t.copy()
        As = Ax.copy()

        ATAs = np.zeros((n, k))
        for i in xrange(k):
            ind_i_m = np.arange(ind[i]+1, ind[i+1]+1)
            ind_i = ind_i_m - 1
            if opts['nFlag'] == 0:
                tt = A[ind_i, :].T * As[ind_i, 0]
            elif opts['nFlag'] == 1:
                tt = A[ind_i, :].T * As[ind_i, 0] - As[ind_i, 0].sum() * mu.T
                tt = tt / nu[:, 0]
            else:
                invNu = As[ind_i, 0] / nu[ind_i, 0]
                tt = A[ind_i, :].T * invNu - invNu.sum() * mu.T
            ATAs[:, i] = tt
            
        ATAx = ATAs.copy()
        for iterStep in xrange(opts['maxIter']):
            ATAxp = ATAx.copy()
            if iterStep != 1:
                As = Ax.copy()
                for i in xrange(k):
                    ind_i_m = np.arange(ind[i]+1, ind[i+1]+1)
                    ind_i = ind_i_m - 1
                    if opts['nFlag'] == 0:
                        tt = A[ind_i, :].T * As[ind_i, 0]
                    elif opts['nFlag'] == 1:
                        tt = A[ind_i, :].T * As[ind_i, 0] - As[ind_i, 0].sum() * mu.T
                        tt = tt / nu[:, 0]
                    else:
                        invNu = As[ind_i, 0] / nu[ind_i, 0]
                        tt = A[ind_i, :].T * invNu - invNu.sum() * mu.T
                    ATAs[:, i] = tt
                ATAx = ATAs.copy()

            #--------Line Search for L begins
            Axnew = np.zeros((m, 1))
            while 1:
                if iterStep != 1:
                    alpha = (-gamma + np.sqrt(gamma ** 2 + 4 * L * gamma)) / (2 * L)
                    beta = (gamma - gamma * alphap) / (alphap * gamma + alphap * L * alpha)
                    s = x + beta * xxp
                    s_t = t + beta * (t - tp)
                    As = Ax + beta * (Ax - Axp)
                    ATAs = ATAx + beta * (ATAx - ATAxp)
                else:
                    alpha = (-1 + np.sqrt(5)) / 2
                    beta = 0
                    s = x.copy()
                    s_t = t.copy()
                    As = Ax.copy()
                    ATAs = ATAx.copy()
                g = ATAs - ATy
                u = s - g / L
                v = s_t - lamda / L
                xnew, tnew = ep21R(np.asarray(u, order='F'), np.asarray(v, order='F'), n, k)
                v = xnew - s
                v_t = tnew - s_t
                for i in xrange(k):
                    ind_i_m = np.arange(ind[i]+1, ind[i+1]+1)
                    ind_i = ind_i_m - 1
                    m_i = ind[i+1] - ind[i]
                    
                    if opts['nFlag'] == 0:
                        Axnew[ind_i, 0:1] = A[ind_i, :] * xnew[:, i]
                    elif opts['nFlag'] == 1:
                        invNu = xnew[:, i] / nu
                        mu_invNu = mu * invNu
                        Axnew[ind_i, 0:1] = A[ind_i, :] * invNu - np.repmat(mu_invNu, m_i, 1)
                    else:
                        Axnew[ind_i, 0:1] = A[ind_i, :] * xnew[:, i] - np.repmat(mu * xnew[:, i], m_i, 1)
                        Axnew[ind_i, 0:1] = Axnew[ind_i, 0] / nu[ind_i, 0]
                Av = Axnew - As
                r_sum = npl.norm(v, 'fro') ** 2 + np.double(v_t.T * v_t)
                l_sum = npl.norm(Av, 'fro') ** 2

                if r_sum <= 1e-20:
                    bFlag = 1
                    break

                if l_sum <= r_sum * L:
                    break
                else:
                    L = np.maximum(2 * L, l_sum / r_sum)
            #--------Line Search for L ends
            gamma = L * alpha * alpha
            alphap = alpha
            ValueL[iterStep] = L

            tao = L * r_sum / l_sum
            if tao >= 5:
                L = L * 0.8

            xp = x.copy()
            x = xnew.copy()
            xxp = x - xp
            Axp = Ax.copy()
            Ax = Axnew.copy()
            tp = t.copy()
            t = tnew.copy()

            Axy = Ax - y
            funVal[iterStep] = Axy.T * Axy / 2 * lamda * t.sum()

            if bFlag != 0:
                break

            if opts['tFlag'] == 0:
                if iterStep >= 1:
                    if np.absolute(funVal[iterStep] - funVal[iterStep-1]) <= opts['tol']:
                        break
            elif opts['tFlag'] == 1:
                if iterStep >= 1:
                    if np.absolute(funVal[iterStep] - funVal[iterStep-1]) <= (opts['tol'] * funVal[iterStep-1]):
                        break
            elif opts['tFlag'] == 2:
                if funVal[iterStep] <= opts['tol']:
                    break
            elif opts['tFlag'] == 3:
                norm_xxp = npl.norm(xxp, 'fro')
                if norm_xxp <= opts['tol']:
                    break
            elif opts['tFlag'] == 4:
                norm_xp = npl.norm(xp, 'fro')
                norm_xxp = npl.norm(xxp, 'fro')
                if norm_xxp <= (opts['tol'] * np.maximum(norm_xp, 1)):
                    break
            elif opts['tFlag'] == 5:
                if iterStep >= opts['maxIter']:
                    break

    # Reformulated problem + Nemirovski's line search
    if (opts['mFlag'] == 1) and (opts['lFlag'] == 0) and (opts['q'] == 2):
        L = 1
        bFlag = 0
        xp = x.copy()
        Axp = Ax.copy()
        xxp = np.zeros((n, k))

        alphap = 0
        alpha = 1
        t = np.sqrt(np.power(x, 2).sum(axis=1))
        tp = t.copy()

        for iterStep in xrange(opts['maxIter']):
            # step 1
            beta = (alphap - 1) / alpha
            s = x + beta * xxp
            s_t = t + beta * (t - tp)

            # step 2
            As = Ax + beta * (Ax - Axp)
            ATAs = np.zeros((n, k))

            for i in xrange(k):
                ind_i_m = np.arange(ind[i]+1, ind[i+1]+1)
                ind_i = ind_i_m - 1
                
                if opts['nFlag'] == 0:
                    tt = A[ind_i, :].T * As[ind_i, 0]
                elif opts['nFlag'] == 1:
                    tt = A[ind_i, :].T * As[ind_i, 0] - As[ind_i, 0].sum() * mu.T
                    tt = tt / nu[:, 0]
                else:
                    invNu = As[ind_i, 0] / nu[ind_i, 0]
                    tt = A[ind_i, :].T * invNu - invNu.sum() * mu.T
                ATAs[:, i] = tt

            g = ATAs - ATy
            xp = x.copy()
            Axp = Ax.copy()
            tp = t.copy()

            while 1:
                u = s - g / L
                v = s_t - lamda / L

                x, t = ep21R(np.asarray(u ,order='F'), np.asarray(v ,order='F'), n, k)
                v = x - s

                for i in xrange(k):
                    ind_i_m = np.arange(ind[i]+1, ind[i+1]+1)
                    ind_i = ind_i_m - 1
                    m_i = ind[i+1] - ind[i]

                    if opts['nFlag'] == 0:
                        Ax[ind_i, 0:1] = A[ind_i, :] * x[:, i]
                    elif opts['nFlag'] == 1:
                        invNu = x[:, i] / nu
                        mu_invNu = mu * invNu
                        Ax[ind_i, 0:1] = A[ind_i, :] * invNu - np.repmat(mu_invNu, m_i, 1)
                    else:
                        Ax[ind_i, 0:1] = A[ind_i, :] * x[:, i] - np.repmat(mu * x[:, i], m_i, 1)
                        Ax[ind_i, 0:1] = Ax[ind_i, 0] / nu[ind_i, 0]
            
                Av = Ax - As
                r_sum = npl.norm(v, 'fro') ** 2
                l_sum = np.double(Av.T * Av)

                if r_sum <= 1e-20:
                    bFlag = 1
                    break

                if l_sum <= r_sum * L:
                    break
                else:
                    L = np.maximum(2 * L, l_sum / r_sum)

            ValueL[iterStep] = L

            # step 3
            alphap = alpha
            alpha = (1 + np.sqrt(4 * alpha ** 2)) / 2

            xxp = x - xp
            Axy = Ax - y
            funVal[iterStep] = Axy.T * Axy / 2 + lamda * t.sum()

            if bFlag != 0:
                break

            if opts['tFlag'] == 0:
                if iterStep >= 1:
                    if np.absolute(funVal[iterStep] - funVal[iterStep-1]) <= opts['tol']:
                        break
            elif opts['tFlag'] == 1:
                if iterStep >= 1:
                    if np.absolute(funVal[iterStep] - funVal[iterStep-1]) <= (opts['tol'] * funVal[iterStep-1]):
                        break
            elif opts['tFlag'] == 2:
                if funVal[iterStep] <= opts['tol']:
                    break
            elif opts['tFlag'] == 3:
                norm_xxp = npl.norm(xxp, 'fro')
                if norm_xxp <= opts['tol']:
                    break
            elif opts['tFlag'] == 4:
                norm_xp = npl.norm(xp, 'fro')
                norm_xxp = npl.norm(xxp, 'fro')
                if norm_xxp <= (opts['tol'] * np.maximum(norm_xp, 1)):
                    break
            elif opts['tFlag'] == 5:
                if iterStep >= opts['maxIter']:
                    break
    if (opts['mFlag'] == 0) and (opts['lFlag'] == 1):
        raise ValueError('The function does not support opts.mFlag=0 & opts.lFlag=1!')

    return (np.asarray(x), np.asarray(funVal), np.asarray(ValueL))

