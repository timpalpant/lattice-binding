import math, numpy

def heaviside(x):
    return x >= 0

def factorial(N):
    if not isinstance(N, int):
        N[N<=0] = 1
        
    if numpy.all(N<=1):
        if isinstance(N,int):
            return 1
        else:
            return N
    else:
        b = N-1
        if not isinstance(b,int):
            b[b<=0] = 1
        return N * factorial(b)

def choose(n, k):
    return factorial(n) / (factorial(k) * factorial(n-k))

def Q(L, N):
    '''
    Canonical partition function for N rods of unit length 
    on a line segment of length L
    '''
    return heaviside(L-N) * (L-N)**N / numpy.array(factorial(N), dtype=float)

def Z(L, u, beta=1):
    '''
    Grand canonical partition function for rods of unit length on 
    a line segment of length L in a bath with chemical potential u
    '''
    if isinstance(L,int) or isinstance(L,float):
        N = numpy.arange(Nmax(L)+1)
        return numpy.sum(numpy.exp(beta*N*u) * Q(L,N))
    else:
        z = numpy.zeros_like(L, dtype=float)
        for i,l in enumerate(L):
            z[i] = Z(l,u,beta)
        return z

def Pn(L, u, N, beta=1):
    '''
    Probability of finding N particles on a line segment of length L
    in a bath with chemical potential u
    '''
    return numpy.exp(beta*N*u) * Q(L,N) / Z(L,u,beta)

def Nmax(L):
    '''
    The maximum number of particles that will fit on a line
    segment of length L
    '''
    return numpy.floor(L)

def Nmean(L, u, beta=1):
    '''
    Mean number of particles on a line segment of length L
    in a bath with chemical potential u
    '''
    if isinstance(L,int) or isinstance(L,float):
        N = numpy.arange(Nmax(L)+1)
        p = Pn(L,u,N,beta)
        return numpy.sum(N*p)
    else:
        nm = numpy.zeros_like(L, dtype=float)
        for i,l in enumerate(L):
            nm[i] = Nmean(l,u,beta)
        return nm

def density(L, u ,beta=1):
    return Nmean(L,u,beta) / L

def r1(L, u, x, beta=1):
    '''
    Probability of finding a particle at x (0.5 <= x <= L-0.5)
    on a line segment of length L at chemical potential u
    (one particle distribution function)
    '''
    return numpy.exp(beta*u) * Z(x-0.5,u,beta) * Z(L-x-0.5,u,beta) / Z(L,u,beta)

def r2(L, u, x1, x2, beta=1):
    '''
    Probability of finding a particle at x2 (x1 <= x2 <= L-0.5),
    given that there is a particle at x1, on a line segment of length L
    at chemical potential u (two particle distribution function)
    '''
    return numpy.exp(2*beta*u) * Z(x1-0.5,u,beta) * Z(x2-x1-1,u,beta) * Z(L-x2-0.5,u,beta) / Z(L,u,beta)

def Prn(L, N, r):
    '''
    Probability of finding the nearest neighbor of a particle
    within r, on a line segment of length L with N particles of unit width
    '''
    rho = N / L
    return 1 - numpy.exp(2*rho*(r-1)/(rho-1))

def rmean(L, N):
    '''
    Mean distance to the nearest neighbor if there are N particles
    on a line segment of length L
    '''
    return (1+L/N) / 2

def Pr(L, u, r, beta=1):
    '''
    Probability of finding the nearest neighbor of a particle
    within r, on a line segment of length L with chemical potential u
    '''
    if isinstance(L,int) or isinstance(L,float):
        N = numpy.arange(Nmax(L)+1)
        return numpy.sum(Prn(L,N,r)*Pn(L,u,N,beta))
    else:
        pr = numpy.zeros_like(L, dtype=float)
        for i,l in enumerate(L):
            pr[i] = Pr(l,u,r,beta)
        return pr

def hn(L, N, r, beta=1):
    '''
    Probability of finding a hole of size r
    i.e. the probability that there is not a particle within r+1
    for N particles on line segment of length L
    '''
    return 1 - Prn(L,N,r+1)

def h(L, u, r, beta=1):
    '''
    Probability of finding a hole of size r
    i.e. the probability that there is not a particle within r+1
    for a line segment of length L with chemical potential u
    '''
    return 1 - Pr(L,u,r+1,beta)

def entropy(v, base=2):
    '''
    Entropy of a vector
    '''
    p = v / numpy.sum(v)
    return -numpy.sum(p*numpy.log(p)) / numpy.log(base)

def bistability(L, u, beta=1):
    '''
    Bistability of an array of length L at chemical potential u.
    The entropy of the probability distribution for number of particles.
    '''
    if isinstance(L,int) or isinstance(L,float):
        N = numpy.arange(Nmax(L)+1)
        return entropy(Pn(L,u,N,beta))
    else:
        B = numpy.zeros_like(L)
        for i,l in enumerate(L):
            B[i] = bistability(l,u,beta)
        return B
        
# H. Davis, 1990
def omega(L, y, u, beta=1):
    N = numpy.arange(int(N))
    deBroglie = h / math.sqrt(3*m*beta)
    return numpy.sum((numpy.exp(N*beta*u)/(factorial(N)*deBroglie**N)) * (y-N)**N)

# T. Chou, 2006
def Q_unwrap(L, N, eps, w):
    '''
    Canonical partition function for N rods of unit length on a 
    line segment of length L, where each unit of adsorbed
    particle contributes eps enthalpy.
    See: Chou, T. (2009) An exact theory of histone-DNA adsorption and wrapping. 
         Europhysics Letters 62(5): 753-759.
    '''
    kstar = min(int((L-N)/float(w)), N)
    total = 0
    for p in xrange(N):
        for k in xrange(kstar+1):
            ak = (-1)**k * choose(N, k)
            bp = factorial(2*(N-1) - p) / (factorial(p) * factorial(N-1-p))
            total += float(ak * bp) * (L - N - k*w)**p * eps**p \
                     * (numpy.exp(-eps*(L-N)) - (-1)**p * numpy.exp(-k*eps*w))
    return ((-1)**N * N / (factorial(N) * eps**(2*N-1))) * total

def Z_unwrap(L, u, eps, w, beta=1):
    fug = numpy.exp(u - eps)
    return sum(fug**N * Q_unwrap(L, N, eps, w) for N in xrange(1,int(L)))
    
def r1_unwrap(L, u, x, eps, w, beta=1):
    return Z_unwrap(x, u, eps, w) * Z_unwrap(L-x, u, eps, w) \
        / Z_unwrap(L, u, eps, w)

##
# Chereji et al. 2013
##

from scipy.misc import logsumexp
def n(L, u, a_min, a_max, beta=1):
    F = numpy.zeros(L+1)
    for i in xrange(1, L+1):
        j = numpy.arange(max(i-a_max,0), max(i-a_min,0)+1)
        s = numpy.zeros(len(j)+1)
        s[1:] = F[j] - F[i-1] + beta*u(j,i-1)
        F[i] = F[i-1] + logsumexp(s)
    R = numpy.zeros(L+1)
    for i in xrange(L-1, -1, -1):
        j = numpy.arange(min(i+a_min,L), min(i+a_max,L)+1)
        s = numpy.zeros(len(j)+1)
        s[1:] = R[j] - R[i+1] + beta*u(i+1,j)
        R[i] = R[i+1] + logsumexp(s)
    #assert numpy.allclose(F[L], R[0]), \
    #   "Z = F[L+1] = %f != R[0] = %f" % (F[L], R[0])
    lnZ = R[0]
    n = numpy.zeros((L, a_max+1))
    for i in xrange(L):
        for a in xrange(a_min, min(a_max+1, L-i)):
            n[i,a] = numpy.exp(F[i+1] + R[i+a+1] - lnZ + beta*u(i,i+a+1))
    return n

def dyads(n):
    L, a_max = n.shape
    d = numpy.zeros(L)
    half_a = numpy.arange(a_max) / 2
    for i in xrange(L-half_a.max()):
        d[i+half_a] += n[i,:]
    for i in xrange(L-half_a.max(), L):
        for a in xrange(2*(L-i)-1):
            d[i+a/2] += n[i,a]
    return d

def occ(n):
    L, a_max = n.shape
    o = numpy.zeros(L)
    for a in xrange(1, a_max):
        for i in xrange(L-a):
            o[i:i+a] += n[i,a]
    return o

def u(i, j):
    if isinstance(i, int):
        if 1000 < i < 1300:
            return 15*numpy.exp(-(i-1150)**2/60**2)
        else:
            return -1 + 0.01*numpy.abs(i-j)
    else:
        v = numpy.zeros(i.shape)
        barrier = numpy.logical_and(i > 1000, i < 1300)
        v[barrier] = 15*numpy.exp(-(i-1150)**2/60**2)
        v[numpy.logical_not(barrier)] = -1 + 0.01*numpy.abs(i-j)
        return v