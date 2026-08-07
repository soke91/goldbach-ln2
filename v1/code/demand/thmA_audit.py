import numpy as np, math, sys

def sieve(N):
    mu = np.ones(N+1, dtype=np.int8)
    primes=[]
    spf = np.zeros(N+1, dtype=np.int64)
    for i in range(2, N+1):
        if spf[i]==0:
            spf[i]=i; primes.append(i)
        for p in primes:
            if p>spf[i] or i*p>N: break
            spf[i*p]=p
    # mobius via spf
    mu = np.zeros(N+1, dtype=np.int8); mu[1]=1
    for i in range(2, N+1):
        p = spf[i]; j = i//p
        if j % p == 0: mu[i]=0
        else: mu[i] = -mu[j]
    # Lambda
    Lam = np.zeros(N+1)
    for p in primes:
        q=p
        while q<=N:
            Lam[q]=math.log(p); q*=p
    return mu, Lam, primes, spf

N = 400000
if len(sys.argv)>1: N=int(sys.argv[1])
if N%2: N+=1
mu, Lam, primes, spf = sieve(N)

theta = 0.56
K = int(N**theta)
print("N=",N,"K=",K,"M=",N//K)

# primes dividing N
Np=[]; n=N
while n>1:
    p=spf[n]; Np.append(p)
    while n%p==0: n//=p
Np=set(Np)

# --- direct computation of T(t) for t = N-1 ---
t = N-1
direct = 0.0
for k in range(1, K):
    if mu[k]==0: continue
    if math.gcd(k,N)!=1: continue
    s=0.0
    n0 = N % k
    if n0==0: n0=k
    for n in range(n0, t+1, k):
        if n>=2 and Lam[n]!=0.0 and mu[N-n]!=0:
            s += Lam[n]*mu[N-n]
    direct += mu[k]*s
print("direct T(t)      =", direct)

# --- complete piece: sum over u with rad(u)|N ---
comp = 0.0
# squarefree divisors of N
divs=[1]
for p in Np:
    divs += [d*p for d in divs]
for u in divs:
    if 1 <= u <= N-2 and u >= N-t:
        comp += Lam[N-u]*mu[u]
print("complete piece   =", comp)

# --- residual: m <= (N-2)/K, k >= K ---
M = (N-2)//K
res = 0.0
for m in range(1, M+1):
    if mu[m]==0: continue
    s=0.0
    kmax = (N-2)//m
    for k in range(K, kmax+1):
        if mu[k]==0: continue          # mu^2(k)=1
        if math.gcd(k,m)!=1: continue
        u = m*k
        if u < N-t: continue
        s += Lam[N-u]
    res += mu[m]*s
print("residual piece   =", res)
print("complete-residual=", comp-res, "  (should equal direct)")

# --- main term prediction for residual ---
# c(m) = A(N) * lambda(m)/m ,  A(N)=prod_{p not| N}(1-1/(p(p-1)))
AN = 1.0
for p in primes:
    if p not in Np:
        AN *= (1 - 1.0/(p*(p-1)))
# tail of the Euler product beyond sieve range is negligible
def lam(m):
    v=1.0; mm=m
    while mm>1:
        p=spf[mm]
        while mm%p==0: mm//=p
        if p not in Np: v /= (1 - 1.0/(p*(p-1)))
    return v
MT=0.0
for m in range(1,M+1):
    if mu[m]==0: continue
    Tm = min(t, N-m*K)
    if Tm<2: continue
    MT += mu[m]*lam(m)*Tm/m
MT *= AN
print("A(N)=",AN)
print("predicted MT     =", MT)
print("residual/N =", res/N, "  MT/N =", MT/N, "  direct/N =", direct/N)

# m=1 term alone, to show the scale that must cancel
m1 = 0.0
for k in range(K,(N-2)+1):
    if mu[k]!=0: m1 += Lam[N-k]
print("m=1 term alone   =", m1, " predicted A(N)*T_1 =", AN*min(t,N-K))
