import numpy as np, math, sys
from thmA_scale import sieve

def run(N, theta):
    if N%2: N+=1
    mu,Lam,primes,spf = sieve(N)
    K=int(N**theta)
    Np=set(); n=N
    while n>1:
        p=int(spf[n]); Np.add(p)
        while n%p==0: n//=p
    C2=1.0
    for p in primes:
        if p>2: C2*=(1-1.0/((p-1)**2))
    S=2*C2
    for p in Np:
        if p>2: S*=(p-1)/(p-2)
    AG=1.0
    for p in primes:
        if p>2 and p not in Np: AG*=(1-1.0/((p-1)**2))
    # C = sum_{n<N} Lam(n) mu(N-n)
    idx=np.arange(1,N)
    Cval = float(np.dot(Lam[1:N], mu[N-idx].astype(np.float64)))
    # LL = sum Lam(n)Lam(N-n)
    LL = float(np.dot(Lam[1:N], Lam[N-idx]))
    # phi
    phi=np.arange(N+1,dtype=np.int64)
    for p in primes:
        phi[p::p] -= phi[p::p]//p
    Tlog=0.0; B=0.0
    for k in range(2,K):
        if mu[k]==0 or math.gcd(k,N)!=1: continue
        s=0.0
        n0=N%k
        if n0==0: n0=k
        for nn in range(n0,N,k):
            if Lam[nn]:
                s+=Lam[nn]*mu[N-nn]
        Tlog += mu[k]*math.log(k)*s
        B += mu[k]*math.log(k)/phi[k]
    E3 = Tlog - Cval*B
    return N,K,S,AG,Cval,LL,B,Tlog,E3

theta=float(sys.argv[1]) if len(sys.argv)>1 else 0.56
print(f"{'N':>8} {'C/N':>8} {'LL/N':>7} {'S':>7} {'A*G1':>7} {'B':>8} {'-S':>7} {'Tlog/N':>8} {'E3/N':>8} {'S-AG':>7}")
for N in [50000,100000,200000,400000,800000]:
    N,K,S,AG,Cval,LL,B,Tlog,E3 = run(N,theta)
    print(f"{N:>8} {Cval/N:>8.4f} {LL/N:>7.4f} {S:>7.4f} {AG:>7.4f} {B:>8.4f} {-S:>7.4f} {Tlog/N:>8.4f} {E3/N:>8.4f} {S-AG:>7.4f}")
