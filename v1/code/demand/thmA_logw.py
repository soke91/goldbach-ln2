import numpy as np, math, sys
from thmA_scale import sieve

def run(N, theta):
    if N%2: N+=1
    mu,Lam,primes,spf = sieve(N)
    K=int(N**theta); M=(N-2)//K; t=N-1
    Np=set(); n=N
    while n>1:
        p=int(spf[n]); Np.add(p)
        while n%p==0: n//=p
    # singular series S(N) = 2 C2 prod_{p|N,p>2}(p-1)/(p-2)
    C2=1.0
    for p in primes:
        if p>2: C2 *= (1-1.0/((p-1)**2))
    S=2*C2
    for p in Np:
        if p>2: S *= (p-1)/(p-2)
    res1=0.0; reslog=0.0
    for m in range(1,M+1):
        if mu[m]==0: continue
        s1=0.0; sl=0.0
        kmax=(N-2)//m
        for k in range(K,kmax+1):
            if mu[k]==0: continue
            if math.gcd(k,m)!=1: continue
            if math.gcd(k,N)!=1: continue
            L=Lam[N-m*k]
            if L!=0.0:
                s1+=L; sl+=L*math.log(k)
        res1+=mu[m]*s1; reslog+=mu[m]*sl
    # complete piece for log-weight: sum_u Lam(N-u) mu(u) * (-Lambda(a(u)))
    complog=0.0
    for u in range(1,N-1):
        if mu[u]==0: continue
        # a(u) = part of u coprime to N
        a=u
        for p in Np:
            while a%p==0: a//=p
        if a==1: continue
        La=Lam[a]
        if La!=0.0:
            complog += Lam[N-u]*mu[u]*(-La)
    return N,K,M,res1,reslog,complog,S

theta=float(sys.argv[1]) if len(sys.argv)>1 else 0.56
print(f"theta'={theta}   (S = singular series)")
print(f"{'N':>8} {'res(w=1)/N':>11} {'res(w=log)/N':>13} {'complete_log/N':>15} {'S':>7} {'S*N/N':>7}")
for N in [50000,100000,200000,400000]:
    N,K,M,res1,reslog,complog,S = run(N,theta)
    print(f"{N:>8} {res1/N:>11.5f} {reslog/N:>13.5f} {complog/N:>15.5f} {S:>7.4f} {S:>7.4f}")
