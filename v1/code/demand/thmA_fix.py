import numpy as np, math, sys
from thmA_scale import sieve

def data(N, theta):
    if N%2: N+=1
    mu,Lam,primes,spf = sieve(N)
    K=int(N**theta); M=(N-2)//K
    Np=set(); n=N
    while n>1:
        p=int(spf[n]); Np.add(p)
        while n%p==0: n//=p
    AN=1.0
    for p in primes:
        if p not in Np: AN*=(1-1.0/(p*(p-1)))
    C2=1.0
    for p in primes:
        if p>2: C2*=(1-1.0/((p-1)**2))
    S=2*C2
    for p in Np:
        if p>2: S*=(p-1)/(p-2)
    def lam(m):
        v=1.0; mm=m
        while mm>1:
            p=int(spf[mm])
            while mm%p==0: mm//=p
            if p not in Np: v/=(1-1.0/(p*(p-1)))
        return v
    res1=reslog=mt1=mtlog=0.0
    for m in range(1,M+1):
        if mu[m]==0: continue
        s1=sl=0.0
        for k in range(K,(N-2)//m+1):
            if mu[k]==0: continue
            if math.gcd(k,m)!=1 or math.gcd(k,N)!=1: continue
            L=Lam[N-m*k]
            if L: s1+=L; sl+=L*math.log(k)
        res1+=mu[m]*s1; reslog+=mu[m]*sl
        if math.gcd(m,N)!=1:      # <-- degenerate class: NO main term
            continue
        c=AN*lam(m)/m; Tm=N-m*K
        mt1 += mu[m]*c*Tm
        lo=m*K; hi=N
        I=(hi*math.log(hi/m)-hi)-(lo*math.log(lo/m)-lo)
        mtlog += mu[m]*c*I
    return N,K,M,res1,mt1,reslog,mtlog,S

theta=float(sys.argv[1]) if len(sys.argv)>1 else 0.56
print("corrected main terms: c(m)=0 unless (m,N)=1")
print(f"{'N':>8} {'res1/N':>8} {'MT1/N':>8} {'reslog/N':>9} {'MTlog/N':>9} {'S(N)':>7} {'(rl-MTl)/N':>11}")
for N in [50000,100000,200000,400000,800000]:
    N,K,M,res1,mt1,reslog,mtlog,S = data(N,theta)
    print(f"{N:>8} {res1/N:>8.4f} {mt1/N:>8.4f} {reslog/N:>9.4f} {mtlog/N:>9.4f} {S:>7.4f} {(reslog-mtlog)/N:>11.4f}")
