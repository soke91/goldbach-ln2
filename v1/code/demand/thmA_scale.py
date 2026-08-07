import numpy as np, math, sys

def sieve(N):
    spf = np.zeros(N+1, dtype=np.int32); primes=[]
    for i in range(2, N+1):
        if spf[i]==0:
            spf[i]=i; primes.append(i)
            for j in range(i*i, N+1, i):
                if spf[j]==0: spf[j]=i
    mu = np.zeros(N+1, dtype=np.int8); mu[1]=1
    for i in range(2, N+1):
        p=int(spf[i]); j=i//p
        mu[i] = 0 if j%p==0 else -mu[j]
    Lam = np.zeros(N+1)
    for p in primes:
        q=p
        while q<=N: Lam[q]=math.log(p); q*=p
    return mu, Lam, primes, spf

def run(N, theta):
    if N%2: N+=1
    mu,Lam,primes,spf = sieve(N)
    K=int(N**theta); M=(N-2)//K; t=N-1
    Np=set(); n=N
    while n>1:
        p=int(spf[n]); Np.add(p)
        while n%p==0: n//=p
    AN=1.0
    for p in primes:
        if p not in Np: AN *= (1-1.0/(p*(p-1)))
    res=0.0; absres=0.0; MT=0.0
    for m in range(1,M+1):
        if mu[m]==0: continue
        s=0.0
        kmax=(N-2)//m
        # k >= K, gcd(k,m)=1, mu(k)!=0, gcd(k,N)=1
        for k in range(K,kmax+1):
            if mu[k]==0: continue
            if math.gcd(k,m)!=1: continue
            if math.gcd(k,N)!=1: continue
            s += Lam[N-m*k]
        lm=1.0; mm=m
        while mm>1:
            p=int(spf[mm])
            while mm%p==0: mm//=p
            if p not in Np: lm /= (1-1.0/(p*(p-1)))
        Tm=min(t,N-m*K)
        mt = AN*lm*Tm/m if Tm>2 else 0.0
        MT += mu[m]*mt
        res += mu[m]*s
        absres += abs(s-mt)
    return N,K,M,res,MT,absres,AN

theta=float(sys.argv[1]) if len(sys.argv)>1 else 0.56
print(f"theta'={theta}")
print(f"{'N':>8} {'K':>6} {'M':>5} {'residual':>12} {'MT':>10} {'res/N':>9} {'res/N^(1-th/2)':>15} {'Sum|err_m|':>12}")
for N in [25000,50000,100000,200000,400000,800000]:
    N,K,M,res,MT,absres,AN = run(N,theta)
    e=1-theta/2
    print(f"{N:>8} {K:>6} {M:>5} {res:>12.1f} {MT:>10.1f} {res/N:>9.5f} {res/N**e:>15.4f} {absres:>12.1f}")
