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
        if p not in Np: AN *= (1-1.0/(p*(p-1)))
    def lam(m):
        v=1.0; mm=m
        while mm>1:
            p=int(spf[mm])
            while mm%p==0: mm//=p
            if p not in Np: v/= (1-1.0/(p*(p-1)))
        return v
    # observed residuals
    res1=0.0; reslog=0.0; mt1=0.0; mtlog=0.0
    for m in range(1,M+1):
        if mu[m]==0: continue
        s1=0.0; sl=0.0
        for k in range(K,(N-2)//m+1):
            if mu[k]==0: continue
            if math.gcd(k,m)!=1: continue
            if math.gcd(k,N)!=1: continue
            L=Lam[N-m*k]
            if L: s1+=L; sl+=L*math.log(k)
        res1+=mu[m]*s1; reslog+=mu[m]*sl
        c = AN*lam(m)/m
        Tm = N-m*K
        mt1  += mu[m]*c*Tm
        # int_{mK}^{N} log(v/m) dv
        hi=N; lo=m*K
        I = (hi*math.log(hi/m)-hi) - (lo*math.log(lo/m)-lo)
        mtlog += mu[m]*c*I
    # G(1) constant  A(N)*G(1) = prod_{p not| N}(1-1/(p-1)^2)
    AG=1.0
    for p in primes:
        if p not in Np and p>2: AG *= (1-1.0/((p-1)**2))
        if p==2 and 2 not in Np: AG=0.0
    C2=1.0
    for p in primes:
        if p>2: C2*= (1-1.0/((p-1)**2))
    S=2*C2
    for p in Np:
        if p>2: S*= (p-1)/(p-2)
    return N,K,M,res1,mt1,reslog,mtlog,AG,S

theta=float(sys.argv[1]) if len(sys.argv)>1 else 0.56
print(f"{'N':>8} {'res1/N':>8} {'MT1/N':>8} {'reslog/N':>9} {'MTlog/N':>9} {'A*G(1)':>8} {'S(N)':>7}")
for N in [50000,100000,200000,400000]:
    N,K,M,res1,mt1,reslog,mtlog,AG,S = data(N,theta)
    print(f"{N:>8} {res1/N:>8.4f} {mt1/N:>8.4f} {reslog/N:>9.4f} {mtlog/N:>9.4f} {AG:>8.4f} {S:>7.4f}")

# limit check: Sum_{m<=x} mu(m) log m / m  ->  -1 ?
X=3000000
mu2=np.ones(X+1,dtype=np.int8)
spf2=np.zeros(X+1,dtype=np.int32)
for i in range(2,X+1):
    if spf2[i]==0:
        for j in range(i,X+1,i):
            if spf2[j]==0: spf2[j]=i
mu2=np.zeros(X+1,dtype=np.int8); mu2[1]=1
for i in range(2,X+1):
    p=int(spf2[i]); j=i//p
    mu2[i]=0 if j%p==0 else -mu2[j]
tot=0.0
marks={10**4:0,10**5:0,10**6:0,3*10**6:0}
for m in range(2,X+1):
    if mu2[m]: tot+=mu2[m]*math.log(m)/m
    if m in marks: marks[m]=tot
print("\nSum_{m<=x} mu(m) log m / m  (should -> -1):")
for k in sorted(marks): print(f"  x={k:>9}: {marks[k]:.5f}")
