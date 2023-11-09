#return true if z is a p-th power modulo N (if z = -1 then z = p)
def pthpowermodN(p,N,z=-1):
	assert p.is_prime and N.is_prime(),"both variables are prime"
	if z == -1:
		z = p
	R=PolynomialRing(Integers(N),'x')
	x = R.gens()[0]
	return len((x^p-z).roots())>0

def discrete_log(x,N,p,g):
	x = x%N
	ans = 0
	while (g^ans) % N != x:
		ans += 1
#		print(g,ans,(g^ans)%N,x,g.parent(),x.parent())
	return ans % p 

def goren_cond(p,N,k):
	assert p.is_prime and N.is_prime(),"both variables are prime"
	assert N%p==1,"N = 1 (mod p)"
	ZN = Integers(N)
	g = ZN.multiplicative_generator()
	z = g^((N-1)/p)

	return sum([i^(k-2) * discrete_log(1-z^i,N,p,g) for i in range(1,p)])%p!=0

def t(i,k,N):
	return sum([(j^(k-1)) % (N-1) for j in range(1,i)])

def rank1(p,N,k):
	assert p.is_prime and N.is_prime(),"both variables are prime"
	assert N%p==1,"N = 1 (mod p)"
	z=prod([(i^t(i,k,N)) % N for i in range(1,N)]) % N
	return not pthpowermodN(p,N,z)

def rank1_merel(p,N):
	assert p.is_prime and N.is_prime(),"both variables are prime"
	assert N%p==1,"N = 1 (mod p)"
	z = prod([i^i % N for i in range(1,(N+1)/2)]) % N 
	return not pthpowermodN(p,N,z)


def regular(p,k):
	return (k % (p-1) != 0) and (bernoulli(k) % p != 0)

def eisen_submodule(M):
	N = M.level()
	p = M.base_ring().characteristic()
	assert p > 0,"Need positive characteristic"
	k = M.weight()
	for q in primes(100):
		if (N*p % q) != 0:			
			Tq = M.hecke_operator(q)
			M = ((Tq-(1+q^(k-1)))^1000).kernel()
	return M 

def analyze(p,N,k):
	gor = goren_cond(p,N,k)
	r1 = rank1(p,N,k)
	ppN = not pthpowermodN(p,N)
	reg = regular(p,k)
	print("Trying (N,p,k) =",N,p,k)
	if not gor:
		print("THIS IS A NON-GORENSTEIN EXAMPLE!!!")
	print("Goren :",gor)
	if not rank1:
		print("THIS IS A RANK > 1 EXAMPLE!!!")		
	print("Rank 1 :",r1)
	print("U_p-1 gen :",ppN)
	print("Regular :",reg)
	M = ModularSymbols(N,k,1).cuspidal_subspace().new_subspace()
	As = M.decomposition()
	eisen_phis = []
	for A in As:
		phi = form_modsym_from_decomposition(A)
		phis = phi.coerce_to_Qp(p,10)
		for t in phis:
			phiQp = t[0]
			q = 2
			eisen = true
			while eisen and q<30:
				if N*p % q != 0:
					eisen = eisen and (phiQp.is_Tq_eigen(q)[1] - (1+q^(k-1))) % p == 0
				q = next_prime(q)
			if eisen:
				eisen_phis += [t]
	if len(eisen_phis)==0:
		print("FAILED --- no eisenstein symbols found")
#	print("Found",len(eisen_phis),"residually eisenstein modular symbol(s)")
	for t in eisen_phis:
		phi = t[0]
		ap = phi.is_Tq_eigen(p)[1]
		phip = phi.p_stabilize_ordinary(p,ap,10)
		alphap = phip.is_Tq_eigen(p)[1]
		vap = (alphap-1).valuation(p)
		print("v_p(a_p-1) :",vap)
		Phi = phip.lift_to_OMS_eigen(p,5,ap=alphap,verbose=false)
		assert Phi.is_Tq_eigen(2)[1]==Infinity,"not an eigensymbol"
		mus = []
		lams = []
		passed = true
		for j in range(0,p-1,2):
			v = Phi.valuation()
			assert v>=0,"why is the valuation negative??"
			Phi = Phi.scale(p^(-v))
			Phi = Phi.change_precision(Phi.num_moments()-v)
			done = false
			coefs = []
			n=0
			while not done and n<Phi.num_moments():
				c = pLfunction_coef(Phi,alphap,n,j,1,1+p)
				coefs += [c]
				if ppN and r1:
					done = c.valuation(p) == vap
				else:
					done = c.valuation(p) == 1
				n += 1
			if not done:
				print("FAILED!!!!!!!!!!!!!!!!!!!!!!!!! for j =",j)
			vs = [a.valuation(p) for a in coefs]
			mu = min(vs)
			lam = vs.index(mu)
			mus += [mu]
			lams += [lam]
			passed = passed and done
		print("mus =",mus)
		print("lambdas =",lams)
		if passed:
			print("Passed")
		else:
			print("FAILED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

def run_loop(maxN,minN=-1):
	if minN == -1:
		minN = 2
	for N in primes(minN,maxN):
		for p in prime_factors(N-1):
			if p > 3:
				for k in range(2,p-2,2):
 					analyze(p,N,k)
 					print()
 					print()



