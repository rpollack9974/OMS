#A is a 1-dimensional subspace of ModularSymbols (so defined over Q)
def ev(A,ell):
	S = A.q_eigenform(ell+1)
	return S[ell]

def next_good_prime(A,p,m,q=-1):
	q = next_prime(q)
	aq = ev(A,q)
	while q % (p^m) != 1 or (aq - 2) % (p^m) != 0 or A.level()%q == 0:
		q = next_prime(q)
		aq = ev(A,q)

	return q

def log_table(ell):
	F = Integers(ell)
	g = F.multiplicative_generator()
	d = {}
	for e in range(1,ell):
		a = (g^e) % ell
		d[a] = e
	return d

def precompute(A,m):
	"""A is a piece of the result from a command like ModularSymbols(---).decomposition()"""
	M=A.ambient_module()
	N=A.level()
	k=A.weight()
	manin=manin_relations(N)
	chi=M.character()
	if chi.order()==1:
		chi = None

	w=A.dual_eigenvector()
	K=w.parent().base_field()
	v=[]
	R.<X,Y>=PolynomialRing(K,2)
	for s in range(0,len(manin.gens)):
		rs=manin.gens[s]
		g=manin.mats[rs]
		a=g[0,0]
		b=g[0,1]
		c=g[1,0]
		d=g[1,1]
		ans=0
		if c!=0:
			r1=a/c
		else:
			r1=oo
		if d!=0:
			r2=b/d
		else:
			r2=oo
		for j in range(k-1):
			coef=w.dot_product(M.modular_symbol([j,r1,r2]).element())
			ans=ans+X^j*Y^(k-2-j)*binomial(k-2,j)*coef
		v=v+[symk(k-2,poly=ans,base_ring=K,chi=chi)]
	return modsym_symk(N,v,manin)


def delta(phi,n,r,p,A,twist):
	assert is_squarefree(n), "need square-free n"
#	assert is_fundamental_discriminant(twist) or twist==1, "need fund disc"

	ells = factor(n)
	ells = [t[0] for t in ells]

	m_prime = min([(ell-1).valuation(p) for ell in ells])
	m_FC = min([(kronecker_symbol(twist,ell)*ev(A,ell)-2).valuation(p) for ell in ells])

	m = min(m_prime,m_FC)

	L = {}
	for ell in ells:
		L[ell] = log_table(ell)

	R.<z> = PolynomialRing(QQ)
	ans = 0
	for a in range(n):
		if gcd(a,n)==1:
			#print(a,n)
			ans += phi.lamb_twist(r-1,a,n,twist) * prod([L[ell][a%ell] for ell in ells])

	return ans % (p^m)

def form_symbol():
	M = ModularSymbols(5,4).cuspidal_subspace()
	Mp = M.plus_submodule()
	Mm = M.minus_submodule()
	phip = form_modsym_from_decomposition(Mp)
	phim = form_modsym_from_decomposition(Mm)
	phi = phip + phim
	return phi,Mm

def do_it(M,p,twist):
	Mp = M.plus_submodule()
	Mm = M.minus_submodule()
	phip = form_modsym_from_decomposition(Mp)
	phim = form_modsym_from_decomposition(Mm)
	phi = phip + phim
	k = phi.weight()+2
	r = k/2
	print("delta_1 has valuation:",phi.lamb_twist(r-1,0,1,twist).valuation(p))
	ell1 = next_good_prime(Mp,p,1)
	print("Using the prime",ell1)
	print("delta_",ell1,":",delta(phi,ell1,k/2,p,Mp,twist))
	ell2 = next_good_prime(Mp,p,1,q=ell1)
	n = ell1 * ell2
	print("Using the primes",ell1,ell2)
	print("delta_",n,":",delta(phi,n,k/2,p,Mp,twist))

def compute_deltas(M,Q):
	phi = form_modsym_from_decomposition(M)
	qs = []
	q = -1
	while q < Q:
		q = next_good_prime(M,3,1,q=q)
		qs += [q]
	print("Good primes:",qs)
	for ell in qs:
		print("Working on",ell)
		d = delta(phi,ell,2,3,M,61)
		print(ell1,":",d)

	for ell1 in qs:
		for ell2 in qs:
			if ell1 != ell2:
				print("Working on",(ell1,ell2))
				d = delta(phi,ell1*ell2,2,3,M,61)
				print((ell1,ell2),d)
	return "Done"
