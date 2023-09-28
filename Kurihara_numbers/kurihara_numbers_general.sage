#A is a 1-dimensional subspace of ModularSymbols (so defined over Q)
def ev(A,ell):
	M = A.ambient_module()
	vec = A.dual_eigenvector()
	return (M.hecke_matrix(ell)*vec)[0]

# A is a simple space of modular symbols
# w is a valuation (over some p) defined over the field of defn of form cutting out A
def good_primes(A,w,m,max,D):
	"""returns the primes ell< max such that w(ell - 1) >= m and w(a_ell(A)-2)>=m"""
	N = A.level()
	p = w.p()
	pi = w.uniformizer()
	e = 1/w(pi)
	q = 2
	ans = []
	k = A.weight()

	while q < max:
		if w(q - 1) >= m/e and N*p % q != 0:
			aq = ev(A,q) * kronecker_symbol(D,q)
#			print(q,aq,w(q-1),w(aq-1-q^(k-1)))
			if w(aq-1-q^(k-1)) >= m/e:
				ans += [q]
		q = next_prime(q)

	return ans

def In(A,w,n,D):
	p = w.p()
	pi = w.uniformizer()
	e = 1/w(pi)
	k = A.weight()

	ells = factor(n)
	ells = [t[0] for t in ells]

	m_prime = min([w(ell-1)*e for ell in ells])
	m_FC = min([w(kronecker_symbol(D,ell) * ev(A,ell)-1-ell^(k-1)) * e for ell in ells])

	m = min(m_prime,m_FC)

	return m

def log_table(ell):
	F = Integers(ell)
	g = F.multiplicative_generator()
	d = {}
	for e in range(1,ell):
		a = (g^e) % ell
		d[a] = e
	return d

# # Computes int_{a/m}^\infty z^j f(z) dz (normalized) for all 1 <= a < m and 0<=j<=r
# def precompute_integrals(A,m,r,magic=-1):
# 	k = A.weight()
# 	M = A.ambient_module()
# 	if magic == -1:
# 		magic = A.dual_eigenvector()

# 	for j in range(r+1):
# 		for a in range(m):
# 			if not (a,m,j) in ints.keys():
# 				ints[(a,m,j)] = magic.dot_product(M.modular_symbol([j,oo,-a/m]).element())

def Phi(A,j,r,magic=-1):
	if magic==-1:
		magic = A.dual_eigenvector()

	M = A.ambient_module()

	return magic.dot_product(M.modular_symbol([j,oo,r]).element())	

def lamb_slow(A,j,b,n):
	return sum([binomial(j,i) * n^i * b^(j-i) * Phi(A,i,-b/n) for i in range(j+1)])

##returns int_infty^{b/n} f(z) z^i dz
def period_integral(A,b,n,i,magic):
	# if not (b,n,i) in ints.keys():
	# 	M = A.ambient_module()
	# 	ints[(b,n,i)] = magic.dot_product(M.modular_symbol([i,oo,-b/n]).element())	
	M = A.ambient_module()

	return magic.dot_product(M.modular_symbol([i,oo,-b/n]).element())	

def lamb(A,j,b,n,magic=-1):
	if magic == -1:
		magic = A.dual_eigenvector()

	return sum([binomial(j,i) * n^i * b^(j-i) * period_integral(A,b,n,i,magic) for i in range(j+1)])

def lamb_twist(A,j,b,n,D,magic=-1):
	assert is_fundamental_discriminant(D) or D==1,"Not a good twist"

	if magic == -1:
		magic = A.dual_eigenvector()

	ans = 0
	for a in range(D):
		c = (D*b - n*a) % (n*D)
		ans += kronecker_symbol(D,a) * sum([binomial(j,i) * (n*D)^i * c^(j-i) * period_integral(A,c,n*D,i,magic) for i in range(j+1)])

	return ans

def delta(A,n,D,magic=-1):
	assert is_squarefree(n), "need square-free n"
	assert is_fundamental_discriminant(D) or D==1, "need fund disc"

	k = A.weight() 
	r = k/2 

	ells = factor(n)
	ells = [t[0] for t in ells]

	L = {}
	for ell in ells:
		L[ell] = log_table(ell)

	if magic == -1:
		magic = A.dual_eigenvector()

	R.<z> = PolynomialRing(QQ)
	ans = 0
	for a in range(n):
		if gcd(a,n)==1:
#			print(a,n)
			log_term = prod([L[ell][a%ell] for ell in ells])
			ans += lamb_twist(A,r-1,a,n,D,magic=magic) * log_term

	return ans

def compute_deltas(A,w,max_ell,depth,D,magic=-1,period_correction=0,filename=-1,vLval=""):
	if magic == -1:
		magic = A.dual_eigenvector()

	pi = w.uniformizer()
	kpi = w.residue_field()

	e = 1/w(pi)
	f = kpi.degree()
	if filename != -1:
		printwritelist(filename,[A])
		printwritelist(filename,["p =",w.p(),"valuation =",w])
		printwritelist(filename,["e =",e,"f =",f])
		if D > 1:
			printwritelist(filename,["Twisting by quadratic character of conductor",D])
	if vLval != "" and filename != -1:
		printwritelist(filename,["Valuation of central L-value:",vLval])



	qs = good_primes(A,w,depth,max_ell,D)
	if filename != -1:
		printwritelist(filename,["Good primes:",qs])	
	print("Good primes:",qs)
	# for ell in qs:
	# 	#print("Working on",ell)
	# 	dn = delta(A,ell,D,magic=magic) * pi^period_correction
	# 	vdn = w(dn)
	# 	m = In(A,w,D)
	# 	if vdn > m:
	# 		print(ell,": vanish")
	# 	else:
	# 		print(ell,": non-vanishing")

	if len(qs) < 2:
		print("Not enough good primes to compute anything")
		if filename != -1:
			printwritelist(filename,["Not enough good primes to compute anything"])
	van = 0 
	nvan = 0
	for ell1 in qs:
		for ell2 in qs:
			if ell1 < ell2:
				#print("Working on",(ell1,ell2))
				dn = delta(A,ell1*ell2,D,magic=magic) * pi^period_correction
				vdn = w(dn)
				m = In(A,w,ell1*ell2,D)
				if vdn > m:
					print((ell1,ell2),": vanish")
					if filename != -1:
						printwritelist(filename,[(ell1,ell2),": vanish"])
					van += 1
				else:
					print((ell1,ell2),": non-vanishing")
					if filename != -1:
						printwritelist(filename,[(ell1,ell2),": non-vanish"])
					nvan += 1

	if van + nvan > 0:
		print("Total vanishing:non-vanishing",van,":",nvan,"=",round(nvan/(van+nvan)*100.0,4),"% non-vanishing")				
		if filename != -1:
			printwritelist(filename,["Total vanishing:non-vanishing,",van,":",nvan,"=",round(nvan/(van+nvan)*100.0,4),"% non-vanishing"])

	if filename != -1:
		printwritelist(filename,[""])

	return "Done"

def form_deltas_in_fixed_weight_and_level(N,k,ps,max_ell,depth,Ds,skip_odd_rank=true,skip_unit_Lval=true,skip_Eisen=true,filename=-1):
	if isinstance(Ds,sage.rings.integer.Integer):
		Ds = [Ds]
	print("Working on level",N,"and weight",k)
	sign = (-1)^((k+2)/2)
	M = ModularSymbols(N,k,sign).cuspidal_subspace()
	As = M.decomposition()
	print("-there are",len(As),"Galois conjugacy classes of forms")
	for A in As:
		print("Working on GC class",As.index(A)+1,"out of",len(As))
		print(A)
		magic = A.dual_eigenvector()
		K = magic.parent().base_ring()
		phi = form_modsym_from_decomposition(A)
		for D in Ds:
			print("    Twisting by",D)
			Lval = lamb_twist(A,(k-2)/2,0,1,D,magic=magic)
			if Lval !=0 or not skip_odd_rank:
				for p in ps:
					if N % p != 0 and D % p != 0:
						print("    Taking p =",p)
						v = QQ.valuation(p)
						ws = v.extensions(K)
						print("       -there are",len(ws),"prime(s) over",p)
						print()
						for w in ws:
							pi = w.uniformizer()
							e = 1/w(pi)
							print("       --Working with valuation",ws.index(w)+1,"/",len(ws),"with e =",e,"and f =",w.residue_field().degree())
							eis = eisenstein(A,w)
							print("       --Is this case Eisenstein?:",eis)
							if not skip_Eisen or not eis:
	#							print("       --Valuation of modular symbol:",phi.valuation(w,remove_binom=true))	
								period_correction = -phi.valuation(w,remove_binom=true)*e 	
								print("       --The central L-value has valuation:",(w(Lval)-phi.valuation(w,remove_binom=true))*e)
								if Lval != 0 or not skip_odd_rank:
									if w(Lval) - phi.valuation(w,remove_binom=true) > 0 or not skip_unit_Lval:
										compute_deltas(A,w,max_ell,depth,D,magic=magic,period_correction=period_correction,vLval=(w(Lval)-phi.valuation(w,remove_binom=true))*e,filename=filename)
									else:
										print("       --Skipping because L-value is a unit")
							else:
								print("       --Skipping because it is Eisenstein")
							print()
			else:
				print("       --Skipping because L-value is 0")
			print()	

	del M

def eisenstein(A,w,max_check=50):
	N = A.level()
	k = A.weight()
	p = w.p()

	bool = true
	q = 2
	while q < max_check and bool:
		if N*p % q != 0:
			bool = bool and w(ev(A,q)-1-q^(k-1)) > 0
		q = next_prime(q)

	return bool

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



def clearfile(filename):
	f = open(filename, 'w')
	f.close()

def writethingtofile(filename, thing, pre = "", post=""):
	f = open(filename, 'a')
	f.write(pre+str(thing)+post)
	f.close()

def writenewlinetofile(filename): 
	f = open(filename, 'a')
	f.write("\n")
	f.close()

def printwritelist(filename, listy, separator = " "):
	#print(liststring(listy))
	if filename != "": 
		writelisttofile(filename, listy, separator)

def printwritestring(filename, stringy):
	print( stringy)
	writelisttofile(filename, [stringy])


def printliststring(listy, separator = " "):
	s = ""
	for thing in listy: 
		s = s + str(thing) + separator
	print(s)


def writelisttofile(filename, listy, separator = " "):
#
# copies entries of a list into a file and moves to the next line. 
#
	f = open(filename, 'a')
	for thing in listy:
		f.write(str(thing)+separator)
	f.write("\n")
	f.close()



def makefilename(fileprefix, mod, accuracy):
#
# filename will look like "deltapowersmod3at10mil"
#
	return fileprefix+"mod"+str(mod)+"at"+makeaccuracystring(accuracy)