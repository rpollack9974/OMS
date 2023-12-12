load("master.sage")
load("Kurihara_numbers/reductions.sage")
load("Kurihara_numbers/LMFDB.sage")

#A is a 1-dimensional subspace of ModularSymbols (so defined over Q)
def ev(A,ell):
	M = A.ambient_module()
	vec = A.dual_eigenvector()
	v = M.hecke_matrix(ell)*vec

	a = 0
	while vec[a] == 0:
		a += 1
	return v[a]/vec[a]

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

def next_prime_1_mod_pr(q,p,r):
	q = next_prime(q)
	while (q-1).valuation(p) < r:
		q = next_prime(q)

	return q

def next_good_prime(q,A,w,m,D,max):
	"""returns the next good prime ell< max such that w(ell - 1) >= m and w(a_ell(A)-2)>=m"""
	N = A.level()
	p = w.p()
	pi = w.uniformizer()
	e = 1/w(pi)
	ans = []
	k = A.weight()
	q = next_prime_1_mod_pr(q,p,m/e)

	aq = ev(A,q) * kronecker_symbol(D,q)
	while (w(aq-1-q^(k-1)) < m/e or N*p % q == 0) and q <= max:
		q = next_prime_1_mod_pr(q,p,m/e)
		aq = ev(A,q) * kronecker_symbol(D,q)

	if q > max:
		return -1
	else:
		return q

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
	assert m>0,"Not a good prime being used!"

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
				if vdn >= m:
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

def find_nonzero_delta(A,w,max_ell,depth,D,magic=-1,period_correction=0,filename=-1,vLval=""):
	if magic == -1:
		magic = A.dual_eigenvector()

	pi = w.uniformizer()
	kpi = w.residue_field()

	e = 1/w(pi)
	f = kpi.degree()
#	if filename != -1:
#		printwritelist(filename,[A])
#		printwritelist(filename,["p =",w.p(),"valuation =",w])
#		printwritelist(filename,["e =",e,"f =",f])
#		if D > 1:
#			printwritelist(filename,["Twisting by quadratic character of conductor",D])
#	if vLval != "" and filename != -1:
#		printwritelist(filename,["Valuation of central L-value:",vLval])


	ell1 = next_good_prime(1,A,w,depth,D,max_ell)
	if ell1 != -1:
		ell2 = next_good_prime(ell1,A,w,depth,D,max_ell)
	else:
		ell2 = -1

	ells = [ell1,ell2]

	if ell1 != -1 and ell2 != -1:
#		if filename != -1:
#			printwritelist(filename,["Trying:",(ell1,ell2)])
		while ells[-1] < max_ell:
			for a in range(len(ells)-1):
				print("Trying",(ells[a],ells[-1]))
				dn = delta(A,ells[a]*ells[-1],D,magic=magic) * pi^period_correction
				vdn = w(dn)
				m = In(A,w,ells[a]*ells[-1],D)
				if vdn >= m:
					print((ells[a],ells[-1]),": vanish")
					if filename != -1:
						printwritelist(filename,[(ells[a],ells[-1]),": vanish"])
				else:
					print((ells[a],ells[-1]),": non-vanishing")
					if filename != -1:
						printwritelist(filename,[(ells[a],ells[-1]),": non-vanish"])
#						printwritelist(filename,[""])
					return "Done"
			ell_new = next_good_prime(ells[-1],A,w,depth,D,max_ell)
			if ell_new == -1:
				break
			else:
				ells = ells + [ell_new]


	if filename != -1:
		if ell1 == -1 or ell2 == -1:
			printwritelist(filename,["**********************Failed to find non-zero delta**********************"])
			printwritelist(filename,["did not find two good primes <",max_ell])
		else:
			printwritelist(filename,["**********************Failed to find non-zero delta**********************"])
			printwritelist(filename,[""])
	return "Done"


def prove_large_image(A,w,max=100):
	N = A.level()
	p = w.p()
	k = A.weight()
	prime_divs = N.prime_factors()
	min_exp = min([N.valuation(q) for q in prime_divs])
	if min_exp > 1:
		return false 
	q = 2 
	while q < max:
		if N * p % q != 0:
			aq = ev(A,q)
			R.<x> = PolynomialRing(w.residue_field())
			g = x^2 - w.reduce(aq)*x + q^(k-1)
			if g.is_irreducible():
				return true
		q = next_prime(q)

	return false

def sign_of_FE(A):
	N = A.level()
	k = A.weight()
	r = (k-2)/2
	Lval = lamb_twist(A,r,0,1,1)
	if Lval != 0:
		return 1
	else:
		Ds = [D for D in range(1,100) if is_fundamental_discriminant(D) and kronecker_symbol(D,N) == -1];
		Lval = lamb_twist(A,r,0,1,Ds[0])
		if Lval != 0:
			return -1 
		else:
			Ds = [D for D in range(1,100) if is_fundamental_discriminant(D)];
			for j in range(len(Ds)):
				Lval = lamb_twist(A,r,0,1,Ds[j])
				if Lval != 0:
					return kronecker_symbol(Ds[j],N) 
			# else:
			# 	Ds = [D for D in range(1,100) if is_fundamental_discriminant(D) and kronecker_symbol(D,N) == 1];
			# 	Lval = lamb_twist(A,r,0,1,Ds[0])
			# 	if Lval != 0:
			# 		return 1
			# 	else:
			# 		Lval = lamb_twist(A,r,0,1,Ds[1])
			# 		if Lval != 0:
			# 			return 1
			# 		else:
			# 			assert Lval!=0,"failed to determine sign of FE"


def lower_bound_from_modsym(phi,w):
	if type(w) == sage.rings.integer.Integer:
		w = QQ.valuation(w)
	e = 1/w(w.uniformizer())
	val = phi.valuation(w)
	k = phi.weight()+2
	r = (k-2)/2

	M1 = min([ min([ (w(P.coef(j))-w(binomial(k-2,j))-val)*e for P in phi.data ]) for j in range(r+1) ])

	if k>2:
		M2 = min([ min([ (w(P.coef(j))-w(binomial(k-2,j))-val)*e+(j-r)*e for P in phi.data ]) for j in range(r+1,k-1) ])
	else:
		M2 = 10^10

	M3 = min([ (w(P.coef(r))-w(binomial(k-2,r))-val)*e for P in phi.data ]) 

	assert M3 == min(M1,M2), "problem with central coh period vs coh period. aborting"

	return min(M1,M2)




def form_deltas_in_fixed_weight_and_level(N,k,ps,max_ell,depth,Ds,require_large_image=true,skip_odd_rank=true,skip_unit_Lval=true,skip_Eisen=true,skip_CM=true,filename=-1,just_nonzero=true,skip_ll=true,level_lowering={},skip_odd_val=true):
	if isinstance(Ds,sage.rings.integer.Integer):
		Ds = [Ds]
	print("Working with newforms of level",N,"and weight",k)
	if require_large_image:
		prime_divs = N.prime_factors()
		if N > 1:
			min_exp = min([N.valuation(q) for q in prime_divs])
		else:
			min_exp = 0
		if min_exp > 1:
			print("Level has no prime factor with multiplicity 1 and so can't prove large image: skipping")
			return 
	sign = (-1)^((k+2)/2)
	M = ModularSymbols(N,k,sign).cuspidal_subspace().new_subspace()
	As = M.decomposition()
	print("-forming LMFDB labels")
	LMFDB_labels = LMFDB_labels_fixed_level(As)
	print("-there are",len(As),"Galois conjugacy classes of forms")
	for A in As:
		j = As.index(A)
		print("Working on GC class",LMFDB_labels[j],":",As.index(A)+1,"out of",len(As))
		print(A)

		print(" (forming reductions data to determine level-lowering)")
		for p in ps:
			if not p in level_lowering.keys():
				level_lowering[p] = {}
			if not N in level_lowering[p].keys():
				level_lowering[p][N] = []
			level_lowering[p][N] += reduce(A,p)
		magic = A.dual_eigenvector()
		K = magic.parent().base_ring()
		print(K)
		phi = form_modsym_from_decomposition(A)
		Lval = lamb_twist(A,(k-2)/2,0,1,1,magic=magic)
		if Lval != 0:
			sFE = 1 
		else:
			sFE = sign_of_FE(A)
		print("   sign of FE =",sFE)
		for p in ps:
			print("    Taking p =",p)
			v = QQ.valuation(p)
			ws = v.extensions(K)
			print("       -there are",len(ws),"prime(s) over",p)
			for w in ws:
				pi = w.uniformizer()
				kpi = w.residue_field()
				e = 1/w(pi)
				print("       --Working with valuation",ws.index(w)+1,"/",len(ws),"with e =",e,"and f =",w.residue_field().degree())
				eis = eisenstein(A,w)
				print("       --Is this case Eisenstein?:",eis)
				cm = CM(A,w)
				print("       --Is this case CM?:",cm)
				large_image = prove_large_image(A,w)
				print("       --Large image?:",large_image)
				ll = is_level_lowerable(A,w,30,level_lowering=level_lowering)
				print("       --Can the form be level-lowered?",ll)
				if N % p == 0:
					llt = is_level_lowerable_quad_twist(A,w,30,level_lowering=level_lowering)
					print("       --Can the form be twist-level-lowered?",llt)
				ap = ev(A,p)
				if N % p != 0:
					if w(ap)==0:
						print("       --good ordinary at this prime")
					else:
						print("       --good non-ordinary at this prime with slope",w(ap)/e)
				elif N.valuation(p) == 1:
					print("       --Steinberg")
				else:
					print("       --supercuspdial")
			header_written = false
			if (not skip_Eisen or not eis) and (not skip_ll or not ll) \
				and (not require_large_image or large_image):
				for D in Ds:
					if (not skip_odd_rank or sFE * kronecker_symbol(D,-N) == 1):
						if D % p != 0:
							print("    Twisting by",D)
							Lval = lamb_twist(A,(k-2)/2,0,1,D,magic=magic)
							period_correction = -lower_bound_from_modsym(phi,w) - phi.valuation(w,remove_binom=true)*e 	
							print("modsym valuation",phi.valuation(w,remove_binom=true))
							print("lower bound is",lower_bound_from_modsym(phi,w))
							print("period correction is",period_correction)
							print("L-val valuation",w(Lval))
							print("       --The central L-value has valuation:",w(Lval)*e+period_correction)
							if w(Lval)*e + period_correction > 0 or not skip_unit_Lval:
								if filename != -1 and not header_written:
									printwritelist(filename,[LMFDB_labels[j],"// p =",w.p()])
									if K.degree() > 1:
										printwritelist(filename,["Defined over number field with defining polynomial:",K.defining_polynomial()])
									else:
										printwritelist(filename,["Defined over QQ"])
									printwritelist(filename,["  choosing prime P over p with e =",e,"and f =",kpi.degree(),"\n  namely:",w])
									printwritelist(filename,["    (which in Sage's ordering is",ws.index(w)+1,"/",len(ws),")"])
									if N.valuation(p) == 0:
										if w(ap)==0:
											printwritelist(filename,["  good ordinary reduction at P"])
										else:
											printwritelist(filename,["  good non-ordinary reduction at P with slope",w(ap)])
									elif N.valuation(p) == 1:
										printwritelist(filename,["  Steinberg at P"])
									else:
										printwritelist(filename,["  supercuspidal at P"])
									if ll:
										printwritelist(filename,["CAN BE LEVEL-LOWERED!"])
									header_written = true									
								if not just_nonzero:
									#this is probably wrong with the periods an valuations and such
									compute_deltas(A,w,max_ell,depth,D,magic=magic,period_correction=period_correction,vLval=(w(Lval)-phi.valuation(w,remove_binom=true))*e,filename=filename)
									printwritelist(filename,[])
								else:
									if filename != -1:
										if Ds != [1]:
											printwritelist(filename,["***Twisting by quadratic character of conductor",D,"***"])
										printwritelist(filename,["The central L-value has valuation:",w(Lval)*e+period_correction])
										if Lval == 0:
											printwritelist(filename,["Note: central L-value vanishes, but sign of FE = 1"])
									if not ll:
										if not skip_odd_val or Lval==0 or (w(Lval)*e+period_correction) % 2 != 1:
											find_nonzero_delta(A,w,max_ell,depth,D,magic=magic,period_correction=period_correction,vLval=(w(Lval)-phi.valuation(w,remove_binom=true))*e,filename=filename)
										else:
											print("Skipping odd valuation L-value: ",end="")
											if llt:
												print("can be level-lowered after twist at p")
											else:
												print("CANNOT be level-lowered after twist at p")
											if llt:
												printwritelist(filename,["Skipping odd valuation L-value: can be level-lowered after twist at p"])
											else:
												printwritelist(filename,["Skipping odd valuation L-value: CANNOT be level-lowered after twist at p"])
									else:
										find_nonzero_delta(A,w,max_ell/10,depth,D,magic=magic,period_correction=period_correction,vLval=(w(Lval)-phi.valuation(w,remove_binom=true))*e,filename=filename)
									printwritelist(filename,[])
							else:
								print("       --Skipping because L-value is a unit")
				if header_written and filename != -1:
					printwritelist(filename,[""])

			else:
				if skip_odd_rank and sFE == -1:
					print("       --Skipping because sign of FE = -1")					
				if skip_Eisen and eis:
					print("       --Skipping because it is Eisenstein")
				if skip_ll and ll:
					print("       --Skipping because of level-lowering")
				if require_large_image and not large_image:
					print("       --Skipping because we can't prove large image")
			print()	

	del M
	return 

def eisenstein(A,w,max_check=50):
	N = A.level()
	k = A.weight()
	p = w.p()

	#checks congruence with E_k
	bool = true
	q = 2
	while q < max_check and bool:
		if N*p % q != 0:
			bool = bool and w(ev(A,q)-1-q^(k-1)) > 0
		q = next_prime(q)

	if bool:
		return bool

	#checks congruence with E_k(chi,chi) for chi quadratic
	for d in N.divisors():
		Gd = DirichletGroup(d,QQ)
		for chi in Gd.gens():
			bool = true
			q = 2
			while q < max_check and bool:
				if N*p % q != 0:
					bool = bool and w(ev(A,q)-chi(q)-chi(q)*q^(k-1)) > 0
				q = next_prime(q)
			if bool:
				return bool

	return bool

def CM(A,w,max_check=50):
	N = A.level()
	k = A.weight()
	p = w.p()

	#checks if f otimes chi = f (mod p)
	for d in N.divisors():
		Gd = DirichletGroup(d,QQ)
		for chi in Gd.gens():
			bool = true
			q = 2
			while q < max_check and bool:
				if N*p % q != 0:
					bool = bool and w(ev(A,q)-chi(q)*ev(A,q))>0
				q = next_prime(q)
			if bool:
				return bool

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