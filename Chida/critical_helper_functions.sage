def beta_root(E,p,M):
	"""returns slope 1 root of hecke polynomial"""
	assert E.ap(p)%p!=0, "supersingular"
	R = PolynomialRing(pAdicField(p,2*M),'x')
	x = R.gen()
	rs = (x^2-E.ap(p)*x+p).roots()
	if rs[0][0].residue() == 0:
		beta = ZZ(rs[0][0])
	else:
		beta = ZZ(rs[1][0])
	return beta

def slope1_lift(E,p,M):
	print("Working with the elliptic curve",E,"and the prime",p)
	print("Forming classical symbol")
	phiE = form_modsym_from_elliptic_curve(E)
	phiE = phiE.minus_part().scale(1/2)
	print("p-stabilizing")
	phi_beta = phiE.p_stabilize_critical(p,E.ap(p),M+2)
	e = phi_beta.valuation(p)
	phi_beta = phi_beta.scale(p^(-e)) # temporarily scales phi_beta to be integral
	beta = beta_root(E,p,M)
	print("Lifting to OMS")
	Phi = phi_beta.lift_to_OMS(p,M)
	phi_beta = phi_beta.scale(p^e) # returns phi_beta to the original correct form
	print("Smoothing by U_p")
	Phi = Phi.hecke(p)
	q = 2
	while (E.conductor()*p) % q == 0:
		q = next_prime(q)
	print("Using T_",q," to kill Eisen part")
	Phi = Phi.hecke(q) - Phi.scale(q+1)
	print(Phi.check_loop().normalize())
	print("Using U_",p," to kill critical Eisen part")
	Phi = Phi.hecke(p) - Phi.scale(p)
	for r in range(M+1):
		print("Applying Up/beta",(r,M+1))
		Phi = Phi.up(p,beta)
		print(Phi.check_loop().normalize())
		print("lots of zeros good")

	return Phi

def hecke_iterates(Psi,q,M,E):
	"""Psi = Phi | U_p - Phi. returns [Psi | T_q^n]"""
#	Phi = Phi.hecke(p) - Phi.scale(p) ## kill of critical Eisenstein contribution
#	Phi = Phi.hecke(p) - Phi.scale(p) ## kill of critical Eisenstein contribution
#	Phi = Phi.hecke(p) - Phi.scale(p) ## kill of critical Eisenstein contribution ## how many times???
	p = Psi.p()
	beta = beta_root(E,p,M)
	#Psi = Psi.hecke(p) - Psi.scale(p) ## kill of critical Eisenstein contribution
	if Psi.valuation() > 0:
		e = Psi.valuation()
		print("valuation",e)
		Psi = Psi.scale(p^(-Psi.valuation()))
		Psi = Psi.change_precision(Psi.num_moments()-e)
	Psis = [Psi]
	for r in range(M):
		print("Applying hecke at",q,(r,M))
		if q != p:
			Psis += [Psis[len(Psis)-1].hecke(q)]
		else:
			Psis += [Psis[len(Psis)-1].up(p,beta)]
	vs = []
	for j in range(len(Psis)):
		vs += [[Psis[j].data[a].moment(1) for a in range(len(Psi.manin.gens))]]
	A = Matrix(vs)
	#print row_reduce_mod_pn(A,p,10)
	return A

def extract_poly_relation(A,row):
	R = PolynomialRing(QQ,'x')
	x = R.gen()
	f = R(0)
	for c in range(A.ncols()):
		f += ZZ(A[row,c]) * x^(c)
	return f

def find_slope1_kernel_relation(E,p,M):
	Phi = slope1_lift(E,p,M)
	beta = beta_root(E,p,M)
	Psi = Phi.up(p,beta) - Phi
	e = Psi.valuation()
	Psi = Psi.scale(1/p^e)
	Psi = Psi.change_precision(Psi.num_moments()-e)
	q = 2
	while E.discriminant() % q == 0:
		q = next_prime(q)
	A = hecke_iterates(Psi,q,M,E)

	return Phi,A,q

def matrix_valuation(A,p):
	v = A.list()
	return min([v[a].valuation(p) for a in range(len(v))])

#f is the polynomial such that f(T_q) kills kernel of specialization 
def find_eigen_slope1(Phi,f,q,E):
	p = Phi.p()

	beta = beta_root(E,p,M)
	if q != p:
		Phi = Phi.hecke_by_poly(q,f,verbose=True)
	else:
		Phi = Phi.up_by_poly(p,beta,f,verbose=True)

	return Phi

def did_it_work(Phi):
	print("Did it work?")
	print(Phi.is_Tq_eigen(2))
	print(Phi.is_Tq_eigen(3))
	print(Phi.is_Tq_eigen(5))		

	return ""

def renormalize(Phi,E):
	phiE = form_modsym_from_elliptic_curve(E)
	phiE = phiE.minus_part().scale(1/2)
	phi_beta = phiE.p_stabilize_critical(p,E.ap(p),Phi.data[0].num_moments()+2)

 	a = 0
	phi = Phi.specialize()
	while phi.data[a].coef(0) == 0:
		a = a + 1
	c = phi_beta.data[a].coef(0) / phi.data[a].coef(0)
	Phi = Phi.scale(c)

	return Phi

def value_at_0(Phi,p,beta):
	M = Phi.num_moments()
	ans = 0
	for a in range(1,p):
		Da = Matrix(2,2,[1,a,0,p])
		for j in range(M):
			ans += (-1)^j * ZZ(a)^(-j-1) * ZZ(p)^j / beta * Phi.eval(Da).moment(j)
	return ans

def row_reduce_mod_pn2(B,p,n):
	A = copy(B)
	rows = A.nrows()
	cols = A.ncols()
	M = MatrixSpace(Integers(p^n),rows,cols)
	A = M(A)
	curr_col = 0
	curr_row = 0
	E = Matrix(Integers(p^n),rows,rows,1)
	while (curr_row < rows) and (curr_col < cols):
#		print A
#		print E
#		print "(curr_row,curr_col)",(curr_row,curr_col)
		v = [A[r][curr_col] for r in range(curr_row,rows)]
#		print v
		vals = [ZZ(v[j]).valuation(p)  for j in range(len(v))]
#		print vals
		m = min(vals) 
#		print "Out: (r,c,m)=",(curr_row,curr_col,m)
		if m < Infinity:
#			print "In: (r,c,m)=",(curr_row,curr_col,m)
			least_ind = curr_row
			while vals[least_ind-curr_row] > m:
				least_ind = least_ind + 1
#			print "least_ind = ",least_ind
			A.swap_rows(curr_row,least_ind)
			F = Matrix(Integers(p^n),rows,rows,1)
			F.swap_rows(curr_row,least_ind)
			E = F * E	
			for r in range(rows):
				if r != curr_row:
#					print "1",A[curr_row,curr_col]
#					print "2",ZZ(A[curr_row,curr_col])
#					print "3",ZZ(A[r,curr_col])
					c = -ZZ(A[r,curr_col])/ZZ(A[curr_row,curr_col])
					if c.valuation(p) >= 0:
						A.add_multiple_of_row(r,curr_row,c)
						F = Matrix(Integers(p^n),rows,rows,1)
						F[r,curr_row] = c
						E = F * E
			curr_row = curr_row + 1
		curr_col = curr_col + 1
		#print A
		#print "----"
		#print E*A
		#print "-------------"
	return A,E,curr_row


def chida(E,p,M):
	Phi,A,q = find_slope1_kernel_relation(E,p,M)
	e = matrix_valuation(A,p)
	d = floor(M/2)
	B,C,r1 = row_reduce_mod_pn2(A/p^e,p,d)
	B,C,r2 = row_reduce_mod_pn2(A/p^e,p,d+1)
	while r1 == r2:
		d = d + 1
		B1,C1,r1 = row_reduce_mod_pn2(A/p^e,p,d)
		B2,C2,r2 = row_reduce_mod_pn2(A/p^e,p,d+1)
	r = r1
	print(B1)
	print(C1)
	print("Using row",r)
	print("Should be zero:",B1[r])
	f = extract_poly_relation(C1,r)

	print("Killing off kernel of specialization")
	beta = beta_root(E,p,M)
	Phi_eigen = find_eigen_slope1(Phi,f,q,E)
	did_it_work(Phi_eigen)
	Phi_eigen = renormalize(Phi_eigen,E)
	print("And here's the value...")
	val = value_at_0(Phi_eigen,p,beta) 
	v = val.valuation(p)
	if v < 0:
		val_norm = ((val * p^(-v)) % (p^(Phi_eigen.num_moments()))) * p^v
	else:
		val_norm = val % (p^(Phi_eigen.num_moments()))
	return val_norm