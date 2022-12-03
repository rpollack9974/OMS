## Enter elliptic curve E, prime p and accuracy M, and this program begins to find the critical slope p-adic L-function of E 

E = EllipticCurve('24a1')
p = 5
M = 12
print "Working with the elliptic curve",E,"and the prime",p
print "Forming classical symbol"
phiE = form_modsym_from_elliptic_curve(E)
phiE = phiE.minus_part().scale(1/2)
print "p-stabilizing"
phi_beta = phiE.p_stabilize_critical(p,E.ap(p),M+2)
e = phi_beta.valuation(p)
phi_beta = phi_beta.scale(p^(-e)) # temporarily scales phi_beta to be integral
R.<x> = PolynomialRing(pAdicField(p,2*M))
rs = (x^2-E.ap(p)*x+p).roots()
if rs[0][0] % p == 0:
	beta = ZZ(rs[0][0])
else:
	beta = ZZ(rs[1][0])
print "Lifting to OMS"
Phi = phi_beta.lift_to_OMS(p,M)
phi_beta = phi_beta.scale(p^e) # returns phi_beta is the original correct form
print "Smoothing by U_p"
Phi = Phi.hecke(p)
q = 2
while (E.conductor()*p) % q == 0:
	q = next_prime(q)
print "Using T_",q," to kill Eisen part"
Phi = Phi.hecke(q) - Phi.scale(q+1)
scale *= E.ap(q)-(q+1)
Phi = Phi.minus_part()
print Phi.check_loop().normalize()
for r in range(M+1):
	print "Applying Up/beta",(r,M+1)
	Phi = Phi.up(p,beta)
	print Phi.check_loop().normalize()
	print "lots of zeros good"
Phi = Phi.hecke(p) - Phi.scale(p) ## kill of critical Eisenstein contribution
Phi = Phi.hecke(p) - Phi.scale(p) ## kill of critical Eisenstein contribution
Phi = Phi.hecke(p) - Phi.scale(p) ## kill of critical Eisenstein contribution ## how many times???
Psi = Phi.up(p,beta) - Phi
#Psi = Psi.hecke(p) - Psi.scale(p) ## kill of critical Eisenstein contribution
if Psi.valuation() > 0:
	e = Psi.valuation()
	print "valuation",e
	Psi = Psi.scale(p^(-Psi.valuation()))
	Psi = Psi.change_precision(Psi.num_moments()-e)
Psis = [Psi]
for r in range(M):
	print "Applying U_",p,"/beta ",(r,M)
	Psis += [Psis[len(Psis)-1].hecke(2)]
vs = []
for j in range(len(Psis)):
	vs += [[Psis[j].data[a].moment(1) for a in range(len(Psi.manin.gens))]]
A = Matrix(vs)
#print row_reduce_mod_pn(A,p,10)


#find (by hand) good relations in the matrix A to kill off Psi
#and record them in the polynomial f

Phi2 = Phi.up_by_poly(p,beta,f,verbose=True)
## Phi2 should now be eigen (check this by hand)

# rescaling now Phi2 to specialize to phi_beta
j = 0
while phi_beta.data[j].coef(0) == 0:
	j = j + 1
ratio = phi_beta.data[j].coef(0) / Phi2.data[j].moment(0) 
Phi2 = Phi2.scale(ratio)



