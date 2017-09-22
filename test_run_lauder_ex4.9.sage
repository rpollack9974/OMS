load "master.sage"

N = 1
p = 23
k = (1234728-2)
r = k % (p-1)
M = 11
R.<w> = PolynomialRing(QQ)

Phi = random_OMS(N,p,r,M)
for j in range(M+1):
	print j
	Phi = Phi.hecke(p)
Phi = Phi - Phi.hecke(p)

print "FINISHED FORMING OMS.  CHECKING IF ITS AN EIGENSYMBOL"

Phi.is_Tq_eigen(5)
Phi.is_Tq_eigen(7)

print "STARTING TO FORM HIDA FAMILY"

Phis=Phi.lift_to_modsym_dist_fam(w)
for j in range(M+1):
	print j
	Phis = Phis.hecke(p)
Phis = Phis - Phis.hecke(p)

ap = Phis.is_Tq_eigen(p)
g = (((1+p)^(k) - 1)/p) % (p^M)
print "The U_p-eigenvalue is", ap.substitute(w=g) % (p^(M-Phis.valuation())), "modulo ",p,"^",M-Phis.valuation()