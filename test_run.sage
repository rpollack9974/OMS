load "master.sage"

N = 5
p = 3
#G = DirichletGroup(p,QQ)
#chi = G.0
r = 0
M = 10

Phi = random_OMS(N,p,r,M,chi)
for j in range(M+1):
	print j
	Phi = Phi.hecke(p)
#Phi = Phi.hecke(2) - Phi.scale(1+chi(2)*2^(r+1))
Phi = Phi.hecke(2) - Phi.scale(3)

print "FINISHED FORMING OMS.  CHECKING IF ITS AN EIGENSYMBOL"

Phi.is_Tq_eigen(2)
Phi.is_Tq_eigen(3)
Phi.is_Tq_eigen(5)
Phi.is_Tq_eigen(7)

print "STARTING TO FORM HIDA FAMILY"

Phis=Phi.lift_to_modsym_dist_fam(w)
Phis=Phis.scale(p^(-Phis.valuation()))
for j in range(M+1):
	print j
	Phis = Phis.hecke(p)
Phis = Phis - Phis.hecke(p)

ap = Phis.is_Tq_eigen(p)
g = (((1+p)^(k) - 1)/p) % (p^M)
print "The U_p-eigenvalue is", ap.substitute(w=g) % (p^(M-Phis.valuation())), "modulo ",p,"^",M-Phis.valuation()