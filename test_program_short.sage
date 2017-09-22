load "master.sage"

p=3
N=5
E=EllipticCurve('15a')
phi=form_modsym_from_elliptic_curve(E)
Phi=phi.lift_to_OMS(3,10)
#creates an overconvergent lifting of phi accurate module Fil^10
if Phi.valuation()<0:
	Phi=Phi.scale(p^(-Phi.valuation()))
#kills off any denominators
Phi=Phi.hecke(p)
#solving the difference equation potentially causes problems that applying Hecke once clears up
Phi=Phi.hecke(p)-Phi
#currently the lifting of classical modular symbols only works up to some Eisenstein error.  Here we apply U_p-1 to  kill off this error
print "Iterating U_p"
for j in range(15):
	print j
	Phi=Phi.hecke(p); Phi
#The result after applying U_p should converge to a Hecke-eigensymbol lifting phi
print "Here's Phi + Phi|U_p"
print Phi+Phi.hecke(p)
#The result should be all zeroes -- i.e. Phi | U_p = - Phi since a_3(E)=-1
R.<w>=PolynomialRing(QQ)

print "Lifting to a family of OMSs"
Phis=Phi.lift_to_modsym_dist_fam(15,w)
print "done!"
if Phis.valuation()<0:
	Phis=Phis.scale(p^(-Phis.valuation()))
Phis=Phis.hecke(p)
#Same deal as before
Phis=Phis.hecke(p)-Phis
#Same deal with Eisenstein stuff as before
print "Iterating U_p in families"
for j in range(15):
	print j
	Phis=Phis.hecke(p); Phi
#Result should be an eigen-family!
print Phis
