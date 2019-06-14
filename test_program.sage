load master.sage

p=3
N=5
E=EllipticCurve('15a')
phi=form_modsym_from_elliptic_curve(E)
Phi=phi.lift_to_OMS(3,20)
#creates an overconvergent lifting of phi accurate module Fil^10
if Phi.valuation()<0:
	Phi=Phi.scale(p^(-Phi.valuation()))
#kills off any denominators
Phi=Phi.hecke(p)
#solving the difference equation potentially causes problems that applying Hecke once clears up
Phi=Phi.hecke(p)-Phi
#currently the lifting of classical modular symbols only works up to some Eisenstein error.  Here we apply U_p-1 to  kill off this error
print "Iterating U_p"
for j in range(20):
	print j
	Phi=Phi.hecke(p); Phi
#The result after applying U_p should converge to a Hecke-eigensymbol lifting phi
print "Here's Phi + Phi|U_p"
print Phi+Phi.hecke(p)
#The result should be all zeroes -- i.e. Phi | U_p = - Phi since a_3(E)=-1
R.<w>=PolynomialRing(QQ)

print "Lifting to a family of OMSs"
Phis=Phi.lift_to_modsym_dist_fam(w)
print "done!"
if Phis.valuation()<0:
	Phis=Phis.scale(p^(-Phis.valuation()))
Phis=Phis.hecke(p)
#Same deal as before
Phis=Phis.hecke(p)-Phis
#Same deal with Eisenstein stuff as before
print "Iterating U_p in families"
for j in range(20):
	print j
	Phis=Phis.hecke(p); Phi
#Result should be an eigen-family!
print Phis
print "forming weight 0 (aka weight 2) specialization)"
Phi0=Phis.specialize(0).normalize()
print "testing Hecke eigen at 2 (lots of zeroes are good)"
Phi0.is_Tq_eigen(2)
print "testing Hecke eigen at 3 (lots of zeroes are good)"
Phi0.is_Tq_eigen(3)
print "testing Hecke eigen at 5 (lots of zeroes are good)"
Phi0.is_Tq_eigen(5)

print "forming weight 6 (aka weight 8) specialization)"
Phi6=Phis.specialize(6).normalize()
#Phi6 is the weight 6 (aka weight 8) specialization of Phis
print "testing Hecke eigen at 2 (lots of zeroes are good)"
Phi6.is_Tq_eigen(2)
print "testing Hecke eigen at 3 (lots of zeroes are good)"
Phi6.is_Tq_eigen(3)
print "testing Hecke eigen at 5 (lots of zeroes are good)"
Phi6.is_Tq_eigen(5)

print "forming weight 1000 (aka weight 1002) specialization)"
Phi1000=Phis.specialize(1000).normalize()
print "testing Hecke eigen at 2 (lots of zeroes are good)"
Phi1000.is_Tq_eigen(2)
print "testing Hecke eigen at 3 (lots of zeroes are good)"
Phi1000.is_Tq_eigen(3)
print "testing Hecke eigen at 5 (lots of zeroes are good)"
Phi1000.is_Tq_eigen(5)
