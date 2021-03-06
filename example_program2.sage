load "master.sage"

p=11
print "Starting with the elliptic curve X_0(11)"
E=EllipticCurve('11a')
print "Forming modular symbols (in sense of Stevens) attached to E"
phi=form_modsym_from_elliptic_curve(E)
print phi
print "This modular symbol is stored by its value on 3 (particular) degree zero divisors which generator over Z[Gamma_0(11)]"

print
print "Lifting to an (11-adic) OMS with approximately 10 moments"
Phi=phi.lift_to_OMS_eigen(11,10)
print "The result is:"
print Phi
print "An OMS of level 11 is represented by 3 distributions -- here each distribution is a vector of length 10 representing the moments"
print
print "The p-adic L-function in the T-variable is:"
print pLfunction(Phi,1)
print "(Note the trivial zero.)"
print
print "Now lifting to a family of modular symbols"
R.<w>=PolynomialRing(QQ)
Phis=Phi.lift_to_modsym_dist_fam(w)
Phis=Phis.normalize()
print Phis
print "Here a family of OMS (of level 11) is represented by 3 distributions -- here each distribution is a vector of length 10 representing the moments and each moment is a power series in w which is the weight variable -- specializing w to (1+p)^k-1 then an OMS of weight k"
print
print "Now we form a Hecke-eigen family:"
print "-Killing Eisenstein part"
Phis=Phis-Phis.hecke(11)
#This will force the family to vanish in weight 0 but this is not a problem because we could then divide the whole OMS family by w
print "-iterating U_p"
for a in range(10):
	Phis=Phis.hecke(11)
	print a,"out of",10
print "The result is:"
print Phis
print
print "Now specializing to weight 12 (this should be the classical Delta function!) -- call the specialization Del"
Del=Phis.specialize(10).normalize()
print Del
print
print "Now computing Del|T_2 - (-24)*T_2 (result should be all zeroes!)"
print (Del.hecke(2) - Del.scale(-24)).normalize()
print
print "Now computing Del|T_3 - (252)*T_3 (result should be all zeroes!)"
print (Del.hecke(3) - Del.scale(252)).normalize()
print
print "computing a_11 as a function of w"
A=Phis.hecke(11)
a=A.data[0].moment(0)
b=Phis.data[0].moment(0)
S.<z>=PowerSeriesRing(QQ)
a11=((S(a)/S(b)).substitute(z=w))
T.<q>=PolynomialRing(Integers(11^10))
print "a11 =",R(T(a11))+O(w^7)



