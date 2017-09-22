load "master.sage"

print walltime()

p = 7
k = 2-2
M = 12
max_k = 12

QQp = Qp(p)
ZmodpM = Zmod(p ** (M + 1))

E = EllipticCurve('27a1') #Has CM by Q(sqrt(-3)) and 7 is the smallest ordinary prime
print "Forming modular symbols (in sense of Stevens) attached to E"
phi = form_modsym_from_elliptic_curve(E)
#print phi

print "p-stabilize"
ap = E.ap(p)
phi_p = phi.p_stabilize_ordinary(p, ap, 25)
#print phi_p

#N = E.conductor()
Px = PolynomialRing(QQp, 'x')
x = Px.gen()
fx = x^2 - ap * x + p
rs = fx.roots()
if rs[0][0].is_unit():
	alpha_p = rs[0][0]
else:
	alpha_p = rs[1][0]

print "Checking if p-stabilisation is U_p-eigenform with unit eigenvalue"
b, c = phi_p.is_Tq_eigen(p, p, M)
print b
print QQp(c)
print alpha_p

print "Lifting to an overconvergent modular symbol with", M, "moments (with some Eisenstein error)"
Phi = phi_p.lift_to_OMS(p, M)
if Phi.valuation()<0:
	Phi = Phi.scale(p^(-Phi.valuation()))
print Phi

Phi = Phi.hecke(p)

Phi=Phi.hecke(5)-Phi.scale(6)

for j in range(M+1):
	print j+1, "-th application of U_p"
	Phi=Phi.hecke(p); print Phi

print "Checking if OMS is a U_p-eigenform with unit eigenvalue"
print QQp(Phi.is_Tq_eigen(p))
print alpha_p
print "Checking if OMS is a T_3 and T-5 eigenform"
print Phi.is_Tq_eigen(3)
print Phi.is_Tq_eigen(5)

print "Lifting to a family"
R.<w> = PolynomialRing(QQ)
Phis = Phi.lift_to_modsym_dist_fam(w)
Phis = Phis.scale(p^(-Phis.valuation()))
for j in range(M+1):
	print j
	Phis = Phis.hecke(p)
	Phis = Phis.scale(p^(-Phis.valuation()))

print "Killing Eisenstein part"
Phis = Phis.hecke(5)-Phis.scale(6)

ap_fam = Phis.is_Tq_eigen(p)
g = ((ZmodpM((1+p))^k - 1).lift()/p) % (p^M)
apk = ap_fam.substitute(w = g)
print "The U_p-eigenvalue is"
print QQp(apk)
print "up to O(p^{0})".format(M - Phis.valuation())
print alpha_p
print ap_fam

#Compare L-invariants of symmetric square
aps = []
ds = []
ks = []
for j in range(max_k):
	delta = (p - 1) * (ZmodpM(p) ** (j + 1)).lift()
	kj = k + delta
	ks.append(kj)
	gj = (((ZmodpM(1 + p) ** (kj)) - 1).lift()/p) % (p^M)
	print "j =", j, "gj =", gj
	aps.append(ap_fam.substitute(w = QQp(gj)))
	ds.append((aps[-1] - QQp(alpha_p)) / QQp(delta))

print "Using the Ferrero-Greenberg/Gross-Koblitz formula, the L-invariant of the symmetric square of", E.label(), "is:"
D = E.cm_discriminant()
L_FGGK = Linvariant(QuadraticField(D), p)
print L_FGGK
print "Using the Hida/Harron derivative formula yields these successive approximations:"
for j in ds:
	print -2 * j / alpha_p
print "Check for constancy of L-invariant (i.e. check that ap_fam is exponential):"
for j in range(max_k):
	print aps[j].log(0)/(ks[j] + 1)