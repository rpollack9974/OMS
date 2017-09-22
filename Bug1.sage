#Form with CM by Gaussian integers
load "master.sage"

f = Newforms(32,4)[1]
D = -4
p = 5
M = 3
k = 4 - 2
R.<w> = PolynomialRing(QQ)

phi = form_modsym_from_decomposition(f.modular_symbols(1))
ap = f.coefficients([p])[0]
QQp = Qp(p)
ZmodpM = Zmod(p ** (M + 1))
phi_p = phi.p_stabilize_ordinary(p, ap, max(25, M+5))
Px = PolynomialRing(QQp, 'x')
x = Px.gen()
fx = x^2 - ap * x + p^(k+1)
rs = fx.roots()
if rs[0][0].is_unit():
	alpha_p = rs[0][0]
else:
	alpha_p = rs[1][0]

b, c = phi_p.is_Tq_eigen(p, p, M)
print b
print QQp(c)
print alpha_p

print "Lifting to an overconvergent modular symbol with", M, "moments (with some Eisenstein error)"
Phi = phi_p.lift_to_OMS(p, M)
if Phi.valuation()<0:
	Phi = Phi.scale(p^(-Phi.valuation()))

Phi = Phi.hecke(p)
Phi=Phi.hecke(3)-Phi.scale(1+3^(k+1))

for j in range(M+1):
	print j+1, "-th application of U_p"
	Phi=Phi.hecke(p)

print "Checking if OMS is a U_p-eigenform with unit eigenvalue"
tempo = QQp(Phi.is_Tq_eigen(p))
print tempo
print alpha_p

tempo3 = Phi.is_Tq_eigen(3)
print tempo3

print "Lifting to a family"
Phis = Phi.lift_to_modsym_dist_fam(w)
Phis = Phis.scale(p^(-Phis.valuation()))
for j in range(M+1):
	print j
	Phis = Phis.hecke(p)
	Phis = Phis.scale(p^(-Phis.valuation()))