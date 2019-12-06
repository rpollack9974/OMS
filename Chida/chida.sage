
attach("master.sage")

def prepare_for_critical_OMS(E,p,M):
	print "forming classical symbol"
	phiE = form_modsym_from_elliptic_curve(E)
	phiE = phiE.minus_part().scale(1/2)
	print "p-stabilizing"
	phi = phiE.p_stabilize_critical(p,E.ap(p),M+2)
	scale = p^-phi.valuation(p)
	phi = phi.scale(p^-phi.valuation(p))
	print "lifting to OMS"
	Phi = phi.lift_to_OMS(p,M)   

	print "killing Eisenstein"
	Phi = (Phi.hecke(p) - Phi).scale((E.ap(p)-1)^(-1))
	Phi = (Phi.hecke(p) - Phi.scale(p)).scale((E.ap(p)-p)^(-1))
	print "projecting to slope 1 subspace"
	R = PolynomialRing(pAdicField(p,2*M),'x')
	x = R.gen()
	f = x^2-E.ap(p)*x+p
	roots = f.roots()
	if roots[0][0].valuation()==1:
		beta = ZZ(roots[0][0])
	else:
		beta = ZZ(roots[1][0])
	for j in range(M):
		print j
		Phi = Phi.up(p,beta)

	return Phi

def value_at_0(Phi,p,beta):
	M = Phi.num_moments()
	ans = 0
	for a in range(1,p):
		Da = Matrix(2,2,[1,a,0,p])
		for j in range(M):
			ans += (-1)^j * ZZ(a)^(-j-1) * ZZ(p)^j / beta * Phi.eval(Da).moment(j)
	return ans
