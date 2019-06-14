def padic_family(p,N,r,ap_list,acc,verbose=False):
	"""computes the OMS family with prime p, tame level N, component r of weight space, character char,
		cut out by maximal ideal described by ap_list = [(ell_1,a_{ell_1}),...]"""
	R = PolynomialRing(QQ,'w')
	w = R.0
	one = trivial_character(1)
	if verbose:
		print "Forming random symbol"
	Phis = random_OMS_fam(p,N,one,acc,r,w)
	if verbose:
		print "Killing Eisenstein family"
	Phis = Phis.hecke(p) - Phis
	if verbose:
		print "Projecting to ordinary subspace"
	for j in range(acc+2):
		Phis = Phis.hecke(p)
		if verbose:
			print "finished ",j+1,"out of",acc+2
	if verbose:
		print "Killing off other families"
	M = ModularSymbols(N*p,r+2,1,GF(p))
	U = M.hecke_operator(p)
	M = (U^M.dimension()).image()
	RZ = PolynomialRing(QQ,'x')
	for t in ap_list:
		q = t[0]
		aq = t[1]
		if verbose:
			print "Using prime ",q
		fq = M.hecke_polynomial(q)
		x = fq.parent().0
		m = fq.valuation(x-aq)
		fq = fq.quo_rem((x-aq)^m)[0]
		fq = RZ(fq)
		for j in range(acc+2):
			if verbose:
				print "finished ",j+1,"out of",acc+2
				Phis = Phis.hecke_by_poly(q,fq,verbose=verbose)
	return Phis
