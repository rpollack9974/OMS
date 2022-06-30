import sys

S = PolynomialRing(PolynomialRing(QQ,'w'),'T')

def collect_padic_Lfunctions(Phis,D,verbose=false):
	Ls = []
	for r in range(p-1):
		if verbose:
			print("twist =",D," component =",r)
		Ls += [Phis.pLfunction(r=r,quad_twist=D)]
	return Ls


def analyze_pLs(D,Phis_list,verbose=true):
	D = ZZ(D)
	p = Phis_list[0].p()
	comp = Phis_list[0].disc()
	Ls = []
	S = PolynomialRing(PolynomialRing(QQ,'w'),'T')
	for i in range(p-1):
		num = 0
		done = false
		while num < len(Phis_list) and not done:
			Phis = Phis_list[num]
			if verbose:
				print("Working with twist ",D,", component ",i,"using",Phis.num_moments(),"moments")
			L = Phis.pLfunction(r=i,quad_twist=D)
			d = S(L).degree()
			mu = mu_inv(L.substitute(w=1)/p^Phis.valuation(),p)
			if mu > 0:
				print("mu positive so probably not enough accuracy")
				num += 1
			else:
				lam = lambda_inv(L.substitute(w=1)/p^Phis.valuation(),p)
				if 2*i % (p-1) == comp and lam % 2 == 1:
					toric_bound = p/(p-1)
				else:
					toric_bound = 1
				if lam == 0:
					print("--lambda = 0 so nothing to check here")
					print("PASSED")
					done = true
				elif lam == 1 and 2*i % (p-1) == comp:
					print("--lambda = 1 because of FE so nothing to check here")
					print("PASSED")
					done = true
				else:
					print("--lambda =",lam)
					done = false
					too_high = false
					n = 1
					while not done and not too_high:
						print("working at level",p^n)
						K = CyclotomicField(p^n)
						z = K.gen()
						pp = K.factor(p)[0][0]
						vals = []
						maxs = []
						for b in range(1,p-1):
							for a in range(p^n):
								t = S(L).substitute(T=z-1).substitute(w=(z^a-1+b*p)/p)
								val = t.valuation(pp)/pp.ramification_index() - Phis.valuation()
								vals += [val]
							m = max(vals)	
							maxs += [m]
						m = max(maxs)
						error_bound = d/(p^(n-1)*(p-1))
						if m < error_bound and m < toric_bound:
							print("Passed! Max valuation is",m,", toric bound is",toric_bound,"and error bound is",error_bound)
							print("PASSED")
							done = true
						elif m >= error_bound:
							print("Failed: not enough accuracy.  Max valuation is",m,"and error bound is",error_bound)
							n += 1
						else:
							print("Failed: max val too high.  Max valuation is",m," and toric bond is ",toric_bound)
							n += 1
						if n > 4:
							num += 1
							if num < len(Phis_list):
								print("Going to more accurate family.")
							else:
								print("giving up!")
							too_high = true
		if not done:
			print("*************************FAILED!!!***************************")
		print("")
	return "done"

def run_me(filename,minD,maxD,Phis_list,level,step=1,log=true):
	old_stdout = sys.stdout

	for d in range(minD,maxD,step):
		if is_fundamental_discriminant(d) and gcd(d,level)==1:
			print("Working on twist d=",d)
			if log:
				log_file = open(filename,"a")
				sys.stdout = log_file			
			analyze_pLs(d,Phis_list)
			print("----------------------------------------")
			if log:
				sys.stdout = old_stdout
				log_file.close()

def OMSfamily_extends_to_unit_disc(Phis):
	for a in range(len(Phis.data)):
		for b in range(Phis.data[a].num_moments()):
			f = Phis.data[a].moment(b)
			if not extends_to_unit_disc(f):
				return false
	return true

def divide_by_pw(Phis):
	data = []
	R = Phis.data[0].moment(0).parent()
	w = R.gen()
	for a in range(len(Phis.data)):
		d = Phis.data[a]
		v = Phis.data[a].moments 
		vv = copy(v)
		for b in range(len(vv)):
			vv[b] = R(vv[b]/(p*w))
		data += [dist_fam(d.p,d.deg,d.disc(),vv,d.char())]

	ans = modsym_dist_fam(Phis.level(),data,Phis.manin)
	ans = ans.change_precision(ans.num_moments()-1)
	ans = ans.change_deg(ans.deg()-1)

	return ans