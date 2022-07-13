import sys

S = PolynomialRing(PolynomialRing(QQ,'w'),'T')

def collect_padic_Lfunctions(Phis,D,verbose=false):
	Ls = []
	for r in range(p-1):
		if verbose:
			print("twist =",D," component =",r)
		Ls += [Phis.pLfunction(r=r,quad_twist=D)]
	return Ls


def analyze_pLs(D,Phis_list,cyc_comp=None,full_search=false,verbose=true,write_to_file=False,filename=None):
	D = ZZ(D)
	p = Phis_list[0].p()
	vp = QQ.valuation(p)
	comp = Phis_list[0].disc()
	Ls = []
	S = PolynomialRing(PolynomialRing(QQ,'w'),'T')
	toroidal_bound = 1
	if cyc_comp == None:
		R = range(p-1)
	else:
		R = [cyc_comp]
	for i in R:
		num = 0
		done = false
		while num < len(Phis_list) and not done:
			Phis = Phis_list[num]
			if verbose:
				print("Working with twist ",D,", component ",i,"using",Phis.num_moments(),"moments")
				print("..computing p-adic L-function")
			L = Phis.pLfunction(r=i,quad_twist=D)
			if verbose:
				print("done!")
			d = S(L).degree()
			mu = mu_inv(L.substitute(w=1)/p^Phis.valuation(),p)
			if mu > 0:
				print("mu positive so probably not enough accuracy")
				num += 1
			else:
				lam = lambda_inv(L.substitute(w=1)/p^Phis.valuation(),p)
				if lam == 0:
					print("--lambda = 0 so nothing to check here")
					print("PASSED")
					if write_to_file:
						o = open(filename,'a')
						o.write(str(D)+": component="+str(i)+" PASSED (nothing to check)\n")
						o.close()
					done = true
				elif lam == 1 and 2*i % (p-1) == comp:
					print("--lambda = 1 because of FE so nothing to check here")
					print("PASSED")
					if write_to_file:
						o = open(filename,'a')
						o.write(str(D)+": component="+str(i)+" PASSED (nothing to check)\n")
						o.close()
					done = true
				else:
					print("--lambda =",lam)
					done = false
					giving_up = false
					## removes 1 from the lambda-invariant when sign of FE = -1
					if 2*i % (p-1) == comp and lam % 2 == 1:
						lam_mod = lam - 1
					else:
						lam_mod = lam 

					n = floor(log(lam_mod * p/(p-1))/log(p))
					if full_search:
						js = range(1,p)
					else:
						js = [1]
					while not done and not giving_up:
						print("----------------------------------------")				
						print("working at level",p^n)
						K = CyclotomicField(p^n)
						v = vp.extensions(K)[0]
						z = K.gen()
						vals = []
						if n > 0:
							error_bound = d / (p^(n-1)*(p-1))
						else:
							error_bound = d
						print("Error bound is",error_bound)
						ji = 0
						while ((full_search and ji < p-1) or (not full_search and ji == 0)) and not done:
							j = js[ji]
							error_bound_violated = false
							t_error = false
							a = 0
							while a < p^n and not error_bound_violated and not t_error:
								if not full_search:
									print("Using a =",a)
								else:
									print("Using (a,j)=",(a,j))
								t1 = S(L).substitute(T=z-1).substitute(w=(z^a-1+j*p)/p)
								val = v(t1) - Phis.valuation()
								if val >= error_bound:
									error_bound_violated = true
									bad_val = val
									print("***This is an error bound violation: value is",val,"error bound is",error_bound)
								elif 2*i % (p-1) == comp and lam % 2 == 1:
									extra_factor = v(z^a-z^2+p)
									print("---->Modification needed. original value:",val,"extra factor:",extra_factor)
									val = val - extra_factor
								print("-- Value has valuation",val)
								if val >=1:
									print("***Possible toroidal factor---aborting")
									t_error = true
								vals += [val]

								t2 = S(L).substitute(T=z^a-1).substitute(w=(z-1+j*p)/p)
								val = v(t2) - Phis.valuation()
								if val >= error_bound:
									error_bound_violated = true
									bad_val = val
									print("***This is an error bound violation: value is",val,"error bound is",error_bound)
								elif 2*i % (p-1) == comp and lam % 2 == 1:
									extra_factor = v(z-z^(2*a)+p)
									print("---->Modification needed. original value:",val,"extra factor:",extra_factor)
									val = val - extra_factor
								print("-- Value has valuation",val)
								if val >=1:
									print("***Possible toroidal factor---aborting")
									t_error = true
								vals += [val]
								print("Current list of values:",vals)
								a = a + 1
							m = max(vals)
							if error_bound_violated:
								print("Failed: not enough accuracy.  There was a valuation of",bad_val,"and the error bound is",error_bound)
							elif m < toroidal_bound:
								print("Passed! Max valuation is",m,"< 1")
								print("PASSED")
								if write_to_file:
									o = open(filename,'a')
									o.write(str(D)+": component="+str(i)+" PASSED (lambda ="+str(lam)+")\n")
									o.close()	

								done = true
							else:
								print("Failed: max val too high.  Max valuation is",m,">= 1")
							ji += 1
						n += 1
						if n > 3:
							num += 1
							if num < len(Phis_list) and not done:
								print("Going to more accurate family.")
								print("------------------------------------------------------------------------")
							elif not done:
								print("giving up!")
							giving_up = true
		if not done:
			print("*************************FAILED!!!***************************")
			if write_to_file:
				o = open(filename,'a')
				o.write(str(D)+": component="+str(i)+" FAILED\n")
				o.close()
		print("")
	return "done"

def run_me(filename,minD,maxD,Phis_list,level,step=1,log=true,full_search=true):
	old_stdout = sys.stdout

	for d in range(minD,maxD,step):
		print(d)
		if is_fundamental_discriminant(d) and gcd(d,level)==1:
			print("Working on twist d=",d)
			if log:
				log_file = open(filename+".details","a")
				sys.stdout = log_file			
			analyze_pLs(d,Phis_list,full_search=full_search,write_to_file=True,filename=filename)
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
