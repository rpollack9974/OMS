class modsym_dist_fam(modsym):
	def ms(self):
		"""demotes to a regular modular symbol"""
		return modsym(self.level,self.data,self.manin)

	def num_moments(self):
		return self.data[0].num_moments()

	def p(self):
		return self.data[0].p

	def deg(self):
		return self.data[0].deg

	## This function returns a number between 0 and p-2 which represents which 	
	## disc in weight the family is supported on -- the integer i corresponds to the
	## disc with characters with finite order part omega^i
	def disc(self):
		return self.data[0].disc()

	def change_deg(self,new_deg):
		v=[self.data[r].change_deg(new_deg) for r in range(self.ngens())]
		return modsym_dist_fam(self.level,v,self.manin)

	def specialize(self,k):
		v=[]
		for j in range(0,len(self.data)):
			v=v+[self.data[j].specialize(k)]
		return modsym_dist_aws(self.level,v,self.manin)
	
	def valuation(self):
		#print [self.data[j].valuation() for j in range(len(self.data))]
		return min([self.data[j].valuation() for j in range(len(self.data))])

	def normalize(self):
		v=[]
		for j in range(0,len(self.data)):
			v=v+[self.data[j].normalize()]
		return modsym_dist_fam(self.level,v,self.manin)
	
	def change_precision(self,M):
		v=[]
		for j in range(0,len(self.data)):
			v=v+[self.data[j].change_precision(M)]
		return modsym_dist_fam(self.level,v,self.manin)

	## This procedure tries to find a power series c(w) such that 
	##      self | T_q = c(w) self
	## It prints to the screen the self | T_q - c(w) self (so all displayed zeroes
	## means it worked, and it returns the power series c(w)
	def is_Tq_eigen(self,q):
		p = self.p()
		M = self.num_moments()
		R = self.data[0].moment(0).parent()
		T = PowerSeriesRing(QQ,'y')

		Phiq = self.hecke(q)
		aq = R(T(T(Phiq.data[0].moment(0))/T(self.data[0].moment(0))).padded_list())

		print (self.scale(aq) - Phiq).normalize()

		return ps_normalize(aq,p,M-self.valuation())

	#@cached_function
	def phi_on_Da(self,a,D):
		p=self.p()
		ans=self.zero_elt()
		for b in range(1,abs(D)+1):
			if gcd(b,D)==1:	
				ans=ans+self.eval(Matrix(2,2,[1,b/abs(D),0,1])*Matrix(2,2,[a,1,p,0])).act_right(Matrix(2,2,[1,b/abs(D),0,1])).scale(kronecker(D,b)).normalize()
		return ans.normalize()

	#@cached_function
	def basic_integral(self,a,j,ap,D):
		"""returns int_{a+pZ_p} (z-{a})^j dPhi(0-infty) -- see formula [PS, sec 9.2] """
		M=self.num_moments()
		p=self.p()
		ap=ap*kronecker(D,p)
		ans=0
		for r in range(j+1):
			ans=ans+binomial(j,r)*(a-teich(a,p,M))^(j-r)*p^r*self.phi_on_Da(a,D).moment(r)
		return ans/ap

	def pLfunction_coef(self,ap,n,r,D,gam,base_ring=QQ,error=None):
		"""Returns the n-th coefficient of the p-adic L-function in the T-variable of a quadratic twist of self.  If error is not specified, then the correct error bound is computed and the answer is return modulo that accuracy.

	Inputs:
		self -- overconvergent Hecke-eigensymbol;
		ap -- eigenvalue of U_p;
		n -- index of desired coefficient
		r -- twist by omega^r
		D -- discriminant of quadratic twist
		gam -- topological generator of 1+pZ_p"""


		S.<z>=PolynomialRing(QQ)
		p=self.p()
		M=self.num_moments()
		R=pAdicField(p,M)
		lb=loggam_binom(p,gam,S.0,n,2*M)
		dn=0
		if n==0:
			err=M
		else:
			if (error==None):
				err=min([j+lb[j].valuation(p) for j in range(M,len(lb))])
			else:
				err=error
			lb=[lb[a] for a in range(M)]
		print "err=",err
		for j in range(min(M,len(lb))):
			cjn=lb[j]
			temp=0
			for a in range(1,p):
				temp=temp+teich(a,p,M)^(r-j)*self.basic_integral(a,j,ap,D)
			dn=dn+cjn*temp
		if base_ring == QQ:
			return dn+O(p^err)
		else:
			v = list(dn)
			v = [v[a] + O(p^err) for a in range(len(v))]
			print v
			t = 0*w
			for j in range(len(v)):
				t += v[j] * w^j
			dn = t
			print dn
			return dn


	def pLfunction(self,ap,r=0,max=Infinity,base_ring=QQ,quad_twist=None):
		"""Returns the p-adic L-function in the T-variable of a quadratic twist of self

	Inputs:
		self -- overconvergent Hecke-eigensymbol;
		ap -- eigenvalue at p;
		r -- twist by omega^r
		quad_twist -- conductor of quadratic character"""

		if quad_twist==None:
			D=1
		else:
			D=quad_twist
		M=self.num_moments()
		p=self.p()
		gam=1+p
	#	for a in range(1,p):
	#		for j in range(M):
	#			basic_integral(self,a,j,ap,D)

		SS.<T>=PowerSeriesRing(base_ring)
		ans=self.pLfunction_coef(ap,0,r,D,gam,base_ring=base_ring)+0*T
		print ans
		S.<z>=PolynomialRing(QQ)
		err=Infinity
		n=1
		while (err>0) and (n<min(M,max)):
			print n
			lb=loggam_binom(p,gam,z,n,2*M)
			err=min([j+lb[j].valuation(p) for j in range(M,len(lb))])
			print "err2",err
			if err>0:
				dn=self.pLfunction_coef(ap,n,r,D,gam,base_ring=base_ring,error=err)
				print n,dn
				print
				ans=ans+dn*T^n
				print "ans",ans
				print
				print

			n=n+1

		return ans


	
#######################################################################################################################
##  This function produces a random family of OMSs.
##
##  p -- prime
##  N -- tame level
##  char -- character of conductor dividing N (or maybe Np?)
##  M -- number of moments
##  r -- integer between 0 and p-2 indicating which disc on weight space we are working
##  w -- variable which is just carried around because I had trouble with symbols being defined over different rings
#####################################################################################################################
def random_OMS_fam(p,N,char,M,r,w):
	deg = M
	v = []

	manin = manin_relations(N*p)

	## this loop runs through each generator (different from D_infty) and picks a random value for that generator
	for j in range(1,len(manin.gens)):
		rj = manin.gens[j]
		mus = random_dist_fam(p,M,deg,r,char,w)
		if (manin.twotor.count(rj) == 0) and (manin.threetor.count(rj) == 0):
			v = v + [mus]
		elif (manin.twotor.count(rj) != 0):
			## Case of two torsion (See [PS] section 4.1)
			gam = manin.gen_rel_mat(j)
			v = v + [(mus - mus.act_right(gam)).scale(1/2)]
		else:
			## Case of three torsion (See [PS] section 4.1)	
			gam = manin.gen_rel_mat(j)
			v = v + [(mus.scale(2) - mus.act_right(gam) - mus.act_right(gam^2)).scale(1/3)]

	t = v[0].zero()
	## This loops adds up around the boundary of fundamental domain except the two verticle lines
	for j in range(1,len(manin.gens)):
		rj = manin.gens[j]
		if (manin.twotor.count(rj) == 0) and (manin.threetor.count(rj) == 0):
			t = t + v[j-1].act_right(manin.gen_rel_mat(j)) - v[j-1]
		else:
			t = t - v[j-1]


	## t now should be sum Phi(D_i) | (gamma_i - 1) - sum Phi(D'_i) - sum Phi(D''_i)
	## (Here I'm using the opposite sign convention of [PS1] regarding D'_i and D''_i)

	## We now need to make some adjustment of Phi(D_i) to make t have total measure 0

	j = 1
	rj = manin.gens[j]
	while (j < len(manin.gens)-1) and ((manin.twotor.count(rj) != 0) or (manin.threetor.count(rj) != 0)):
		j = j + 1
		rj = manin.gens[j]
	assert j < len(manin.gens) - 1, "Everything is 2 or 3 torsion!  NOT YET IMPLEMENTED IN THIS CASE"

	gam = manin.gen_rel_mat(j)
	a = gam[0,0]
	c = gam[1,0]
	K = aut(p,deg,M,a,c,r,char,w)
	K0 = K[0]  ## K0 is the coefficient of z^0 in K
	K1 = K[1]  ## K1 is the coefficient of z^1 in K
	t0 = t.moment(0)
	T = PowerSeriesRing(QQ,'ww')
	err = mus.zero()

	k = r
	if r != 0:
		## The following code simply computes -t0/(K0-1)
		temp = T(T(-t0)/T(K[0]-1))
		temp = temp.truncate(deg)
		R = w.parent()
		temp = R(temp.padded_list())

		print "The result will be a lifting modulo p^",valuation(temp.substitute(w=((1+p)^k-1)/p),p)
		err.moments[0] = temp
	else:
		## The following code simply computes -t0/(K1)
		temp=T(T(-t0)/T(K1))
		temp=temp.truncate(deg)
		R = w.parent()
		temp=R(temp.padded_list())
		print "The result will be a lifting modulo p^",valuation(temp.substitute(w=((1+p)^k-1)/p),p)
		err.moments[1] = temp

	v[j-1] = v[j-1] + err
	t = t + err.act_right(gam)-err
	print "Checking that t has total measure 0: ",t.normalize().moment(0)

	mus = t.solve_diff_eqn()
	print (mus.act_right(Matrix(2,2,[1,1,0,1]))-mus-t).normalize()

	v = [mus.scale(-1)] + v

	Phis = modsym_dist_fam(N*p,v,manin)

	return Phis


@cached_function
def teich(a,p,M):
	R=pAdicField(p,M)
	return ZZ(R.teichmuller(a))

#@cached_function
def logp(p,z,M):
	"""returns the truncation sum_{j=1}^{M-1} (-1)^(j+1)/j z^j of log_p(1+z)"""
	ans=0
	for j in range(1,M):
		ans=ans+(-1)^(j+1)/j*z^j
	return ans

@cached_function
def loggam_binom(p,gam,z,n,M):
	L=logp(p,z,M)
	logpgam=L.substitute(z=(gam-1))
	loggam=L/logpgam
#	print loggam
#	print binomial(loggam,n)

	return binomial(loggam,n).list()
