def ram(w):
	pi = w.uniformizer()
	return 1/w(pi)

def inertia(w):
	kw = w.residue_field()
	return kw.degree()

def ev(A,ell):
	M = A.ambient_module()
	vec = A.dual_eigenvector()
	v = M.hecke_matrix(ell)*vec

	a = 0
	while vec[a] == 0:
		a += 1
	return v[a]/vec[a]

#max: reduces a_q for q <= max
def reduce(A,p,max=30,include_p=false):
	N = A.level()
	vec = A.dual_eigenvector()
	K = vec.parent().base_field()
	v = QQ.valuation(p)
	ws = v.extensions(K)

	evs = {}
	for q in primes(max+1):
		if N*p % q != 0 or (q == p and include_p):		
			evs[q] = ev(A,q)

	rhobars = []
	for w in ws:
		evs_modp = {}
		for q in primes(max+1):
			if N*p % q != 0 or (q == p and include_p):		
				evs_modp[q] = w.reduce(evs[q])
#		rhobars += [(evs_modp,ram(w),inertia(w))]
		rhobars += [evs_modp]

		## adds galois conjugates of reduction
		f = w.residue_field().degree()
		for a in range(1,f):
			evs_frob = {}
			for q in primes(max+1):
				if N*p % q != 0:		
					evs_frob[q] = evs_modp[q]^(p^a)
#			rhobars += [(evs_frob,ram(w),inertia(w))]
			rhobars += [evs_frob]

	return rhobars

def equal_on_common_keys(dict1, dict2):
	common_keys = set(dict1.keys()) & set(dict2.keys())

	for key in common_keys:
		if dict1[key] != dict2[key]:
			return False

	# If the loop completes without returning False, dictionaries are equal on common keys
	return True

#def contains(list_of_dicts,dict):
#	for d in list_of_dicts:
#		if equal_on_common_keys(d,dict):
def is_dict_in_list(dict_list, target_dict):
    for current_dict in dict_list:
        common_keys = set(current_dict.keys()) & set(target_dict.keys())

        # Check if the dictionaries are equal on common keys
        if all(current_dict[key] == target_dict[key] for key in common_keys):
            return True

    # If no match is found in the loop, the dictionary is not in the list
    return False


def is_level_lowerable(A,w,max):
	level_lowering = {}
	p = w.p()
	N = A.level()
	k = A.weight()
	evs_modp = {}
	for q in primes(max+1):
		if N*p % q != 0:
			evs_modp[q] = w.reduce(ev(A,q))

	bool = false
	for d in N.divisors():
		if d < N:
			if not p in level_lowering.keys():
				level_lowering[p] = {}
			if not d in level_lowering[p].keys():
				level_lowering[p][d] = []
				M = ModularSymbols(d,k,1).cuspidal_subspace().new_subspace()
				As = M.decomposition()
				for A in As:
					level_lowering[p][d] += reduce(A,p)
			bool = bool or is_dict_in_list(level_lowering[p][d],evs_modp)
			if bool:
				return bool

	return bool


def is_level_lowerable_quad_twist(A,w,max,level_lowering={}):
	p = w.p()
	N = A.level()
	assert N % p == 0,"Level not divisible by p"
	k = A.weight()
	evs_modp = {}
	G = DirichletGroup(p,QQ)
	chip = G.0
	for q in primes(max+1):
		if N*p % q != 0:
			evs_modp[q] = w.reduce(ev(A,q))*GF(p)(chip(q))

	bool = false
	for d in N.divisors():
		if d < N:
			if not p in level_lowering.keys():
				level_lowering[p] = {}
			if not d in level_lowering[p].keys():
				level_lowering[p][d] = []
				M = ModularSymbols(d,k,1).cuspidal_subspace().new_subspace()
				As = M.decomposition()
				for A in As:
					level_lowering[p][d] += reduce(A,p)
			bool = bool or is_dict_in_list(level_lowering[p][d],evs_modp)
			if bool:
				return bool

	return bool






