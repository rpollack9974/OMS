load("master.sage")
load("ordinary_subspace/ordinary_subspace.sage")

p = 3
N = 23
G = DirichletGroup(N*p,QQ)
chi = G.1 #this character should have conductor N and level Np
k = -1 #this corresponds to weight 1
acc = 8 #this makes the programs work mod p^acc
sign = 1
d = dimension_of_ordinary_subspace(chi,p,k,sign)
Phis = form_basis2(N,p,k,acc,chi,d,sign)

aqs = {2:-1,3:-1,5:0,7:0,11:0,13:-1}
for q in primes(15):
	T = hecke_matrix(Phis,q) 
	f = charpoly(T)
	print "Char poly of T_",q," is: ",f
	print "a_",q," of classical weight 1 form is: ",aqs[q]
	print "Char poly at a_",q," has valuation: ",f.substitute(aqs[q]).valuation(p)
	print ""
