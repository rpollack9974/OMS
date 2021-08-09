print "Forming modular symbol space of level 67 and weight 2"

M=ModularSymbols(67,2,1)
S = M.cuspidal_submodule()
N = S.new_subspace()

print "Decomposing it"
D = N.decomposition(); D

print "here is the decomposition"
print D

print "Only the 2nd and 3rd are defined over a quadratic extension"

print "Let's start with the 3rd one:"

phi=form_modsym_from_decomposition(D[2])

print "Here's it's value on {0}-{infty}"
print phi.eval(Id)

print "So L(f,1) doesn't vanish for this form."

print "Let's try the 2nd one:"

phi=form_modsym_from_decomposition(D[1])

print "Here's it's value on {0}-{infty}"
print phi.eval(Id)

print "So L(f,1) vanishes for this form -- this is the one we want."

print "This is how the symbol is stored -- by its values on a finite list of divisors"

print phi

p=11
print "Taking p=",p

print "Embedding the modular symbol into Q_p in the two possible ways"
phis=phi.coerce_to_Qp(p,10)

phi1=phis[0][0]
phi2=phis[1][0]

print "Here is our first Q_p symbol:"

print phi1

print "Here is our second Q_p symbol:"

print phi2

psi1=phis[0][1]
psi2=phis[1][1]

print "They were obtained by applying the maps ",psi1,psi2

print "p-stabilizing to level 67*p (modulo p^10)"

f=D[1].q_eigenform(20,'alpha')
ap=f[p]

phi1p = phi1.p_stabilize_ordinary(p,ZZ(psi1(ap)),10)
phi2p = phi2.p_stabilize_ordinary(p,ZZ(psi2(ap)),10)

print "Lifting first symbol"

Phi1=phi1p.lift_to_OMS_eigen(p,10)

print "Lifting second symbol"

Phi2=phi2p.lift_to_OMS_eigen(p,10)

print "Now we compute the special value of the first symbol"

print "We first need the unit root of x^2-a_p*x+p in Q_p (via the two embdedings of a_p"

h1=x^2-ZZ(psi1(ap))*x+p
alpha1=h1.roots()[0][0]
h2=x^2-ZZ(psi2(ap))*x+p
alpha2=h2.roots()[0][0]

v1= pLfunction_coef(Phi1,alpha1,1,1,20)
print v1

print "Here's the special value of the first symbol"

v2= pLfunction_coef(Phi2,alpha2,1,1,20)
print v2

print "Here's their product"
print v1*v2
