load("master.sage")

N=17
k=4
M = ModularSymbols(N,k,1).cuspidal_subspace()
Ds = M.decomposition()
## randomly picking one of the modular forms of this weight and level
A = Ds[0]
phi = form_modsym_from_decomposition(A)
acc = 10  ## accuracy to compute Teichmuller lifts
p = 3
n = 3
print "The Mazur-Tate element for p=",p,"and level",n,"is:"
print phi.MazurTate(p,n,acc)