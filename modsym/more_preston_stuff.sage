G=DirichletGroup(3,ZZ)
chi=G.0
p=331
k=149
acc=5
print "Forming modular symbol space"
M=ModularSymbols(chi,k,1).cuspidal_subspace(); M
print "Forming modular symbols"
wp=load("dual_vector_padic.sobj")
phi=form_modsym_from_decomposition_padic_enhanced(M,p,acc+5,wp)
phi=phi.scale(p)
ap=1 + 303*331^2 + 254*331^3 + 278*331^4 + 94*331^5 + 29*331^6 + 325*331^7 + 294*331^8 + 107*331^9
convert_modsym_padic_to_integral(phi)
print "p-stabilizing"
phip=phi.p_stabilize_ordinary(p,ap,acc+5)
print "lifting"
Phi=phip.lift_to_OMS(p,acc)
 #Phi=Phi.hecke(2)-Phi.scale(1+chi(2)*2^(k-1))
Phi=Phi.hecke(2)-Phi.scale(chi(2)+2^(k-1))
t=cputime()
for j in range(acc+1):
 	print j,cputime(t)
 	Phi=Phi.hecke(p)


## (p,k) is an irregular pair and this function returns the OMS
## of a cuspform congruent to E_k.  Probably assuming rank 1
#def Phi_congruent_to_eisen(p,k):
	














