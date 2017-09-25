G1=DirichletGroup(5^2)  
chi1=(G1.0)^4
M1=ModularSymbols(chi1,2,1).cuspidal_subspace()
f1=M1.hecke_polynomial(5)
Qmu1=f1.parent().base_ring()
Qp=pAdicField(5,5)
Qpmu1.<z1>=Qp.ext((x+1)^4+(x+1)^3+(x+1)^2+(x+1)+1)
a51=f1.roots()[0][0]
i1=Qmu1.hom([1+z1],Qpmu1)
phi1=form_modsym_from_decomposition(M1)
phi1_padic = phi1.map2(i1)
print "Checking",phi1_padic.hecke(5)-phi1_padic.scale(i1(a51))
L2 = phi1_padic.MazurTate(5,2,10)
print str(L2.newton_slopes())

G2=DirichletGroup(5^3)  
chi2=(G2.0)^4
M2=ModularSymbols(chi2,2,1).cuspidal_subspace()
f2=M2.hecke_polynomial(5)
Qmu2=f2.parent().base_ring()
Qp=pAdicField(5,5)
Qpmu2.<z2>=Qp.ext((x+1)^20+(x+1)^15+(x+1)^10+(x+1)^5+1)
phi2=form_modsym_from_decomposition(M2)

G3=DirichletGroup(2*3^2)  
chi3=(G3.0)^2
M3=ModularSymbols(chi3,2,1).cuspidal_subspace()
f3=M3.hecke_polynomial(3)
Qmu3=f3.parent().base_ring()
Qp=pAdicField(3,5)
Qpmu3.<z3>=Qp.ext((x+1)^2+(x+1)+1)
a53=f3.roots()[0][0]
i3=Qmu3.hom([-1-z3],Qpmu3)
phi3=form_modsym_from_decomposition(M3)
phi3_padic = phi3.map2(i3)
print "Checking",phi3_padic.hecke(3)-phi3_padic.scale(i3(a53))
L2 = phi3_padic.MazurTate(3,2,10)
print str(L2.newton_slopes())
L3 = phi3_padic.MazurTate(3,3,10)
print str(L3.newton_slopes())
L4 = phi3_padic.MazurTate(3,4,10)
print str(L4.newton_slopes())

G4=DirichletGroup(2*3^3)  
chi4=(G4.0)^2
M4=ModularSymbols(chi4,2,1).cuspidal_subspace()
A4=M4.decomposition()[1]
f4=A4.hecke_polynomial(3)
Qmu4=f4.parent().base_ring()
Qp=pAdicField(3,5)
Qpmu4.<z4>=Qp.ext((x+1)^6+(x+1)^3+1)
i4=Qmu4.hom([-1-z4],Qpmu4)
phi4=form_modsym_from_decomposition(A4)
K4=phi4.base_ring()
k4=K4.defining_polynomial()
L4.<a>=Qmu4.extension(f4)
S.<y>=PolynomialRing(L4)
r=S(k4).roots()[0][0]
j4=K4.hom([r],L4)



phi3_padic = phi3.map2(i3)
print "Checking",phi3_padic.hecke(3)-phi3_padic.scale(i3(a53))
L2 = phi3_padic.MazurTate(3,2,10)
print str(L2.newton_slopes())
L3 = phi3_padic.MazurTate(3,3,10)
print str(L3.newton_slopes())
L4 = phi3_padic.MazurTate(3,4,10)
print str(L4.newton_slopes())


G1=DirichletGroup(5^2)  
chi1=(G1.0)^4
M1=ModularSymbols(chi1,4,1).cuspidal_subspace()
f1=M1.hecke_polynomial(5)
phi1=form_modsym_from_decomposition(M1)



G=DirichletGroup(11^2)  
chi=(G.0)^10
M=ModularSymbols(chi,2,1).cuspidal_subspace()
f=M.hecke_polynomial(11)
phi=form_modsym_from_decomposition(M)
Qmu=f.parent().base_ring()
Qp=pAdicField(11,10)
Qpmu1.<z1>=Qp.ext((x+1)^4+(x+1)^3+(x+1)^2+(x+1)+1)
a51=f1.roots()[0][0]
i1=Qmu1.hom([1+z1],Qpmu1)
phi1_padic = phi1.map2(i1)
print "Checking",phi1_padic.hecke(5)-phi1_padic.scale(i1(a51))
L2 = phi1_padic.MazurTate(5,2,10)
print str(L2.newton_slopes())

p=3
r=2
N=1
k=4
G=DirichletGroup(N*p^r)
chi=(G.0)^(p-1)
print "Forming modsym space"
M=ModularSymbols(chi,k,1).cuspidal_subspace()
print M
print "Hecke poly"
f=M.hecke_polynomial(p)

print "Decomposing"
A=M.decomposition()
print A
for b in range(len(A)):
	B=A[b]
	print "forming modular symbol"
	phi=form_modsym_from_decomposition(B)
	K=phi.base_ring()
	print "Finding primes over p in",K
	pp=K.primes_above(p)
	bool,ap=phi.is_Tq_eigen(p)
	slopes = [ap.valuation(pp[j])/pp[j].absolute_ramification_index() for j in range(len(pp))]
	print "Slopes:",slopes
	h = max(slopes)
	print "Max slope is",h
	print "Computing MT"
	L=phi.MazurTate(p,4,10,h=h,verbose=false)
	v=L.padded_list()
	for j in range(len(pp)):
		print "Slope",ap.valuation(pp[j])/pp[j].absolute_ramification_index()
		NP=NewtonPolygon([(a,v[a].valuation(pp[j])/pp[j].absolute_ramification_index()) for a in range(len(v))]).slopes()
		ans=[]
		sNP = Sequence(Set(NP))
		sNP.sort()
		print [[NP.count(sNP[a]),-sNP[a]] for a in range(len(sNP))]

