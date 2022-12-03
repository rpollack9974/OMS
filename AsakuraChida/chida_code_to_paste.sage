#---------------------
#determine E
#determine p
#determine M

Phi,A,q = find_slope1_kernel_relation(E,p,M)
e = matrix_valuation(A,p)
B,C,r = row_reduce_mod_pn2(A/p^e,p,1)



#using row_reduce_mod_pn2 with matrix B and row r to get polynoial relation f

print "Killing off kernel of specialization"
f = extract_poly_relation(C,r)
beta = beta_root(E,p,M)
Phi_eigen = find_eigen_slope1(Phi,f,q,E)
did_it_work(Phi_eigen)
Phi_eigen = renormalize(Phi_eigen,E)
print "And here's the value..."
val = value_at_0(Phi_eigen,p,beta) 
v = val.valuation(p)
if v < 0:
	val_norm = ((val * p^(-v)) % (p^(Phi_eigen.num_moments()))) * p^v
else:
	val_norm = val % (p^(Phi_eigen.num_moments()))
val_norm