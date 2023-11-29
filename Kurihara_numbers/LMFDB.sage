# LMFDB LABELS

def trace(M,prec):
	f = M.q_eigenform(prec+1,'a')	
	return [f[i].trace() for i in range(1,prec+1)]


def is_all_distinct(lst):
    seen = set()

    for elem in lst:
        # Convert elements to tuples for nested lists
        elem_tuple = tuple(elem) if isinstance(elem, list) else elem

        if elem_tuple in seen:
            return False
        seen.add(elem_tuple)

    return True



alph = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"]

#returns the list of digits of k base n
def digits_base_n(k, n):
	v = k.digits(n)
	v.reverse()
	return v 

def number_to_letters(k, alph):
    v = digits_base_n(k, 26)

    # Replace all n0... occurrences with (n-1)26 --- this is a Z and 0 problem
    for i in range(len(v) - 1, 0, -1):
        if v[i] == 0:
            v[i] += 26
            v[i - 1] -= 1

    # Remove leading zeros
    while v and v[0] == 0:
        v = v[1:]

    ans = ''.join(alph[val - 1] for val in v)

    return ans

#takes in a complete list of newforms of a fixed weight and level
#and returns a list of their LMFDB labels
def LMFDB_labels_fixed_level(lst):
	if not lst:
		return []

	labels = []
	prec = 1
	done = False

	while not done:
		traces = [trace(f, prec) for f in lst]
		if is_all_distinct(traces):
			done = True
		else:
			prec += 1

	sorted_traces = sorted(traces)

	N = lst[0].level()
	k = lst[0].weight()
	labels = [str(N) + "." + str(k) + ".a." + number_to_letters(sorted_traces.index(t) + 1, alph) for t in traces]

	return labels

