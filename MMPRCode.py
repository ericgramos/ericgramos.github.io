#Defines ring of symmetric functions over Q[q]
base_ring.<q> = PolynomialRing(QQ)
Sym = SymmetricFunctions(base_ring)
Sym.inject_shorthands()

def lieN(n):
	"""
	Return the symmetric function `\\ell_n`, which is the
	Frobenius image of `\\mathrm{Lie}_{(n)}` in [HR17].
	"""
	div = divisors(n)
	return (1/n)*sum(moebius(d)*p[d]^(n//d) for d in div)

def lieLambda(lamb):
	"""
	Return the symmetric function 
	`\\mathrm{ch}(\\mathrm{Lie}_{\\lambda})` from [HR17].
	"""
	expPartition = lamb.to_exp()
	return prod(h[expPartition[i]](lieN(i+1)) for i in range(0,len(expPartition)))

def wLambda(lamb):
	"""
	Return the symmetric function 
	`\\mathrm{ch}(\\mathrm{W}_{\\lambda})` from [HR17].
	"""
	expPartition = lamb.to_exp()
	return prod(h[expPartition[i]](lieN(i+1).omega()) for i in range(0,len(expPartition))[::2])*prod(e[expPartition[i]](lieN(i+1).omega()) for i in range(1,len(expPartition))[::2])

def rankOfPartition(lamb):
	"""
	Return the rank of a partition `\\lambda` as in [HR17, Definition 2.5].
	"""
	expPartition = lamb.to_exp()
	return sum(i*expPartition[i] for i in range(0,len(expPartition)))

def rep_dimension(f):
	"""
	Return the dimension of the representation
	associated to the symmetric function `f`.
	"""
	f = s(f)
	d = f.monomial_coefficients()
	result = 0
	for mu in d:
		result += Partition(mu).dimension() * d[mu]
	return result

#note to myself
#C_n^i = theConf(n,3,i)
#A_n^i = theConf(n,2,i)

memoize_A = {}
def A(n,i):
	"""
	Return the `i`-th cohomology of the configuration space FINISH
	"""
	if (n,i) not in memoize_A:
		ourPartitions = []
		for k in Partitions(n):
			if rankOfPartition(k) == i:
				ourPartitions.append(k)
		memoize_A[(n,i)] = sum(wLambda(lamb) for lamb in ourPartitions)
	return memoize_A[(n,i)]

def qA(n):
	return sum(A(n,j) * q**j for j in range(n))

memoize_B = {}
def B(n,i):
	"""
	Return the symmetric function corresponding to the 
	`\\mathfrak{S}_n`-representation `B_n^i` from [MMPR21].
	"""
	if (n,i) not in memoize_B:
		if n == 1:
			print('do not give me n = 1')
		if i == 0:
			memoize_B[(n,i)] = s([n])
		else:
			memoize_B[(n,i)] = A(n,i)-B(n,i-1)
	return memoize_B[(n,i)]

def qB(n):
	"""
	Return the symmetric function, with `q`-coefficients,
	corresponding to the graded
	`\\mathfrak{S}_n`-representation `B_n` from [MMPR21].
	"""
	return sum(B(n,j) * q**j for j in range(n))

memoize_C = {}
def C(n,i):
	"""
	Return the `i`-th cohomology of the configuration space FINISH
	"""
	if (n,i) not in memoize_C:
		ourPartitions = []
		for k in Partitions(n):
			if rankOfPartition(k) == i:
				ourPartitions.append(k)
		memoize_C[(n,i)] = sum(lieLambda(lamb) for lamb in ourPartitions)
	return memoize_C[(n,i)]

'''
memoize_theConf = {}
def theConf(n,d,i):
	"""
	Return the `i`-th cohomology of the configuration space FINISH
	"""
	if (n,d,i) not in memoize_theConf:
		ourPartitions = []
		for k in Partitions(n):
			if rankOfPartition(k) == i:
				ourPartitions.append(k)
		if (d % 2) == 0:
			memoize_theConf[(n,d,i)] = sum(wLambda(lamb) for lamb in ourPartitions)
		else:
			memoize_theConf[(n,d,i)] = sum(lieLambda(lamb) for lamb in ourPartitions)
	return memoize_theConf[(n,d,i)]

def qtheConf(n,d):
	return sum(theConf(n,d,j) * q**j for j in range(n))
'''

def qC(n):
	return sum(C(n,j) * q**j for j in range(n))

memoize_D = {}
def D(n,i):
	"""
	Return the symmetric function corresponding to the 
	`\\mathfrak{S}_n`-representation `D_n^i` from [MMPR21].
	"""
	if (n,i) not in memoize_D:
		if n == 1:
			print('do not give me n = 1')
		if i == 0:
			memoize_D[(n,i)] = s([n])
		else:
			memoize_D[(n,i)] = theConf(n,3,i)-(D(n,i-1).itensor(s([n-1,1])))
	return memoize_D[(n,i)]

def qD(n):
	"""
	Return the symmetric function, with `q`-coefficients,
	corresponding to the graded
	`\\mathfrak{S}_n`-representation `D_n` from [MMPR21].
	"""
	return sum(D(n,j) * q**j for j in range(n))

def ceiling_division(n, d):
	"""
	Return the ceiling of `n` divided by `d`.
	"""
    return -(n // -d)


'''
#maybe delete some of this stuff and put it in another file for reference

def getSubLists(n):
	listOfNumbers = []
	for i in range(1,ceiling_division(n,2)):
		for j in range(i+1,ceiling_division(n,2)+1):
			for k in range(j,n-j+i):
				listOfNumbers.append([i,j,k,j+k-i])
	return listOfNumbers

def eqLC(n,d,i,j,k,l):
	return (theConf(n,d,j).itensor(theConf(n,d,k))-theConf(n,d,i).itensor(theConf(n,d,l))).is_schur_positive()

def eqLC_all(n,d):
	ourSubLists = getSubLists(n)
	bool = True
	for aList in ourSubLists:
		i = aList[0]
		j = aList[1]
		k = aList[2]
		l = aList[3]
		if (eqLC(n,d,i,j,k,l) == False):
			print('Counterexample:','n=',n,' d=',d,' i=',i,' j=',j,' k=',k,' l=',l)
			bool = False
	if (bool == True):
		print('All good!')
	else:
		print('We failed, man...')

def ACheck(n,d,i):
	R.<t> = QQ[]
	f = prod((1+m*t) for m in range(1,n))
	return (rep_dimension(theConf(n,d,i)) == f.list()[i])

def BCheck(n,i):
	R.<t> = QQ[]
	f = prod((1+m*t) for m in range(2,n))
	return (rep_dimension(B(n,i)) == f.list()[i])

def DCheck(n,i):
	R.<t> = QQ[]
	f = prod((1+m*t) for m in range(1,n-1))
	return (rep_dimension(D(n,i)) == f.list()[i])

def stabilize(f):
	f = s(f)
	addOneBox = lambda mu: mu.add_cell(0)
	return f.map_support(addOneBox)

def whenStab(d,i,j):
	for n in range(3*(i+j),4*(i+j)+1):
		if (theConf(n,d,i)) != 0 and (theConf(n,d,j)) != 0 and (theConf(n-1,d,i)) != 0 and (theConf(n-1,d,j)) != 0:
			if s(theConf(n,d,i).itensor(theConf(n,d,j))) == stabilize(theConf(n-1,d,i).itensor(theConf(n-1,d,j))):
				print(n-1,'works')
				if (n-1==3*(i+j) and d==3):
					print('and',n-1,'= 3(i+j)')
				if (n-1==3*(i+j)+2 and d==2):
					print('and',n-1,'= 3(i+j)+2')
				break

'''

memoize_Atensor = {}
def Atensor(n,m,i):
	"""
	Return the internal (Kronecker) product of the symmetric functions 
	`A_n^i` and `A_n^{m-i}`.
	"""
	if (n,m,i) not in memoize_Atensor:
		if (A(n,i)) != 0 and (A(n,m-i) != 0):
			memoize_Atensor[(n,m,i)] = A(n,i).itensor(A(n,m-i))
		else:
			memoize_Atensor[(n,m,i)] = s.zero()
	return memoize_Atensor[(n,m,i)]

memoize_Btensor = {}
def Btensor(n,m,i):
	"""
	Return the internal (Kronecker) product of the symmetric functions 
	`B_n^i` and `B_n^{m-i}`.
	"""
	if (n,m,i) not in memoize_Btensor:
		if (B(n,i)) != 0 and (B(n,m-i) != 0):
			memoize_Btensor[(n,m,i)] = B(n,i).itensor(B(n,m-i))
		else:
			memoize_Btensor[(n,m,i)] = s.zero()
	return memoize_Btensor[(n,m,i)]

memoize_Ctensor = {}
def Ctensor(n,m,i):
	"""
	Return the internal (Kronecker) product of the symmetric functions 
	`C_n^i` and `C_n^{m-i}`.
	"""
	if (n,m,i) not in memoize_Ctensor:
		if (C(n,i)) != 0 and (C(n,m-i) != 0):
			memoize_Ctensor[(n,m,i)] = C(n,i).itensor(C(n,m-i))
		else:
			memoize_Ctensor[(n,m,i)] = s.zero()
	return memoize_Ctensor[(n,m,i)]

memoize_Dtensor = {}
def Dtensor(n,m,i):
	"""
	Return the internal (Kronecker) product of the symmetric functions 
	`D_n^i` and `D_n^{m-i}`.
	"""
	if (n,m,i) not in memoize_Dtensor:
		if (D(n,i)) != 0 and (D(n,m-i) != 0):
			memoize_Dtensor[(n,m,i)] = D(n,i).itensor(D(n,m-i))
		else:
			memoize_Dtensor[(n,m,i)] = s.zero()
	return memoize_Dtensor[(n,m,i)]

def Achecker(n,m):
	"""
	Write win or lose into A_ELC.txt, depending on whether or not
	`A_n` is strongly equivariantly log concave in degree `m`.
	"""
	ourBool = True
	f = open('A_ELC.txt', 'a', buffering=1)
	ourTensors = []
	for i in range(1,m//2+1):
		ourTensors.append(Atensor(n,m,i))
	for i in range(0,len(ourTensors)-1):
		if (ourTensors[i+1]-ourTensors[i]).is_schur_positive() == False:
			ourBool = False
			break
	if (ourBool == False):
		print('lose',file=f)
	else:
		print('win',file=f)

def Bchecker(n,m):
	"""
	Write win or lose into B_ELC.txt, depending on whether or not
	`B_n` is strongly equivariantly log concave in degree `m`.
	"""
	ourBool = True
	f = open('B_ELC.txt', 'a', buffering=1)
	ourTensors = []
	for i in range(1,m//2+1):
		ourTensors.append(Btensor(n,m,i))
	for i in range(0,len(ourTensors)-1):
		if (ourTensors[i+1]-ourTensors[i]).is_schur_positive() == False:
			ourBool = False
			break
	if (ourBool == False):
		print('lose',file=f)
	else:
		print('win',file=f)

def Cchecker(n,m):
	"""
	Write win or lose into C_ELC.txt, depending on whether or not
	`C_n` is strongly equivariantly log concave in degree `m`.
	"""
	ourBool = True
	f = open('C_ELC.txt', 'a', buffering=1)
	ourTensors = []
	for i in range(1,m//2+1):
		ourTensors.append(Ctensor(n,m,i))
	for i in range(0,len(ourTensors)-1):
		if (ourTensors[i+1]-ourTensors[i]).is_schur_positive() == False:
			ourBool = False
			break
	if (ourBool == False):
		print('lose',file=f)
	else:
		print('win',file=f)

def Dchecker(n,m):
	"""
	Write win or lose into D_ELC.txt, depending on whether or not
	`D_n` is strongly equivariantly log concave in degree `m`.
	"""
	ourBool = True
	f = open('D_ELC.txt', 'a', buffering=1)
	ourTensors = []
	for i in range(1,m//2+1):
		ourTensors.append(Dtensor(n,m,i))
	for i in range(0,len(ourTensors)-1):
		if (ourTensors[i+1]-ourTensors[i]).is_schur_positive() == False:
			ourBool = False
			break
	if (ourBool == False):
		print('lose',file=f)
	else:
		print('win',file=f)

def doWeHaveAThm(m):
	"""
	Checks strong equivariant log concavity of `A_n` in degree `m`
	for all `n`, using representation stability.
	"""
	f = open('A_ELC.txt', 'a', buffering=1)
	count = 1
	for n in range(1,3*m+3):
		print('iteration',count,'of',3*m+2,'total iterations',file=f)
		count = count + 1
		Achecker(n,m)

def doWeHaveBThm(m):
	"""
	Checks strong equivariant log concavity of `B_n` in degree `m`
	for all `n`, using representation stability.
	"""
	f = open('B_ELC.txt', 'a', buffering=1)
	count = 2
	for n in range(2,3*m+3):
		print('iteration',count,'of',3*m+2,'total iterations',file=f)
		count = count + 1
		Bchecker(n,m)

def doWeHaveCThm(m):
	"""
	Checks strong equivariant log concavity of `C_n` in degree `m`
	for all `n`, using representation stability.
	"""
	f = open('C_ELC.txt', 'a', buffering=1)
	count = 1
	for n in range(1,3*m+1):
		print('iteration',count,'of',3*m,'total iterations',file=f)
		count = count + 1
		Cchecker(n,m)

def doWeHaveDThm(m):
	"""
	Checks strong equivariant log concavity of `D_n` in degree `m`
	for all `n`, using representation stability.
	"""
	f = open('D_ELC.txt', 'a', buffering=1)
	count = 2
	for n in range(2,3*m+1):
		print('iteration',count,'of',3*m,'total iterations',file=f)
		count = count + 1
		Dchecker(n,m)

def MequalsD(n):
	"""
	Returns `True` or `False` depending on whether or not
	the two symmetric functions `M_n` and `qD_n`, with `q`-coefficients,
	are equal.
	"""
	return (M(n) == qD(n))

def A_ELC(last):
	"""
	Checks strong equivariant log concavity of `A_n` in degrees `1`
	through `last` for all `n`, using representation stability.
	"""
	import time
	totalStart = time.time()
	f = open('A_ELC.txt', 'a', buffering=1)
	for counter in range(1,last+1):
		subStart = time.time()
		print('Checking m =',counter,file=f)
		doWeHaveAThm(counter)
		subEnd = time.time()
		print('m =',counter,'took',subEnd - subStart,'seconds',file=f)
		print('total time so far is',subEnd - totalStart,'seconds',file=f)
		print('',file=f)
	totalEnd = time.time()
	print('the total computation took',totalEnd - totalStart,'seconds',file=f)
	print('',file=f)
	f.close()

def B_ELC(last):
	"""
	Checks strong equivariant log concavity of `B_n` in degrees `1`
	through `last` for all `n`, using representation stability.
	"""
	import time
	totalStart = time.time()
	f = open('B_ELC.txt', 'a', buffering=1)
	for counter in range(1,last+1):
		subStart = time.time()
		print('Checking m =',counter,file=f)
		doWeHaveBThm(counter)
		subEnd = time.time()
		print('m =',counter,'took',subEnd - subStart,'seconds',file=f)
		print('total time so far is',subEnd - totalStart,'seconds',file=f)
		print('',file=f)
	totalEnd = time.time()
	print('the total computation took',totalEnd - totalStart,'seconds',file=f)
	print('',file=f)
	f.close()

def C_ELC(last):
	"""
	Checks strong equivariant log concavity of `C_n` in degrees `1`
	through `last` for all `n`, using representation stability.
	"""
	import time
	totalStart = time.time()
	f = open('C_ELC.txt', 'a', buffering=1)
	for counter in range(1,last+1):
		subStart = time.time()
		print('Checking m =',counter,file=f)
		doWeHaveCThm(counter)
		subEnd = time.time()
		print('m =',counter,'took',subEnd - subStart,'seconds',file=f)
		print('total time so far is',subEnd - totalStart,'seconds',file=f)
		print('',file=f)
	totalEnd = time.time()
	print('the total computation took',totalEnd - totalStart,'seconds',file=f)
	print('',file=f)
	f.close()

def D_ELC(last):
	"""
	Checks strong equivariant log concavity of `D_n` in degrees `1`
	through `last` for all `n`, using representation stability.
	"""
	import time
	totalStart = time.time()
	f = open('D_ELC.txt', 'a', buffering=1)
	for counter in range(1,last+1):
		subStart = time.time()
		print('Checking m =',counter,file=f)
		doWeHaveDThm(counter)
		subEnd = time.time()
		print('m =',counter,'took',subEnd - subStart,'seconds',file=f)
		print('total time so far is',subEnd - totalStart,'seconds',file=f)
		print('',file=f)
	totalEnd = time.time()
	print('the total computation took',totalEnd - totalStart,'seconds',file=f)
	print('',file=f)
	f.close()

def MD(last):
	"""
	Checks whether `M_n = qD_n` for all `n` up to `n = last`, and prints
	the result in MD.txt.
	"""
	import time
	totalStart = time.time()
	f = open('MD.txt', 'a', buffering=1)
	for counter in range(1,last+1):
		subStart = time.time()
		print('Checking n =',counter,file=f)
		print(MequalsD(counter),file=f)
		subEnd = time.time()
		print('n =',counter,'took',subEnd - subStart,'seconds',file=f)
		print('total time so far is',subEnd - totalStart,'seconds',file=f)
		print('',file=f)
	totalEnd = time.time()
	print('the total computation took',totalEnd - totalStart,'seconds',file=f)
	print('',file=f)
	f.close()

memoize_M = {
    1: s([1]),
}
def M(n):
	"""
	Return the symmetric function, with `q`-coefficients,
	corresponding to the graded
	`\\mathfrak{S}_n`-representation `M_n` from [MMPR21].
	"""
	if n not in memoize_M:
		result = Sym(0)
		for i in range(n-1):
			result += q^i * M_coeff_from_OT(n, i)
		memoize_M[n] = result
	return memoize_M[n]

memoize_M_coeff_from_OT = {}
def M_coeff_from_OT(n, i):
	"""
	Returns 
	"""
	if (n,i) not in memoize_M_coeff_from_OT:
	    if n == 1:
	        if i==0: 
	            memoize_M_coeff_from_OT[(n,i)] = s([1])
	        else:
	            memoize_M_coeff_from_OT[(n,i)] = Sym(0)
	    elif i > n-2:
	        memoize_M_coeff_from_OT[(n,i)] = Sym(0)
	    else:
	        result = extract_coeff(fake_OT(n), i)
	        for k in range(1, i+1):
	            left_piece = M_coeff_from_OT(n, i-k)
	            right_piece = extract_coeff(R(n,k), k)
	            result -= left_piece.inner_tensor(right_piece)
	        memoize_M_coeff_from_OT[(n,i)] = result
	return memoize_M_coeff_from_OT[(n,i)]

def extract_coeff(symm_func, i):
	"""
	Returns the `i`th coefficient of `symm_func`.
	"""
    result = Sym(0)
    for (part, old_poly) in symm_func:
        result += s(part) * old_poly[i] 
    return result

def fake_OT(n):
	"""
	a
	"""
    result = Sym(0)
    for llambda in Partitions(n):
    	llambdaDegBound = ((llambda.to_exp())[0])-2
    	if llambdaDegBound >= 0:
	        if llambda != Partition([n]):
	            L = decompose_partition_into_rectangles(llambda)
	            M = [u[1] for u in L]
	            C_rep = qC(len(llambda))    
	            for partition_list in iterate_over_partition_tuples(M):
	                reference_rep = prod([s(mu) for mu in partition_list])
	                term = reference_rep.scalar(C_rep)
	                for i in range(len(partition_list)):
	                    mu = partition_list[i]
	                    r_i = L[i][0]
	                    term *= pleth(mu,r_i)
	                result += term
    return result

def decompose_partition_into_rectangles(p):
	"""
	Returns a list consisting of the decomposition
	of the partition `p` into rectangles.
	"""
    return [(i, list(p).count(i)) for i in reversed(range(max(p)+1)) if i in p]   

def iterate_over_partition_tuples(M):
	"""
	Returns an iterator for the Cartesian product of the partitions
	of the integers in the list `M`.  In other words, this iterates
	over lists of the form 
	[partition of M[0], ..., partition of M[len(M)]].
	"""
	return xmrange_iter([Partitions(u) for u in M])

def M_compact_supp(n):
	"""
	Returns `M_n^c` from [MPY17, Proposition 3.6], which is just
	`M_n` backward.
	"""
    result = Sym(0)
    Mn = M(n)
    for i in range(n):
        coeff = extract_coeff(Mn, i)
        result += coeff * q**(2*(n-1) - i)
    return result

def truncator(gradedRep,d):
	"""
	Return the truncation of the `q`-symmetric function gradedRep
	at degree `d`, which is a polynomial of degree `d`.
	"""
	result = Sym(0)
	for i in range(0,d+1):
		result += q**i * extract_coeff(gradedRep,i)
	return result

memoize_pleth = {}
def pleth(nu,i,n=24):
	"""
	Returns the plethysm of `s_{\\nu}` with the internal (Kronecker)
	product of `M_i^c` and `R_i` (with `R_i` appropriately truncated).
	"""
	if (nu,i) not in memoize_pleth:
		if i == 1:
			memoize_pleth[(nu,i)] = s(nu)
		else:
			if nu == [1]:
				memoize_pleth[(nu,i)] = M_compact_supp(i).itensor(R(i,n-2-i))
			else:
				memoize_pleth[(nu,i)] = s(nu).plethysm(M_compact_supp(i).itensor(R(i,n-2-i)))
	return memoize_pleth[(nu,i)]

memoize_R = {}
def R(n, degree_bound):
	"""
	Returns `R_n` up to the `degree_bound` coefficient 
	(see [MPY17, Section 3.3]).
	"""
	if (n,degree_bound) not in memoize_R:
		if degree_bound >= n-2:
		    result = 0
		    for llambda in Partitions(n):
		        coeff = base_ring(1)
		        for (i, part) in enumerate(llambda):
		            coeff *= q**(i * part)
		        for (i,j) in llambda.cells():
		            coeff *= (1-q**(degree_bound + 1 + j - i))
		            hook = llambda.hook_length(i,j)
		            coeff *= sum(q**(t*hook) for t in range(int((degree_bound + 1)/hook + 2)))
		        coeff = ((1-q)*coeff).truncate(degree_bound + 1)  
		        result += coeff * s(llambda)
		    memoize_R[(n,degree_bound)] = result
		else:
			memoize_R[(n,degree_bound)] = truncator(R(n,n-2),degree_bound)
	return memoize_R[(n,degree_bound)]

"""
Bibliography:

[HR17] - 

[MPY17] -

[MMPR21] -
"""