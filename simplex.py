import numpy as np
from sympy import *
import sys

from fractions import Fraction as f
np.set_printoptions(precision=3,suppress=True)

#def printm(a):
#    """Prints the array as strings
#    :a: numpy array
#    :returns: prints the array
#    """
#    def p(x):
#        return str(x)
#    p = vectorize(p,otypes=[str])
#    print p(a)
#def roundto(A,v=2):
#    """Returns numpy array with fraction values
#    :A: numpy array
#    :returns: numpy array with fraction values
#    """
#    return np.array([map(lambda x: round(float(str(x)),v), S) for S in A])

def tofrac(A):
    """Returns numpy array with fraction values
    :A: numpy array
    :returns: numpy array with fraction values
    """
    return np.array([map(lambda x: f(str(x)), S) for S in A])

def tofloat(A, decimals=-1):
    """Returns numpy array with float values
    :A: numpy array
    :r: rounds down to r decimals
    :returns: numpy array with float values
    """
    return np.array([map(lambda x: round(float(x), decimals),S) for S in A])


def tableau(a,W=7):
    """Returns a string for verbatim printing
    :a: numpy array
    :returns: a string
    """
    if len(a.shape) != 2:
        raise ValueError('verbatim displays two dimensions')
    rv = []
    rv+=[r'|'+'+'.join('{:-^{width}}'.format('',width=W) for i in range(a.shape[1]))+"+"]
    rv+=[r'|'+'|'.join(map(lambda i: '{0:>{width}}'.format("x"+str(i+1)+" ",width=W), range(a.shape[1]-2)) )+"|"+
         '{0:>{width}}'.format("-z ",width=W)+"|"
         '{0:>{width}}'.format("b ",width=W)+"|"]
    rv+=[r'|'+'+'.join('{:-^{width}}'.format('',width=W) for i in range(a.shape[1]))+"+"]
    for i in range(a.shape[0]-1):
        rv += [r'| '+' | '.join(['{0:>{width}}'.format(str(a[i,j]),width=W-2) for j in range(a.shape[1])])+" |"]
    rv+=[r'|'+'+'.join('{:-^{width}}'.format('',width=W) for i in range(a.shape[1]))+"+"]
    i = a.shape[0]-1
    rv += [r'| '+' | '.join(['{0:>{width}}'.format(str(a[i,j]),width=W-2) for j in range(a.shape[1])])+" |"]
    rv+=[r'|'+'+'.join('{:-^{width}}'.format('',width=W) for i in range(a.shape[1]))+"+"]
    print '\n'.join(rv)

def bmatrix(a):
    """Returns a LaTeX bmatrix
    :a: numpy array
    :returns: LaTeX bmatrix as a string
    """
    if len(a.shape) > 2:
        raise ValueError('bmatrix can at most display two dimensions')
    lines = str(a).replace('[', '').replace(']', '').splitlines()
    rv = [r'\begin{bmatrix}']
    rv += [r'  ' + ' & '.join(l.split()) + '\\\\' for l in lines]
    rv +=  [r'\end{bmatrix}']
    return '\n'.join(rv)

def findpivot(tabl, basiscol):
    """
    :tabl: numpy array
    :basiscol: column index for which to do the pivot operation
    :returns the pivot row 
    """
    tabl = np.array(tabl, dtype="float64")
    pivot = float(sys.maxint)
    pivot_index = -1
    bi = (tabl[0].size)-1
    aHeight = (tabl.T[0].size)-1
    for i in range(0,aHeight):
        a = tabl[i][basiscol]
        #print a 
        if a > 0:
            # bi / ai
            t = float(tabl[i][bi])/float(a)
    #        print tabl[i][bi], " divided by ", a, " is ", t
            if t < pivot and t > 0:
                pivot = t
                pivot_index = i
    #print "returning row ",pivot_index
    return pivot_index, pivot

def simplex(tabl, latex=False, frac=True, decimals=-1, verbose=True):
    b = tabl.T[tabl[0].size-1][0:-1]
    if np.any(b < 0):
        print "Infeasible, use dual simplex"
        return
    if np.any(b == 0):
        print "Warning: Tableau has degeneracies, ie. a value in b is zero."

    pivot_row = -1
    pivot_col = -1
    pivot_value = sys.maxint
    
    #height = tabl.T[0].size
    c = tabl[-1,:-2] #cost/optimisation function

    #print c
    for col in range(0,tabl[0].size-2):
        if(c[col] > 0):
            r,v = findpivot(tabl, col)
            if v < pivot_value:
                pivot_value = v
                pivot_col = col
                pivot_row = r

    if(pivot_row == -1 or pivot_col == -1):
        print "couldnt find pivot element"
        return

    print "(largest increase?) pivot, (i: ",pivot_row," j: ",pivot_col, ") = ",pivot_value
    return enterbasis(tabl, pivot_row, pivot_col, latex=latex, frac=frac, decimals=decimals, verbose=verbose)

def dualsimplex(tabl, latex=False, frac=True, decimals=-1, verbose=True):
    tabl = np.array(tabl, dtype="float64")
    #height = tabl.T[0].size
    pivot = float(sys.maxint)
    pivot_i = -1
    pivot_j = -1

    #find infeasable rows
    b = tabl.T[tabl[0].size-1][0:-1]
    c = tabl[-1][0:-2] #cost/optimisation function
    infrows = []
    for i in range(0, b.size):
        if b[i] < 0:
            #print i,b[i]
            infrows.append(i)

    # find pivot, from infeasible rows, with negative a_ij
    for i in infrows:
        for j in range(0,c.size):
            aij = tabl[i][j]
            if aij < 0:
                t = abs(c[j]/aij)
                if t < pivot:
                    pivot = t
                    pivot_i = i
                    pivot_j = j
    if(pivot_j == -1 or pivot_i == -1):
        print "couldnt find pivot element"
        return 
	print "pivot: (i: "+str(pivot_i)+", j: "+str(pivot_j)+") = "+str(pivot)
    return enterbasis(tabl, pivot_i, pivot_j, latex=latex, frac=frac, decimals=decimals, verbose=verbose)
    


def enterbasis(tabl, row, col, latex=False, frac=True, decimals=-1, verbose=True):
    """ 
    :tabl: Matrix of the full tableau, numpy array
    :row: The row index which was chosen by pivot operation, int
    :col: The column(variable) index to enter basis, int
    :latex: Print to latex barray, optional, bool
    :frac: Convert values to fractions, optional, bool
    :returns the new tableau as a numpy array 
    """
    tabl = np.array(tabl, dtype="float64")
    if frac:
        tabl = tofrac(tabl)
    print "x",col," entering basis"
        
    basisrow = tabl[row]
    #print "basisrow: ", basisrow
    print "R",row,"= R",row,"/",basisrow[col]
    tabl[row] = basisrow = basisrow/basisrow[col]    
    #print "basisrow now: ", tabl[row]
    height = tabl.T[0].size
    
    for i in range(0,height):
        if i != row:
           # if frac:
           #     print "R",i,"= R",row,"-",f(str(float(tabl[i][col]))),"*R",row 
           # else:
            print "R",i,"= R",row,"-",float(tabl[i][col]),"*R",row 
            tabl[i] = tabl[i]-(float(tabl[i][col])*basisrow)            
            if(verbose and frac):
                tableau(tofrac(tabl))
            elif(verbose and decimals != -1):
                tableau(tofloat(tabl, decimals=decimals))
            
    if latex:
        print bmatrix(tabl)
    if frac:
        return tofrac(tabl)
    elif decimals != -1:
        return tofloat(tabl, decimals=decimals)
    return tabl

#printvalues example:
#A = np.array([ 1, 1, -1 , 0 , 1 , -1 , 0 , 1,2 , 0 , 3 , 1 , -1 , 2 , 0 , 3,-3 , 0 ,-2 , 0 , -1 , -1, 1 , -9]).reshape(3,-1)
#s.printvalues(A)
#Solution: feasable and optimal.
#x_1 = 0
#x_2 = 1
#x_3 = 0
#x_4 = 3
#x_5 = 0
#objval = 9
#y_1 = 1
#y_2 = 1
def printvalues(tabl):
	""" Call this on an optimal tableau to print the value of the variables and objective function
	:tabl: optimal tableau
	"""
	objval = tabl[-1][-1]
	b = tabl.T[-1][0:-1]
	c = tabl[-1][0:-2] #cost/optimisation function
	M = tabl[0:-1,0:-2]
	MT = M.T
	n = M[0].size-MT[0].size
	text = ""
	if np.any(b < 0):
		text += "Solution: not feasible and "
	else:
		text += "Solution: feasable and "
	if np.any(c > 0):
		text += "not optimal.\n"
	else:
		text += "optimal.\n"

	for i in range(0,n+1):
		t = list(MT[i])
		#if column/variable in basis, print corresponding b value
		if t.count(1) == 1 and t.count(0) == len(t)-1:
			j = t.index(1)
			text += "x_"+str(i+1)+" = "+str(b[j])+"\n"
		else:
			text += "x_"+str(i+1)+" = 0\n"
	text += "objval = "+str(-objval)+"\n"

	j=1
	for i in range(n, M[0].size):
		text += "y_"+str(j)+" = "+str(-c[i])+"\n"
		j += 1
	print(text)
	return

def variablesinbasis(A):
	""" Get the indices of the variables in basis
	:A: coefficient matrix
	:Returns: list of indices of variables in basis
	"""
	AT = A.T
	inbasis = []
	for i in range(0, A.shape[1]):
		t = list(AT[i])
		if t.count(1) == 1 and t.count(0) == len(t)-1:
			inbasis.append(i)
	return inbasis


def revsimplex(tabl, basis):
	""" Returns reduced costs and dual variable values only
	:tabl: tableau
	:basis: list of variables in basis
	"""
	A = tabl[0:-1,0:-2]
	c = tabl[-1][0:-2] #cost/optimisation function
	m,n = A.shape
	non_basis = []
	for i in range(n):
		if i not in basis:
			non_basis.append(i)
	A_B, A_N = (A[:,basis], A[:,non_basis])
	c_B, c_N = (c[basis], c[non_basis])
	print "TEST"
	y = np.linalg.solve(A_B,c_B)
	RES = c_N-(y.T.dot(A_N))
#	print "c_N:"
#	print c_N
#	print "y.T:"
#	print y.T
#	print "A_N:"
#	print A_N.shape
#	print A_N
#	print "RES:"
#	print RES
	for i in range(0, y.size):
		print("y_"+str(i+1)+" = "+str(y[i]))
	print "Reduced cost:"
	for i,v in zip(non_basis, RES):
		print("x_"+str(i)+" = "+str(v))
	return y, RES


def getdual(tabl):
    #height = tabl.T[0].size
    objval = tabl[-1][-1]
    b = tabl.T[-1][0:-1]
    c = tabl[-1][0:-2] #cost/optimisation function
    M = tabl[0:-1,0:-2]

    # add -z
    M = M.T
    #get coeffient  matrix, add slack vars
    M = np.hstack([M, np.identity(M.T[0].size)])
    M = np.hstack([M, np.zeros(M.T[0].size).reshape((M.T[0].size, -1))])
    # add new b from c
    M = np.hstack([M, c.reshape((-1,1))])
    # add new c from b
    z = np.zeros(M[0].size)
    for i in range(0,b.size):
        z[i]=b[i]
    z[-2] = 1
    z[-1] = objval
    return np.vstack([M,z])
    
    #M = np.hstack([M,c.reshape((-1,1))])
    #MBC = np.vstack([MB, c])
    #np.vstack(MB, 

#get the dual of a problem

##eksempel task 5 exam 2014
#equalities = ["<=", "<=", "="]
#lastConst = ["inR","inR",">=0",">=0"]
#A = np.array([[1,2,1,1,5], [3,1,-1,0,8], [0,1,1,1,1], [6,1,-1,-1,0]])

def printdual(M,E,LC, obj="max"):
	"""
	:M: Matrix
	:E: equalities pr constraint
	:LC: Variable constraints
	"""
	M = M.T
	for i in range(0,len(M)):
		text = ""
		counter = 1
		for j in range (0,len(M[i])):
			if j == len(M[i])-2 and i == len(M)-1: #when reaching the last row and column do this
					text += (str(M[i][j])+"y_"+str(counter+j))
					if obj=="max":
						text = "objfunc: min " + text
					else:
						text = "objfunc: max " + text
					for k in range(0,len(M[i])-1):
						print lastCon(E,k, obj)
					break # just so i wont run til end
			if j == len(M[i])-1:
				text += addEquality(M,LC,i,j, obj)
			else :
				if j == len(M[i])-2:
					text += (str(M[i][j])+"y_"+str(counter+j))
				else :
					text += (str(M[i][j])+"y_"+str(counter+j)+" + ")
		print text
	# return A.T, new b, new objfunc
	#return M[:-1,:-1],M.T[-1,:],M[-1,:-1]  
	return M

def addEquality(M,LC,i,j, obj):
	if obj == "max":
		if LC[i] == "<=0" or LC[i] == "=<0":
			return " <= " + str(M[i][j])
		if LC[i] == "=>0" or LC[i] == ">=0":
			return " >= " + str(M[i][j])
		if LC[i] == "inR":
			return " = " + str(M[i][j])
	elif obj == "min":
		if LC[i] == "<=0" or LC[i] == "=<0":
			return " >= " + str(M[i][j])
		if LC[i] == "=>0" or LC[i] == ">=0":
			return " <= " + str(M[i][j])
		if LC[i] == "inR":
			return " = " + str(M[i][j])


def lastCon(E,i, obj):
	#print E, " and ", i
	if obj == "max":
		if E[i] == "<=" or E[i] == "=<":
			return "y" + str(i+1) + " >= 0"
		elif E[i] == ">=" or E[i] == "=>":
			return "y" + str(i+1) + " <= 0"
		elif E[i] == "=":
			return "y" + str(i+1) + " in R"
	elif obj == "min":
		if E[i] == "<=" or E[i] == "=<":
			return "y" + str(i+1) + " <= 0"
		elif E[i] == ">=" or E[i] == "=>":
			return "y" + str(i+1) + " >= 0"
		elif E[i] == "=":
			return "y" + str(i+1) + " in R"

def pad_to_square(a, pad_value=0):
  m = a.reshape((a.shape[0], -1))
  padded = pad_value * np.ones(2 * [max(m.shape)], dtype=m.dtype)
  padded[0:m.shape[0], 0:m.shape[1]] = m
  return padded

def printLHS(res,b):
	counter = 0
	print res
	rank = np.linalg.matrix_rank(res)
	res = res[:rank,rank:]
	res = res*-1
	res = np.vstack([res,np.identity(res.shape[1])])
	text = ""
	for i in range(0,len(res.T)):
		if (i==0):
			text += "a"+str(i+1)+" * " + str(res.T[i])
		else:
			text += " + a"+str(i+1)+" * " + str(res.T[i])
	print b,"+", text

def printRHS(res,b):
	counter = 0
	for i in range(0,len(res)):
		text = ""
		curA = [elem for elem in res[i] if elem != 0 and elem !=1]
		for j in range(0,len(curA)):
			if(j == 0):
				text += "x"+str(i+1)+" = "+str(b[i])+"-"+str(curA[j])+"a_"+str(j+1)+"-"
			if(j == len(curA)-1): #reaching last
				text += str(curA[j])+"a_"+str(j+1)
			else:
				text += str(curA[j])+"a_"+str(j+1)+"-"
		if curA == []:
			text += "x"+str(i+1)+ " = a_"+str(counter+1)
			counter += 1
		print text

def nonSingularScalarForm(M):
	res, bres = np.array(Matrix(tofrac(M)).rref())
	b = res.T[-1]
	res = np.delete(res,-1,axis=1)
	iden = np.identity(res.shape[0])
	
	if np.array_equal(res,iden) or np.array_equiv(res,iden):
		print "square matrix, vector x ="
		print bres 
		return	
	res = pad_to_square(res)

	while(len(b)!= len(res[0])):
		b = np.append(b,0)

	print "x: "
	print "RHS:"
	printRHS(res,b)
	print "LHS:"
	printLHS(res,b)
	return res

#nonSingularScalarForm(A)

def getRank(M):
	# without b vector
	return np.linalg.matrix_rank(np.delete(M,-1,axis=1))

def lin(M):
  #set given by linear span of the column of the coefficient matrix A, from rref, since rank is 2, 
  #hence the space of solutions is spanned by two vectors, column of the leading ondes ing the rref
  #lin(A)=lin({v_1,..v_n})
	res = np.array(Matrix(tofrac(M)).rref()[0])
	#rank = np.linalg.matrix_rank(res)
	rank = getRank(M)
	K = variablesinbasis(res)
	R = np.array([])
	for x in K:
		R = np.append(R,M.T[x])
	print "Linearly independent columns:"
	print R.reshape((rank,-1)).T

def geteigen(A):
	""" Returns a tuple of eigenvalues and eigenvectors
	:A: matrix
	:returns evals, evecs
	"""
	# remember you probably should transpose the eigenvectors array
	evals, evecs = np.linalg.eig(A)
	return evals, np.array([ev/max(ev) for ev in evecs.T]).T

def addslackandz(M):
	""" Adds slack variables and the -z vector
    """
	return np.concatenate([M[:,:-1],np.identity(M.shape[0]), M[:,-1].reshape(M.shape[0],-1)], axis=1)

def autosimplex(tabl_orig, verbose=False, frac=True, decimals=-1):
	""" Returns optimal tnd feasible tableau"""
	tabl = np.array(tabl_orig)
	#objval = tabl[-1][-1]
	b = tabl.T[-1][0:-1]
	c = tabl[-1][0:-2] #cost/optimisation function
    #M = tabl[0:-1,0:-2]
	while(1): 
		if np.any(0>b): #infeasible -> feasible : PHASE 1
			if verbose:
				print "PHASE 1: Infeasible, using dual simplex"
			tabl = dualsimplex(tabl, verbose=verbose, frac=frac, decimals=decimals )
		elif np.any(c>0):# -> optimal, PHASE 2
			if verbose:
				print "PHASE 2: not optimal, using simplex"
			tabl = simplex(tabl, verbose=verbose, frac=frac, decimals=decimals)
		else: break
		if tabl == None:
			print "probably unbounded, or ??? :( "
			break
		b = tabl.T[-1][0:-1]
		c = tabl[-1][0:-2] #cost/optimisation function
	#printvalues(tabl)
	return tabl






##example of complementary slackness check:
#DK, yvec, deq, dlc, obj = s.complslacknesscheck(K,[3,-1,0,2], ["<=","<=","="], ["inR","inR",">=0",">=0"], obj="max")
## We then check if the values hold in the dual also, by calling again on the result:
#s.complsclacknesscheck(DK,yvec,deq,dlc,obj)
## in this case the dual is infeasible,indicating that the primal is unbounded.

def complslacknesscheck(A,xvec,eq,lc, obj="max"):
	""" Check if given x-values makes solution feasible
	:A: problem without slackvariables
	:xvec: x values as vector
	:eq: equalities for constraints
	:lc: variable constraints
	"""
	
	objval = A[-1][-1]
	b = A.T[-1][0:-1] #primal b
	c = A[-1][0:-1] #cost/optimisation function
	M = A[0:-1,:-1] # primal coefficient matrix
	DM = M.T
	
	M1 = M.dot(xvec)
	print "constraints:"
	for x,e,rhs in zip(M1,eq,b):
		text = str(x)+e+str(rhs)+" is "
		if e == ">=" or e == "=>":
			if x >= rhs: 
				text += "OK"
				print text
				continue
		elif e == "<=" or e == "=<":
			if x <= rhs:
				text += "OK"
				print text
				continue
		elif e == "=" or e == "==":
			if x == rhs:
				text += "OK"
				print text
				continue
		text += "NOT OK!   <---- infeasible"
		print text
	
	print ""
	for i, x,con in zip(range(1,len(xvec)+1),xvec, lc): 
		if con == ">=0" or con == "=>":
			if x>=0:
				print "x_"+str(i)+": ",x,con, " OK"
			else:
				print "x_"+str(i)+": ",  x,con, "NOT  OK <-----"
		elif con == "<=0" or con == "=<0":
			if x<=0:
				print "x_"+str(i)+": ",x,con, " OK"
			else:
				print "x_"+str(i)+": ",  x,con, "NOT  OK <-----"
		elif con == "inR":
				print "x_"+str(i)+": ",x,con, " OK"
	print ""
		

	#for i,m,v in zip(np.array(range(1, M1.size+1)), M1,b):
	#	print "x_"+str(i)+" = "+str(m-v)

	print "The complementary slackness conditions are:"

	#y part
	for i,row in zip(range(1, DM.shape[0]+1), DM):
		text =  "y_"+str(i)+" * ( "
		res = ""
		res_val = 0
		for j,x,v in zip(range(1,len(xvec)+1),xvec,row):
			if v < 0:
				text += str(v)+"x_"+str(j)+" "
			else:
				text += "+"+str(v)+"x_"+str(j)+" "
			if v!=0 and x!=0:
				res_val += v*x
				
		
		text += ") = 0"
		text +=" = "+str(res_val)
		if res_val == 0:
			text += "  ==> "+"y_"+str(i)+" = 0"
		print text
	print ""
		
	
	#x part
	basis = []
	for i, x,row,v in zip(range(1,len(xvec)+1),xvec, DM, c):
		print "x_"+str(i)+":"
		text = str(x)+" * ( " 
		for j,q in zip(range(1,row.size+1),row):
			if q < 0:
				text += str(q)+"y_"+str(j)+" "
			else:
				text += "+"+str(q)+"y_"+str(j)+" "
			# result
		if (-v < 0):
			text += str(-v)+" ) = 0"
		else:
			text += "+"+str(-v)+" ) = 0"

		#print text

		text += " <=> "
		if x != 0:
			basis.append(i-1)
			for j,q in zip(range(1,row.size+1),row):
				if q < 0:
					text += str(q)+"y_"+str(j)+" "
				elif q == 0:
					text = text
				else:
					text += "+"+str(q)+"y_"+str(j)+" "
			
			text += "= "+str(v)
		else:
			text += "0"
		
		print text
		print ""
	yvec = np.array(Matrix(A.T[basis]).rref()[0][:,-1]).reshape(1,-1)[0]
	#yvec = np.array(yvec).reshape(1,-1)

	for i,y in zip(range(1,len(yvec)+1),yvec):
		print "y_"+str(i)+" = "+str(y)

	deq, dlc, obj = getdualconstraints(eq,lc, obj)
	return A.T, yvec, deq, dlc, obj 


def getdualconstraints(eq,lc, obj="max"):
	dlc = []
	deq	= []
	for e in eq:
		if obj=="max":
			if e == "<=" or e == "=<":
				dlc.append(">=0")
			elif e == ">=" or e == "=>":
				dlc.append("<=0")
			elif e== "=":
				dlc.append("inR")
		elif obj=="min":
			if e == "<=" or e == "=<":
				dlc.append("<=0")
			elif e == ">=" or e == "=>":
				dlc.append(">=0")
			elif e== "=":
				dlc.append("inR")
	
	for c in lc:
		if obj=="max":
			if c == "<=0" or c == "=<0":
				deq.append("<=")
			elif c == ">=0" or c == "=>0":
				deq.append(">=")
			elif c == "inR":
				deq.append("=")
		elif obj=="min":		
			if c == "<=0" or c == "=<0":
				deq.append(">=")
        	elif c == ">=0" or c == "=>0":
        		deq.append("<=")
        	elif c == "inR":
        		deq.append("=")

	if obj=="max":
		return deq, dlc, "min"
	else:
		return deq, dlc, "max"



#for row operations, do f.ex.:
# A[i,:]=A[i,:]+A[k,:]
# or explicitly
# A[0,:]=A[0,:]+A[3,:]
## Row operations and dual simplex example:
#T4 = np.array([[ 1.   ,  0.   ,  0.05 ,  0.25 ,  0.   ,  0.   ,  0.   ,  0.   ,
	#         1.5  ],
	#       [ 0.   ,  0.   , -0.075,  0.125,  0.   ,  1.   ,  0.   ,  0.   ,
	#         0.25 ],
#       [ 0.   ,  1.   ,  0.075, -0.125,  0.   ,  0.   ,  0.   ,  0.   ,
#         1.75 ],
#       [ 0.   ,  0.   ,  0.275, -0.125,  1.   ,  0.   ,  0.   ,  0.   ,
#         5.75 ],
#       [ 0.   , -1.   ,  0.   ,  0.   ,  0.   ,  0.   ,  1.   ,  0.   , -2.   ],
#       [ 0.   ,  0.   , -0.175, -0.375,  0.   ,  0.   ,  0.   ,  1.   ,
#        -4.75 ]])
#T4 = toFrac(T4)
#tableau(T4)
#T4b = T4
#T4b[-2,:]=T4b[2,:]+T4b[-2,:]
#tableau(T4b)
#use simplex to find optimal solution
#tableau(autosimplex(T4b))
