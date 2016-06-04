import numpy as np
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
    
    height = tabl.T[0].size
    c = tabl[height-1][0:-2] #cost/optimisation function

    for col in range(0,tabl[0].size-2):
        if(c[col] > 0):
            r,v = findpivot(tabl, col)
            if v < pivot_value:
                v = pivot_value
                pivot_col = col
                pivot_row = r

    if(pivot_row == -1 or pivot_col == -1):
        print "couldnt find pivot element"
        return

    print "(largest increase?) pivot, i: ",pivot_row," j: ",pivot_col
    return enterbasis(tabl, pivot_row, pivot_col, latex=latex, frac=frac, decimals=decimals, verbose=verbose)

def dualsimplex(tabl, latex=False, frac=True, decimals=-1, verbose=True):
    tabl = np.array(tabl, dtype="float64")
    height = tabl.T[0].size
    pivot = float(sys.maxint)
    pivot_i = -1
    pivot_j = -1

    #find infeasable rows
    b = tabl.T[tabl[0].size-1][0:-1]
    c = tabl[height-1][0:-2] #cost/optimisation function
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
                print "*** HERE ***"
                tableau(tofrac(tabl))
            elif(verbose and decimals != -1):
                print "*** HERE2 ***"
                tableau(tofloat(tabl, decimals=decimals))
            
    if latex:
        print bmatrix(tabl)
    if frac:
        return tofrac(tabl)
    elif decimals != -1:
        return tofloat(tabl, decimals=decimals)
    return tabl


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
#use the dualsimplex to find optimal solution
#tableau(dualsimplex(T4b))
