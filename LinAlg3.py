
import matplotlib.pyplot as plt     # optional, just for plotting some error data
from random import uniform
from datetime import datetime
from math import cos, pi, sqrt, log


class Matrix:
    'fields: body, rows, cols'
    def __init__(self, x):
        r = len(x)
        if isinstance( x[0], (float,int) ):
            self.body = [[j] for j in x]
            self.rows = r
            self.cols = 1
        else:
            c = []
            for i in range(r):
                c.append(len(x[i]))
            if not len(set(c))==1:
                raise ValueError('The entries of the main argument must all be of equal length.')
            self.body = x
            self.rows = r
            self.cols = c[0]
#
#   end Matrix
#



def ZeroMat(m,n):
    x = [[0] * n for i in range(m)]
    return Matrix(x)
#
#   end Zero_Mat
#


def EyeBall(n):
    mat = ZeroMat(n,n)
    for i in range(1,n+1):
        Assign(mat, i, i, 1)
    return(mat)
#
#   end EyeBall
#


#
##
###     print a matrix in a readable format
##
#

def Print(mat, roun=3):
    if isinstance(mat, Matrix):
        m = mat.rows
        n = mat.cols
        lengths = [0 for i in range(m*n)]
        for i in range(m):
            for j in range(n):
                lengths[i*n+j] = len(str( round(mat.body[i][j], roun) ))
        space = 1+max(lengths)
        z = [[0] * n for i in range(m)]
        for i in range(m):
            for j in range(n):
                z[i][j] = round(mat.body[i][j], roun)
        print('')
        for i in range(m):
            x = str(z[i])
            x = x.lstrip('[')
            x = x.rstrip(']')
            x = x.split(', ')
            y = '   '.join([' ', x[0].rjust(space)])
            for j in range(1, len(x)):
                y = ' '.join([y, x[j].rjust(space)])
            print(y)
        print('')
    else:
        raise ValueError('Argument passed is not of class Matrix.')
        
#
#   end Print
#



#
##
###     assign a value the the (i,j)-th entry of a matrix
##
#

def Assign(mat,i,j,val):
    if not isinstance(mat, Matrix):
        raise ValueError('Argument passed is not of class Matrix.')
    elif not isinstance( val, (float,int) ):
        raise ValueError('Please make val numeric (float or int).')
    elif ( i-1<mat.rows and j-1<mat.cols and 0<i and 0<j and isinstance(i,int) and isinstance(j,int) ):
        mat.body[i-1][j-1] = val
    else:
        raise ValueError('(i,j) out of bounds. Pleas make i and j positive integers, at most', mat.rows, 'and', mat.cols, 'respectively.')
#
#   end Assign
#





#
##
###     multiply a matrix Y by X on the left, i.e., XY
##
#

def Multiply(X,Y):
    if not isinstance(X, Matrix):
        raise ValueError('Please make X a matrix in "Multiply(X,Y)."')
    if not isinstance(Y, Matrix):
        raise ValueError('Please make Y a matrix in "Multiply(X,Y)."')
    if not X.cols == Y.rows:
        raise ValueError('Dimensions incompatible. Multiply(X,Y) requries X.cols == Y.rows.')
    common = X.cols
    dummy = ZeroMat(X.rows, Y.cols)
    for i in range(X.rows):
        for j in range(Y.cols):
            val = [0]*common
            for k in range(common):
                val[k] = X.body[i][k] * Y.body[k][j]
            Assign(dummy, i+1, j+1, sum(val))
    return(dummy)



def t(mat):
    if not isinstance(mat, Matrix):
        ValueError('The argument passed is not of class "Matrix."')
    ret = ZeroMat(mat.cols, mat.rows)
    for i in range(mat.cols):
        for j in range(mat.rows):
            Assign(ret, i+1, j+1, mat.body[j][i])
    return(ret)
#
#   end t
#




#
##
###     check if a matrix is upper triangular
##
#

def isUT(M):
    
    cc = []
    if not isinstance(M, Matrix):
        ValueError('The argument passed is not of class "Matrix."')
    m = M.rows
    n = M.cols
    if m==1 and n==1:
        ValueError('Question is moot. The matrix is 1-by-1.')
    
    for i in range(1,min(m,n)):
        for j in range(i):
            cc.append(M.body[i][j])
    
    if m>n:
        for i in range(n,m):
            for j in range(n):
                cc.append(M.body[i][j])
    
    if len(set(cc))==1 and cc[0]==0:
        return(True)
    else:
        return(False)

#
#   end isUT
#




def FirstNonZero(vector):
    which = -1
    for q in range(len(vector)):
        if vector[q]==0:
            which = which+1
        else:
            break
    return which+1


def RevRows(M):
    if not isinstance(M, Matrix):
        raise ValueError('Please make M in "RevRows(M)" an object of class Matrix.')
    reved = []
    m = M.rows
    for q in range(m):
        reved.append( M.body[m-q-1] )
    reved = Matrix(reved)
    return reved


def NonRedundantCopy(M):
    A = ZeroMat(M.rows, M.cols)
    for i in range(M.rows):
        for j in range(M.cols):
            Assign(A, i+1, j+1, M.body[i][j])
    return A

#
##
###
####    convert a triangular matrix to RREF
###
##
#

def UpperTriRREF(M):
    if not isUT(M):
        raise ValueError('Please make M an upper trianguar in "UpperTriRREF(M)".')
    m = M.cols
    n = M.rows
    for i in range(n):
        i=n-i-1
        q = FirstNonZero( M.body[i] )
        if q<m:
            aiq = M.body[i][q]
            for j in range(q+1,m):
                if not M.body[i][j]==0:
                    M.body[i][j] = M.body[i][j]/aiq
            M.body[i][q] = 1
            for ii in range(i):
                scale = M.body[ii][q]
                M.body[ii][q] = 0
                for j in range(q+1,m):
                    if not M.body[i][j]==0:
                        M.body[ii][j] = M.body[ii][j] - scale*M.body[i][j]





def TriInvert(M):
    A = NonRedundantCopy(M)
    reved = False
    if not isUT(A):
        if isUT(t(A)):
            reved = True
            A = RevRows(t(RevRows(t(A))))
        else:
            raise ValueError('Pleas make M a triangular matrix in "TriIvert(M)".')
    m = M.rows
    n = M.cols
    if not m==n:
        raise ValueError('Cannot invert a non-square matrix.')
    if A.body[0][0]==0:
        raise ValueError('Matrix is not invertible.')
    for i in range(1,n):
        q = FirstNonZero(A.body[i])
        if not q==i:
            raise ValueError('Matrix is not invertible.')
    dummy = ZeroMat(n, 2*n)
    for i in range(n):
        dummy.body[i][n+i] = 1
        for j in range(n):
            dummy.body[i][j] = A.body[i][j]
    UpperTriRREF(dummy)
    A = dummy.body
    for i in range(n):
        A[i] = A[i][n:]
    A = Matrix(A)
    if reved:
        A = RevRows(t(RevRows(t(A))))
    return A

#
##
###
####
#####   matrix factorization
####
###
##
#

def Factorize(A, method='LU', OnesOnL=True, DoPrint=True, Return=False, both=False):
    supported = ['LU', 'Cholesky']
    if both:
        method = 'Cholesky'
        DoPrint = False
    if not method in supported:
        print('Currently supported methods of matrix factorization are as follows')
        print(*supported)
        raise ValueError('The method ' + method + ' is not supported.')
    if not isinstance(A, Matrix):
        A = Matrix(A)
        print('Warning: M was attempted to be coerced to a matrix. The user should verify that the output is sensible.')
    if not A.rows==A.cols:
        raise ValueError('The matrix M must be square.')
    n = A.rows
    if (method=='Cholesky'):
        if not A.body==t(A).body:
            raise ValueError('Cannot perform Cholesky factorization because A is not symmetric.')
        if not OnesOnL:
            OnesOnL = True
            print('Warning: the argument OnesOnL was set to true for Cholesky factorization.')
    if True:    
        #
        #
        #   LU factorization
        #
        #
        L = ZeroMat(n,n)
        U = ZeroMat(n,n)
        a11 = A.body[0][0]
        if a11 == 0:
            raise ValueError('LU factorization not possible.')
        if OnesOnL:
            Assign(L,1,1,1)
            Assign(U,1,1,a11)
            for j in range(2,n+1):
                thing = A.body[j-1][0]/a11
                Assign(U, 1, j, A.body[0][j-1])
                Assign(L, j, 1, thing)
        else:
            Assign(L,1,1,a11)
            Assign(U,1,1,1)
            for j in range(2,n+1):
                thing = A.body[0][j-1]/a11
                Assign(U, 1, j, thing)
                Assign(L, j, 1, A.body[j-1][0])
        for i in range(2, n):
            diagonalentry = A.body[i-1][i-1]
            for j in range(i-1):
                diagonalentry = diagonalentry - L.body[i-1][j]*U.body[j][i-1]
            if diagonalentry==0:
                raise ValueError('LU factorization not possible.')
            if OnesOnL:
                Assign(L, i, i, 1)
                Assign(U, i, i, diagonalentry)
                for j in range(i+1, n+1):
                    uij = A.body[i-1][j-1]
                    lji = A.body[j-1][i-1]
                    for h in range(i-1):
                        uij = uij - L.body[i-1][h]*U.body[h][j-1]
                        lji = lji - L.body[j-1][h]*U.body[h][i-1]
                    Assign(U, i, j, uij)
                    thing = lji/diagonalentry
                    Assign(L, j, i, thing)    
            else:
                Assign(L, i, i, diagonalentry)
                Assign(U, i, i, 1)
                for j in range(i+1, n+1):
                    uij = A.body[i-1][j-1]
                    lji = A.body[j-1][i-1]
                    for h in range(i-1):
                        uij = uij - L.body[i-1][h]*U.body[h][j-1]
                        lji = lji - L.body[j-1][h]*U.body[h][i-1]
                    thing = uij/diagonalentry
                    Assign(U, i, j, thing)
                    Assign(L, j, i, lji)
        ann = A.body[n-1][n-1]
        for h in range(n-1):
            ann = ann - L.body[n-1][h]*U.body[h][n-1]
        if ann==0:
            raise ValueError('Factorization not possible.')
        if OnesOnL:
            Assign(L, n, n, 1)
            Assign(U, n, n, ann)
        else:
            Assign(L, n, n, ann)
            Assign(U, n, n, 1)
    if method=='Cholesky':
        #
        #
        #   obtain Cholesky factorization from the LU factorization
        #
        #
        D = ZeroMat(n,n)
        for i in range(n):
            D.body[i][i] = pow(U.body[i][i], 0.5)
        C = Multiply(L,D)
        for i in range(n):      # optional
            for j in range(n):
                if C.body[i][j] == 0.0:
                    C.body[i][j] = 0
    #
    #
    #   what to do with the results
    #
    #
    if DoPrint:
        if method=='LU':
            print('L:')
            Print(L)
            print('U:')
            Print(U)
        if method=='Cholesky':
            print('The given matrix factors into CC^T where C is:')
            Print(C)
    if Return:
        if method=='LU':
            return(L, U)
        if method=='Cholesky':
            return(C)
    if both:
        return(L, U, C)



def HomeworkMatrix(m):      # definitely not my best work
    A = ZeroMat(m**2, m**2)
    Assign(A, 1, 1, 4)
    Assign(A, 1, 2, -1)
    Assign(A, m, m, 4)
    Assign(A, m, m-1, -1)
    for i in range(2,m):
        Assign(A, i, i-1, -1)
        Assign(A, i, i, 4)
        Assign(A, i, i+1, -1)
    
    Assign(A, m*(m-1)+1, m*(m-1)+1, 4)
    Assign(A, m*(m-1)+1, m*(m-1)+2, -1)
    Assign(A, m*m, m*m, 4)
    Assign(A, m*m, m*m-1, -1)
    for i in range(2,m):
        Assign(A, m*(m-1)+i, m*(m-1)+i-1, -1)
        Assign(A, m*(m-1)+i, m*(m-1)+i, 4)
        Assign(A, m*(m-1)+i, m*(m-1)+i+1, -1)
    
    for i in range(m):
        A.body[i][m+i] = -1
        A.body[m*(m-1)+i][m*(m-2)+i] = -1
    
    thing = [-1,4,-1]
    for i in range(1,m):
        for j in range(1,m):
            j = j
            for q in range(2):
                Assign(A, i*m+j, (i-1+q)*m+j, thing[q])
        for h in range(2,m):
            Assign(A, i*m+h, i*m+h-1, -1)
            Assign(A, i*m+h, i*m+h+1, -1)
        Assign(A, m*(i+1), m*i, -1)
        Assign(A, m*(i+1), m*(i+1), 4)
        Assign(A, m*(i+1), m*(i+1)-1, -1)
        Assign(A, m*i+1, m*i+2, -1)
    for h in range(2,m):
        for p in range(1,m+1):
            Assign(A, (h-1)*m+p, h*m+p, -1)
#        Assign(A, i*m, i*m+1, -1)
#       Assign(A, i*(m+1)-1, i*m-1, -1)
    return A



for m in [2,3]:
    ('   m=' + str(m) + '   ').center(30, '-')
    Factorize( HomeworkMatrix(m) )


print('Press enter to continue.')
f = input()


for m in [2,3]:
    ('   m=' + str(m) + '   ').center(30, '-')
    Factorize( HomeworkMatrix(m), 'Cholesky' )



print('Press enter to continue.')
f = input()

for m in [2,3]:
    ('   m=' + str(m) + '   ').center(30, '-')
    LU = Factorize( HomeworkMatrix(m), DoPrint=False, Return=True )
    x = Multiply( TriInvert(LU[1]) , TriInvert(LU[0]) )
    print('The inverse matrix in question is:')
    Print(x)


print('Press enter to continue.')
f = input()



def InfNorm(A,B):
    if not (isinstance(A, Matrix) and isinstance(B, Matrix)):
        raise ValueError('A and B should be objects of class Matrix.')
    if not (A.rows==B.rows and A.cols==B.cols):
        raise ValueError('The dimensions of A and B must agree.')
    measurements = []
    for i in range(A.rows):
        take = 0
        for j in range(A.cols):
            take = take+abs(A.body[i][j]-B.body[i][j])
        measurements.append( take )
    return max(measurements)
    



def AccuracyTest(m):
    A = HomeworkMatrix(m)
    x = Factorize(A, both=True)
    measurement = []
    measurement.append( InfNorm(A, Multiply(x[0],x[1])) )
    measurement.append( InfNorm(A, Multiply(x[2],t(x[2]))) )
    measurement.append( InfNorm(EyeBall(m**2), Multiply(A, Multiply(TriInvert(x[1]), TriInvert(x[0])))) )
    return measurement

data = []
teston = [w for w in range(2,11)]
#for w in [12, 14, 16, 18, 20, 25]:
#    teston.append(w)



for m in teston:
    start = datetime.now()
    data.append( AccuracyTest(m) )
    print('m='+str(m) + ' completed. Compuation time: ' + str(datetime.now()-start) + '.')



print(
    ''.center(10,' ')
    + '||A-LU||'.rjust(17)
    + '||A-CC^T||'.rjust(17)
    + '||AA^{-1}-I||'.rjust(17)
    )
for q in range(len(teston)):
    print(
        ('m='+str(teston[q])).center(10,' ')
        + str(round(data[q][0],20)).rjust(17)
        + str(round(data[q][1],20)).rjust(17)
        + str(round(data[q][2],20)).rjust(17)
    )



print('Press enter to continue.')
f = input()

def ExtractCol(mat, col):
    return [mat[i][col-1] for i in range(len(mat))]
    

def loglist(vector):
    return [log(q) for q in vector]




fig, ax = plt.subplots()
fig.set_tight_layout(True)
ax.plot(teston[1:], ExtractCol(data[1:],1), c='C1')
ax.plot(teston[1:], ExtractCol(data[1:],2), c='C3')
ax.plot(teston[1:], ExtractCol(data[1:],3), c='C4')
ax.set(xlabel='m', ylabel='Errors', title='The Error Data')
ax.grid()
plt.show()


fig, ax = plt.subplots()
fig.set_tight_layout(True)
ax.plot(loglist(teston[1:]), loglist(ExtractCol(data[1:],1)), c='C1')
ax.plot(loglist(teston[1:]), loglist(ExtractCol(data[1:],2)), c='C3')
ax.plot(loglist(teston[1:]), loglist(ExtractCol(data[1:],3)), c='C4')
ax.set(xlabel='ln(m)', ylabel='ln(Errors)', title='The Same Data: Log-Log')
ax.grid()
plt.show()





#
##  solve Ax=b by backward substitution when A is upper triangular
#

def BackSubUT(A,b):
    if not isinstance(A, Matrix):
        A = Matrix(A)
        print('Warning: M was attempted to be coerced to a matrix. The user should verify that the output is sensible.')
    if not isUT(A):
        raise ValueError('Please make A upper triangular.')
    C = NonRedundantCopy(A).body
    n = A.cols
    if not len(b)==len(C):
        raise ValueError('b must have length equal to the number of rows of A.')
    for q in range(len(b)):
        C[q].append( b[q] )
    C = Matrix(C)
    UpperTriRREF(C)
    result = []
    for q in range(len(b)):
        result.append( C.body[q][n] )
    return result



#
##  solve Ax=b by forward substitution when A is lower triangular
#

def ForSubLT(A,b):
    if not isinstance(A, Matrix):
        A = Matrix(A)
        print('Warning: M was attempted to be coerced to a matrix. The user should verify that the output is sensible.')
    C = NonRedundantCopy(A)
    C = RevRows(t(RevRows(t(C))))
    if not isUT(C):
        raise ValueError('Please make A lower triangular.')
    n = A.cols
    C = RevRows(C).body
    if not len(b)==len(C):
        raise ValueError('b must have length equal to the number of rows of A.')
    for q in range(len(b)):
        C[q].append( b[q] )
    C = RevRows(Matrix(C))
    UpperTriRREF(C)
    C = RevRows(C)
    result = []
    for q in range(len(b)):
        result.append( C.body[q][n] )
    return result


#
##  solve Ax=b using LU factorization
#

def MySolve(A,b):
    LU = Factorize(A, Return=True, DoPrint=False)
    y = ForSubLT(LU[0],b)
    x = BackSubUT(LU[1],y)
    return x
    

A = Matrix([ [-0.5,9,-2,1], [-1.5,30,-12,0], [1,-15,0,-4], [0,-6,18,8] ])
b = [3,3,2,-4]
x = Matrix(MySolve(A,b))
Print(x)
Print(Multiply(A,x))
Print(Matrix(b))


print('Press enter to continue.')
f = input()




def MaxAbs(vector):
    return max([ abs(q) for q in vector ])



data = []
teston = [w for w in range(2,16)]
#teston = [25,30]
for w in teston:
    if w==2:
        print('')
        print(
            ''.center(10,' ')
            + 'Residual'.rjust(17)
            + 'Corrected'.rjust(17)
        )
    IV = [1 for q in range(w**2)]
    H = HomeworkMatrix(w)
    soln = MySolve(H,IV)
    thing = Multiply( H, t(Matrix([soln])) )
    residual = []
    for q in range(w**2):
        residual.append( IV[q]-thing.body[q][0] )
    corrector = MySolve(H,residual)
    corrected = []
    for q in range(w**2):
        corrected.append( soln[q]+corrector[q] )
    thing = Multiply( H, t(Matrix([corrected])) )
    reresidual = []
    for q in range(w**2):
        reresidual.append( IV[q]-thing.body[q][0] )
    print(
        ('m='+str(w)).center(10,' ')
        + str(round(MaxAbs(residual),20)).rjust(17)
        + str(round(MaxAbs(reresidual),20)).rjust(17)
    )
    if w==max(teston):
        print('')


print('Press enter to continue.')
f = input()







#
##
###
####    iterative methods
###
##
#



def ModSqr(vector):
    return sum( x**2 for x in vector )


def meg(m):
    return 2./( 1 + sqrt(1-cos(pi/(m+1))**2) )


def Mes(A,x,b):      # ||b-Ax||
    return InfNorm( Matrix(b), Multiply(A, Matrix(x)) )



def Jacobi(A, b, naught=None, ShowProgress=True, lo=0, hi=0, times=None):
    #
    #   preamble
    #
    if not isinstance(A, Matrix):
        A = Matrix(A)
        print('Warning: A was attempted to be coerced to class Matrix.')
    if naught is None:
        naught = [uniform(lo,hi) for w in range(len(b))]
    x = naught
    Measure = Mes(A,x,b)
    #
    #   define the iterative step
    #
    def identical():
        oldx = [w for w in x]
        for i in range(len(b)):
            thing = [w for w in range(len(b))]
            thing.pop(i)
            thing = sum(
                    [A.body[i][j]*oldx[j] for j in thing]
                )
            x[i] = (b[i]-thing)/A.body[i][i]
    #
    #   iterate
    #
    k = 0
    if times is None:
        while (k<20000 and Measure>MaxAbs(x)*1e-8):
            identical()
            k += 1
            Measure = Mes(A,x,b)
            if ShowProgress:
                print('k=' + str(k) + '; ' + str(Measure))
    else:
        while k<times:
            identical()
            k += 1
            Measure = Mes(A,x,b)
            if ShowProgress:
                print('k=' + str(k) + '; ' + str(Measure))
    #
    #   what to do with the results
    #
    return {'soln' : x, 'times' : k}






##    successive over-relaxation, including Gauss-Siedel as the case omega==1

def SOR(A, b, N=100, naught=None, Forget=True, Plot=False, While=True, lo=0, hi=0, ShowProgress=True, omega=1):
    #
    #
    #   preamble
    #
    #
    if While and (not Forget):
        raise ValueError("The settings 'Forget=False' and 'While=True' are currently incompatible; sorry.")
    if Forget and Plot:
        Plot = False
        print('Waring: Plot set to False automatically because Forget is set to True.')
    if not isinstance(A, Matrix):
        A = Matrix(A)
        print('Warning: A was attempted to be coerced to class Matrix.')
    if naught is None:
        naught = [uniform(lo,hi) for w in range(len(b))]
    if Forget:
        x = naught
    else:
        x = [naught]
    if ShowProgress and (not Forget):
        ShowProgress = False
        print('Waring: ShowProgress set to False automatically because Forget is set to False.')
    if not (0<omega and omega<2):
        raise ValueError("Are you sure about that value of omega? Please, see Kahan's theorem.")
    #
    #
    #   define the iterative step
    #
    #
    def identical():
        if Forget:
            for i in range(len(b)):
                thing = [w for w in range(len(b))]
                thing.pop(i)
                thing = sum(
                        [A.body[i][j]*x[j] for j in thing]
                    )
                if omega==1:
                    x[i] = (b[i]-thing)/A.body[i][i]
                else:
                    x[i] = (1-omega)*x[i] + omega*(b[i]-thing)/A.body[i][i]
        else:
        #
        #   in this case, the data is saved in a matrix-like structure
        #
            dummy = [w for w in x[k]]
            for i in range(len(b)):
                thing = [w for w in range(len(b))]
                thing.pop(i)
                thing = sum(
                        [A.body[i][j]*dummy[j] for j in thing]
                    )
                if omega==1:
                    dummy[i] = (b[i]-thing)/A.body[i][i]
                else:
                    dummy[i] = (1-omega)*dummy[i] + omega*(b[i]-thing)/A.body[i][i]
            x.append( dummy )
    #
    #
    #   iterate
    #
    #
    if While:
        Measure = Mes(A,x,b)
        k = 0
        while (k<20000 and Measure>MaxAbs(x)*1e-8):
            identical()
            k += 1
            Measure = Mes(A,x,b)
            if ShowProgress:
                print('k=' + str(k) + '; ' + str(Measure))
    else:
        for k in range(N):
            identical()
            if ShowProgress:
                if Forget:
                    Measure = Mes(A,x,b)
                else:
                    Measure = Mes(A,x[k],b)
                print('k=' + str(k) + '; ' + str(Measure))
    #
    #
    #   postamble
    #
    #
    if Plot:
        domain = [w for w in range(N+1)]
        fig, ax = plt.subplots()
        fig.set_tight_layout(True)
        for w in range(len(b)):
            ax.plot(domain, ExtractCol(x,w+1))
        ax.set(xlabel='Iterations', ylabel='Coordinate Values', title='Component-Wise Convergence of the Approximate Solution Vector')
        ax.grid()
        plt.show()
    if While:
        return {'soln' : x, 'times' : k}
    else:
        return x






#data = []
#for k in range(2,4):
if False:
    m = 2**k
    A = HomeworkMatrix(m)
    b = [1 for w in range(m**2)]
    yJ = Jacobi(A, b, naught=[2 for w in range(len(b))], ShowProgress=True)['times']
    print('')
    print('')
    print('Jacobi completed for m='+str(m) + '.')
    print('')
    print('')
    yGS = SOR(A, b, naught=[2 for w in range(len(b))], ShowProgress=True)['times']
    print('')
    print('')
    print('Gauss-Seidel completed for m='+str(m) + '.')
    print('')
    print('')
    ySOR = SOR(A, b, ShowProgress=True, naught=[2 for w in range(len(b))], omega=meg(m))['times']
    print('')
    print('')
    print('SOR completed for m='+str(m) + '.')
    print('')
    print('')
    data.append([yJ, yGS, ySOR])



A = HomeworkMatrix(15)
b = [1 for w in range(15**2)]
y = SOR(A, b, While=False, Forget=False, Plot=True, ShowProgress=False, lo=0, hi=20)





def SteepestDescent(A, b, naught=None, ShowProgress=True, tol=1e-8, lo=0, hi=0):
    #
    #   preamble
    #
    if not isinstance(A, Matrix):
        A = Matrix(A)
        print('Warning: A was attempted to be coerced to class Matrix.')
    if not A.body==t(A).body:
        raise ValueError('Gradient descent can only be used to solve Ax=b if A is strictly positive definite. The given matrix is not even symmetric.')
    if naught is None:
        naught = [ uniform(lo,hi) for w in range(len(b)) ]
    x = naught
    modr = 1e8
    #
    #   iterate
    #
    k = 0
    while (k<20000 and modr>tol*MaxAbs(x)):
        r = [                               # b-Ax
                b[i] - sum(A.body[i][j]*x[j] for j in range(len(x)))
            for i in range(len(b))
                ]
        modr = ModSqr(r)
        if modr==0:
            break
        alpha = modr / sum(                 # ||r||^2 / \lange r, Ar \rangle
                    r[i]*sum(A.body[i][j]*r[j] for j in range(len(r)))
                for i in range(len(b))
                    )
        for i in range(len(x)):             # x = x+alpha*r
            x[i] +=alpha*r[i]
        k += 1
        if ShowProgress:
            print('k=' + str(k) + '; ' + str(modr))
    #
    #   postamble
    #
    return {'soln' : x, 'times' : k}



def ConjugateGradient(A, b, naught=None, ShowProgress=True, tol=1e-8, lo=0, hi=0, check=False):
    #
    #   preamble
    #
    if not isinstance(A, Matrix):
        A = Matrix(A)
        print('Warning: A was attempted to be coerced to class Matrix.')
    if check:
        if not A.body==t(A).body:
            raise ValueError('The CG method can only be used to solve Ax=b if A is strictly positive definite. The given matrix is not even symmetric.')
    if naught is None:
        naught = [uniform(lo,hi) for w in range(len(b))]
    x = naught
    r = [       # b-Ax the initial residual
            b[i] - sum(A.body[i][j]*x[j] for j in range(len(x)))
        for i in range(len(b))
            ]
    v = [ j for j in r ]
    Measure = Mes(A,x,b)
    modr = ModSqr(r)
    #
    #   iterate
    #
    k = 0
    while (k<20000 and Measure>tol*MaxAbs(x)):
        Av = [
            sum(A.body[i][j]*v[j] for j in range(len(v)))
                for i in range(len(b))
            ]
        modv = sum( v[i]*Av[i] for i in range(len(v)) )
        alpha = float(modr)/modv
        for i in range(len(b)):
            x[i] += alpha*v[i]              # and x = x+alpha*v
            r[i] -= alpha*Av[i]             # r = b-alphaAv   (note r != b-Ax)
        modv = modr                         # a dumb attempt to conserve RAM
        modr = ModSqr(r)                    # ||b-alphaAv||^2, not ||residual||^2
        alpha = float(modr)/modv            # ||new.r||^2 / ||old.r||^2
        v = [ r[i] + alpha*v[i] for i in range(len(v)) ]
        k += 1
        Measure = Mes(A,x,b)
        if ShowProgress:
            print('k=' + str(k) + '; ' + str(Measure))
    #
    #   postamble
    #
    return {'soln' : x, 'times' : k}




data = []   # test how many iterations are necessary for the two method
for k in range(1,5):
    m = 2**k
    print('')
    print('')
    print('Working on m='+str(m) + '.')
    print('')
    print('')
    A = HomeworkMatrix(m)
    b = [1 for w in range(A.rows)]
    data.append([
        SteepestDescent(A, b, ShowProgress=True, lo=2, hi=2)['times'],
        ConjugateGradient(A, b, ShowProgress=True, lo=2, hi=2)['times'],
    ])





#
##
### ok, now back to numpy
##
#

from numpy import array as arr
from numpy import inf as inf
from numpy.linalg import norm as norm
from numpy import matmul as mul


def ConjugateGradient2(A, b, naught=None, ShowProgress=True, tol=1e-8, lo=0, hi=0, max_iter=20000):
    #
    #   preamble
    #
    if naught is None:
        naught = [uniform(lo,hi) for w in range(len(b))]
    if isinstance(A, Matrix):       # optional
        A = arr(A.body)             # optional
    x = arr([naught]).T
    b = arr([b]).T
    r = b - mul(A,x)
    v = arr([ j for j in r ])
    Measure = norm((mul(A,x)-b).T[0], ord=inf)
    modr = mul(r.T,r)[0][0]
    #
    #   iterate
    #
    k = 0
    while (k<max_iter and Measure>tol*max(x.max(),x.min())):
        Av = mul(A,v)
        modv = mul(v.T,Av)[0][0]
        alpha = float(modr)/modv
        x += alpha*v
        r -= alpha*Av
        modv = modr                         # a dumb attempt to conserve RAM
        modr = mul(r.T,r)[0][0]             # ||b-alphaAv||^2, not ||residual||^2
        alpha = float(modr)/modv            # ||new.r||^2 / ||old.r||^2
        v = r + alpha*v
        k += 1
        Measure = norm((mul(A,x)-b).T[0], ord=inf)
        if ShowProgress:
            print('k=' + str(k) + '; ' + str(Measure))
    #
    #   postamble
    #
    return {'soln' : x, 'times' : k}


data = []
for k in range(1,7):
    print('')
    print('')
    print('Working on k='+str(k) + '.')
    print('')
    print('')
    m = 2**k
    A = HomeworkMatrix(m)
    b = [1 for w in range(A.rows)]
    data.append( ConjugateGradient2( A, b, naught=[1 for w in range(A.cols)])['times'] )





