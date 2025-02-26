import numpy as np
import timeit
import sys
import os
import matplotlib.pyplot as plt

# Importing data
data=np.genfromtxt(os.path.join(sys.path[0],"Vandermonde.txt"),comments='#',dtype=np.float64)
x=data[:,0]
y=data[:,1]
xx=np.linspace(x[0],x[-1],1001) # x values to interpolate at

# Functions used
def LU_decomposition(matrix:np.ndarray):
    """Takes a matrix and computes its LU-decomposition.
    Returns the same matrix which now holds the LU-decomposition: the upper triangle and diagonal hold the beta_ij (U),
    while the lower triangle hold the alpha_ij (L). The alpha_ii are chosen to be equal to 1."""
    # Choose alpha_ii = 1
    # Loop over the columns j
    for j in range(matrix.shape[0]):
        # Keep beta_0j = a_0j
        # Update beta_ij where 0<i<=j
        for i in range(1, j+1):
            term = np.float64(0)
            for k in range(0, i):
                term += matrix[i,k]*matrix[k,j]
            matrix[i,j] -= term
        # Update alpha_ij where i>j
        beta_jj = 1/matrix[j,j]
        for i in range(j+1, matrix.shape[0]):
            term = np.float64(0)
            for k in range(0, j):
                term += matrix[i,k]*matrix[k,j]
            matrix[i,j] = (matrix[i,j]-term)*beta_jj
    return np.float64(matrix)

def solving_system_with_LU(LU:np.ndarray, b:np.array) -> np.array:
    """Solves the system Ax=b by taking the LU decomposition of A and b as input and applying forward and backward substitution.
    Returns the vector b which now holds the solution x."""
    # Knowing that alpha_ii = 1
    # Forward substitution
    # Keep y_0 = b_0 as alpha_00 = 1
    for i in range(1, LU.shape[0]):
        term = np.float64(0)
        for j in range(i):
            term += LU[i,j]*b[j]
        b[i] -= term
    # Backward substitution
    b[-1] /= LU[-1,-1]
    for i in range(LU.shape[0]-1):
        i = LU.shape[0]-1 - i
        term = np.float64(0)
        for j in range(i+1, LU.shape[0]):
            term += LU[i,j]*b[j]
        b[i] = (b[i]-term)/LU[i,i]
    return np.float64(b)

def bisection(M:np.int64, x_samp:np.array, x_interp:np.array) -> np.ndarray:
    """Takes as input a number of points M, a sample array x and an interpolation values array x_interp.
    Uses bisection to find the closest one or two points in x_samp for each point in x_interp. Subsequently finds the closest M
    points in x_samp and returns the index in x_samp of the lowest (leftmost) point."""
    # Assumes x_samp are uniformly spaced
    if (M > x_samp.shape[0]):
        print("M is too large for the array x_samp")
        return
    # Loop over all points x in x_interp
    for (i,x) in enumerate(x_interp):
        # Set left and right edge indices within x_samp
        l_idx = np.int64(0)
        r_idx = np.int64(x_samp.shape[0]-1)
        # Half the searching area, check in which half x lies, update edge indices for this half,
        # until the left and right edge index differ by 0 or 1
        while (r_idx - l_idx > 1):
            m_idx = np.int64(np.floor((l_idx + r_idx)*0.5))
            if (x <= x_samp[m_idx]):
                r_idx = m_idx
            elif (x > x_samp[m_idx]):
                l_idx = m_idx
        # Save the middle value of the current interval to later check whether x is closer to the left or right edge
        m = (x_samp[l_idx] + x_samp[r_idx])*0.5
        # Computing the left boundary index (ignoring range)
        l_idx -= np.int64(np.ceil(0.5*M-1))
        # Correcting left boundary if x is closer to the right boundary 
        # (in which case the right boundary is seen as the first point instead of the second)
        if ((x > m) & (M % 2 == 1)):
            l_idx += np.float64(1)
        # Checking range (ensuring the left and right edge do not fall outside x_samp)
        if l_idx < 0:
            l_idx = 0
        if (l_idx + M >= x_samp.shape[0]):
            l_idx -= l_idx + M - x_samp.shape[0]
        x_interp[i] = l_idx
    return np.int64(x_interp)

def neville(M:np.int64, x_samp:np.array, y_samp:np.array, x_interp:np.array) -> np.ndarray:
    """Takes an integer number of points M, sample array x_samp with corresponding y values y_samp and interpolation values x_interp.
    Computes for each x in x_interp the value of the unique Lagrange polynomial through the M points in x_samp closest to x, using Neville's algorithm.
    Returns the M-1 order Lagrange polynomial values corresponding to the x in x_interp and the corresponding error estimates."""
    if (M < 1):
        print("M is too small")
        return
    # Finding left index of M tabulated points using bisection
    l_indices = bisection(M, x_samp.copy(), x_interp.copy())
    # Initializing matrix P for Neville's algorithm (each row corresponds to a x in x_interp, each column corresponds to one of the M closest sample points)
    Ps = np.zeros((x_interp.shape[0], M), dtype=np.float64)
    for i in range(M):
        Ps[:,i] = y_samp[l_indices+i]
    # Looping over orders k and current intervals
    for k in range(1, M):
        for i in range(M-k):
            j = i+k # note this j is only for the x_samp, as Ps is overwritten so just j=i+1 as index
            # Ps[:,i] are overwritten to hold the linear interpolation between the previous values
            Ps[:,i] = ((x_samp[l_indices+j]-x_interp)*Ps[:,i]+(x_interp-x_samp[l_indices+i])*Ps[:,i+1])/(x_samp[l_indices+j]-x_samp[l_indices+i])
    if (M < 2): # No error estimate
        return Ps[:,0], 0
    error_estimate = np.absolute(Ps[:,0]-Ps[:,1])
    return Ps[:,0], error_estimate

def polynomial(c:np.array, x:np.array) -> np.array:
    """Constructing the values y of the polynomial with coefficients given by c at points x.
    Returns y values."""
    y = np.float64(0)
    for j in range(c.shape[0]):
        y += c[j] * np.power(x,j)
    return np.float64(y)

# 2a

def q2a(x:np.array, y:np.array, xx:np.array) -> np.ndarray:
    """Constructing Vandermonde matrix V from data x, LU-decomposing V and solving Vc=y for c.
    Returns c and the values of the polynomial, yya and ya, corresponding to xx and x."""
    # Constructing V
    V = np.zeros((x.shape[0], x.shape[0]), dtype=np.float64)
    for j in range(V.shape[0]):
        V[:,j] = np.power(x,j)
    # LU-decomposing V
    LU_V = LU_decomposition(V.copy())
    # Solving for c
    c0 = solving_system_with_LU(LU_V, y.copy())
    # Finding the polynomial values y
    yya = polynomial(c0, xx)
    ya = polynomial(c0, x)
    return V, LU_V, c0, yya, ya

V, LU_V, c0, yya, ya = q2a(x, y, xx)

# Printing the values of c
print("The coefficients of the polynomial found with LU-decomposition are")
print(c0)

# Plot of points with absolute difference shown on a log scale (question 2a)
fig=plt.figure()
gs=fig.add_gridspec(2,hspace=0,height_ratios=[2.0,1.0])
axs=gs.subplots(sharex=True,sharey=False)
axs[0].plot(x,y,marker='o',linewidth=0)
plt.xlim(-1,101)
axs[0].set_ylim(-400,600)
axs[0].set_ylabel('$y$')
axs[1].set_ylim(1e-18,1e4)
axs[1].set_ylabel('$|y-y_i|$')
axs[1].set_xlabel('$x$')
axs[1].set_yscale('log')
line,=axs[0].plot(xx,yya,color='orange')
line.set_label('Via LU decomposition')
axs[0].legend(frameon=False,loc="lower left")
axs[1].plot(x,abs(y-ya),color='orange')
plt.savefig('my_vandermonde_sol_2a.png',dpi=600)

# 2b

def q2b(x:np.array, y:np.array, xx:np.array, M:np.int64):
    """Computing the lagrange polynomial through the given data points (x,y) using Neville's algorithm (iff M=x.shape[0]).
    Returns values yyb and yb of the polynomial corresponding to the points in xx and x respectively."""
    yyb = neville(M, x, y, xx)[0]
    yb = neville(M, x, y, x)[0]
    return yyb, yb

yyb, yb = q2b(x, y, xx, 20)

# For questions 2b and 2c, add this block
line,=axs[0].plot(xx,yyb,linestyle='dashed',color='green')
line.set_label('Via Neville\'s algorithm')
axs[0].legend(frameon=False,loc="lower left")
axs[1].plot(x,abs(y-yb),linestyle='dashed',color='green')
plt.savefig('my_vandermonde_sol_2b.png',dpi=600)

# 2c

def q2c(LU_V:np.ndarray, V:np.ndarray, c:np.array, y:np.array, x:np.array, xx:np.array, num:np.int64) -> np.ndarray:
    """Computing num iterative improvements on the solution c of the system Vc=y solved by LU-decomposition.
    Returns values yyc and yc of the polynomial corresponding to the points in xx and x respectively."""
    # Iterative improvements
    for i in range(num):
        c -= solving_system_with_LU(LU_V, (np.matmul(V,c)-y).copy())
    # Computing polynomial values
    yyc = polynomial(c, xx)
    yc = polynomial(c, x)
    return yyc, yc

yyc1, yc1 = q2c(LU_V, V, c0.copy(), y, x, xx, 1)
yyc10, yc10 = q2c(LU_V, V, c0.copy(), y, x, xx, 10)

# For question 2c, add this block too
line,=axs[0].plot(xx,yyc1,linestyle='dotted',color='red')
line.set_label('LU with 1 iteration')
axs[1].plot(x,abs(y-yc1),linestyle='dotted',color='red')
line,=axs[0].plot(xx,yyc10,linestyle='dashdot',color='purple')
line.set_label('LU with 10 iterations')
axs[1].plot(x,abs(y-yc10),linestyle='dashdot',color='purple')
axs[0].legend(frameon=False,loc="lower left")
plt.savefig('my_vandermonde_sol_2c.png',dpi=600)

# 2d

# Timing the previous exercises
num = np.int64(1000)
print(f"Running 2a {num} times takes", timeit.timeit(lambda: q2a(x,y,xx), number=num), "seconds.")
print(f"Running 2b {num} times takes", timeit.timeit(lambda: q2b(x,y,xx,20), number=num), "seconds.")
print(f"Running 2c {num} times takes", timeit.timeit(lambda: q2c(LU_V,V,c0.copy(),y,x,xx,10), number=num), "seconds.")