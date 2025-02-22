#This script is to get you started with reading the data and plotting it
#You are free to change whatever you like/do it completely differently,
#as long as the results are clearly presented

import numpy as np
import sys
import os
import matplotlib.pyplot as plt


# gecheckt, werkt voor A, geeft geen nans meer voor V
def LU_decomposition(matrix:np.ndarray):
    # Assume alpha_ii = 1
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


# TODO Neville's algorithm

def bisection(M:np.int64, x_samp:np.array, x_interp:np.array) -> np.ndarray:
# Assumes x_samp are uniformly spaced
    if (M > x_samp.shape[0]):
        print("M is too large for the array x_samp")
        return False
    for (i,x) in enumerate(x_interp):
        l_idx = np.int64(0)
        r_idx = np.int64(x_samp.shape[0]-1)
        while (r_idx - l_idx > 1):
            m_idx = np.int64((l_idx + r_idx)*0.5)
            if (x <= x_samp[m_idx]):
                r_idx = m_idx
            elif (x > x_samp[m_idx]):
                l_idx = m_idx
        # Computing the left boundary index (ignoring range)
        l_idx -= np.int64(np.floor(0.5*M))
        # Checking range
        if l_idx < 0:
            l_idx = 0
        if (l_idx + M >= x_samp.shape[0]):
            l_idx -= l_idx + M - x_samp.shape[0]
        x_interp[i] = l_idx
    return np.int64(x_interp)

def neville(M, x_samp, y_samp, x_interp):
    # Finding left index of M tabulated points using bisection
    l_indices = bisection(M, x_samp.copy(), x_interp.copy())
    #y_samp = y_samp[np.newaxis,:]
    # Initializing matrix P
    Ps = np.zeros((x_interp.shape[0], M), dtype=np.float64)
    for i in range(M):
        Ps[:,i] = y_samp[l_indices+i]
    print(Ps)
    # Looping over orders k and current intervals
    for k in range(1, M):
        for i in range(M-1-k):
            j = i+k
            Ps[:,i] = ((x_samp[l_indices+j]-x_interp)*Ps[:,i]+(x_interp-x_samp[l_indices+i])*Ps[:,j])/(x_samp[l_indices+j]-x_samp[l_indices+i])
    error_estimate = np.absolute(Ps[:,0]-Ps[:,1])
    return Ps[:,0], error_estimate

# TODO something wrong with neville (bisection seems fine)
# For lower M, result is better than for larger M (shouldnt be right)

data=np.genfromtxt(os.path.join(sys.path[0],"Vandermonde.txt"),comments='#',dtype=np.float64)
x=data[:,0]
y=data[:,1]
xx=np.linspace(x[0],x[-1],1001) #x values to interpolate at

#Insert your own code to calculate interpolated y values here!

# 2a

# Constructing V
V = np.zeros((x.shape[0], x.shape[0]), dtype=np.float64)
for j in range(V.shape[0]):
    V[:,j] = np.power(x,j)
#print(V) this is correct TODO

# check, LU is correct
# A = np.array([[1,2,4],[3,8,14],[2,6,13]])
# print(A)
# print(LU_decomposition(A))

# LU-decomposing V
LU_V = LU_decomposition(V.copy())
#print(LU_V)

# Solving for c
c = solving_system_with_LU(LU_V, y.copy())
#print(c)

# Constructing the polynomial y
def polynomial(c, x):
    y = np.float64(0)
    for j in range(c.shape[0]):
        y += c[j] * np.power(x,j)
    return np.float64(y)

# 2b

# 3b


#The plotting code below assumes you've given the interpolated
#values for 2a suffix "a", those for 2b "b", and those for 2c 
#"c1" and "c10" – feel free to change these
#Note that you interpolate to xx for the top panel but to
#x for the bottom panel (since we can only compare the accuracy
#to the given points, not in between them)
# yya=np.zeros(1001,dtype=np.float64) #replace!
# ya=np.zeros(len(x),dtype=np.float64) #replace!
yya = polynomial(c, xx)
ya = polynomial(c, x)
yyb = neville(3, x, y, xx)[0]
yb = neville(3, x, y, x)[0]
yyc1=yya #replace!
yc1=ya #replace!
yyc10=yya #replace!
yc10=ya #replace!

#Don't forget to output the coefficients you find with your LU routine

#Plot of points with absolute difference shown on a log scale (question 2a)
fig=plt.figure()
gs=fig.add_gridspec(2,hspace=0,height_ratios=[2.0,1.0])
axs=gs.subplots(sharex=True,sharey=False)
axs[0].plot(x,y,marker='o',linewidth=0)
plt.xlim(-1,101)
axs[0].set_ylim(-400,400)
axs[0].set_ylabel('$y$')
axs[1].set_ylim(1e-16,1e1)
axs[1].set_ylabel('$|y-y_i|$')
axs[1].set_xlabel('$x$')
axs[1].set_yscale('log')
line,=axs[0].plot(xx,yya,color='orange')
line.set_label('Via LU decomposition')
axs[0].legend(frameon=False,loc="lower left")
axs[1].plot(x,abs(y-ya),color='orange')
plt.savefig('my_vandermonde_sol_2a.png',dpi=600)

#For questions 2b and 2c, add this block
line,=axs[0].plot(xx,yyb,linestyle='dashed',color='green')
line.set_label('Via Neville\'s algorithm')
axs[0].legend(frameon=False,loc="lower left")
axs[1].plot(x,abs(y-yb),linestyle='dashed',color='green')
plt.savefig('my_vandermonde_sol_2b.png',dpi=600)

#For question 2c, add this block too
line,=axs[0].plot(xx,yyc1,linestyle='dotted',color='red')
line.set_label('LU with 1 iteration')
axs[1].plot(x,abs(y-yc1),linestyle='dotted',color='red')
line,=axs[0].plot(xx,yyc10,linestyle='dashdot',color='purple')
line.set_label('LU with 10 iterations')
axs[1].plot(x,abs(y-yc10),linestyle='dashdot',color='purple')
axs[0].legend(frameon=False,loc="lower left")
plt.savefig('my_vandermonde_sol_2c.png',dpi=600)

#Don't forget to caption your figures to describe them/
#mention what conclusions you draw from them!