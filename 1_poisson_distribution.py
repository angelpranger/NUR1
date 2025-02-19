import numpy as np

# overflow caused by the factorial in the traditional way for large k

# using log conversion
def poisson1(lamda, k):
    j = np.float32(0)
    for i in range(np.int32(k)):
        j += np.log(k-i)
    P = np.float32(np.exp(k*np.log(lamda)-lamda-j))
    return P

# using different order of operations
def poisson2(lamda, k):
    P = np.float32(np.exp(-lamda))
    for i in range(np.int32(k)):
        P *= lamda/(k-i)
    return P

# using different order of operations
def poisson3(lamda, k):
    P = np.float32(1)
    for i in range(np.int32(k)):
        P *= lamda/(k-i)
    P = np.exp(np.log(P) - lamda)
    return P

lamda_k = np.array([[1,0],[5,10],[3,21],[2.6,40],[100,5],[101,200]], dtype=np.float32)
for values in lamda_k:
    print(np.exp(-values[0]))
    P1 = poisson1(values[0], values[1])
    P2 = poisson2(values[0], values[1])
    P3 = poisson3(values[0], values[1])
    print(P1, P2, P3)

# TODO: why do the results for the last two values differ? seems like the exponent gives a large round off error
# version 3 is the most accurate according to wolframalpha, but why? Im not sure