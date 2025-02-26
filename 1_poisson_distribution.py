#!/usr/bin/env python
import numpy as np

# Using log conversion
def poisson1(lamda:np.float32, k:np.float32) -> np.float32:
    """Computing the poisson distribution value for parameters (lamda, k) using log conversion."""
    j = np.float32(0)
    for i in range(np.int32(k)):
        i = np.float32(i)
        j += np.float32(np.log(k-i))
    P = np.float32(np.exp(k*np.log(lamda)-lamda-j))
    return P

# Using different order of operations
def poisson2(lamda:np.float32, k:np.float32) -> np.float32:
    """Computing the poisson distribution value for parameters (lamda, k) using an unconventional order of operations."""
    P = np.float32(np.exp(-lamda))
    for i in range(np.int32(k)):
        i = np.float32(i)
        P *= np.float32(lamda/(k-i))
    return P

# Using different order of operations and log conversion, avoiding np.exp(-lamda)
def poisson3(lamda:np.float32, k:np.float32) -> np.float32:
    """Computing the poisson distribution value for parameters (lamda, k) using a different order of operations and a log conversion."""
    P = np.float32(1)
    for i in range(np.int32(k)):
        i = np.float32(i)
        P *= np.float32(lamda/(k-i))
    P = np.float32(np.exp(np.log(P) - lamda))
    return P

# Computing and printing results
lamda_k = np.array([[1,0],[5,10],[3,21],[2.6,40],[100,5],[101,200]], dtype=np.float32)
print("[lamda k] P1 P2 P3")
for values in lamda_k:
    P1 = poisson1(values[0], values[1])
    P2 = poisson2(values[0], values[1])
    P3 = poisson3(values[0], values[1])
    print(f"[{values[0]:.1f} {values[1]:.0f}] {P1:.6} {P2:.6} {P3:.6}")