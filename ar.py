import numpy as np
import matplotlib.pyplot as plt
import time
import array
import cmath
import math

#MZ interferometer parametrs
delta = 0
k = 0.005 * np.pi* pow(10, (6))         # in meter^-1
ldc = 50 * pow(10,(-6))                 # in meter
larm = 0 * pow(10, (-6))                # in meter
q = 2
beta = 5


#angle phi in radian (0 - 2 pi)
phi = np.linspace(0,2*np.pi,num=100)

# 2 x 2 T_DC matrix elements
a11 =       (np.cos(q*ldc) +1J*(delta/q) * np.cos(q*ldc)) * np.exp(-1J*delta*ldc)
a12 = -1J *(k/q) *np.sin(q*ldc) * np.exp(-1J*delta*ldc)
a21 = -1J *(k/q) *np.sin(q*ldc) * np.exp(1J*delta*ldc)
a22 =       (np.cos(q*ldc) - 1J*(delta/q) * np.cos(q*ldc)) * np.exp(1J*delta*ldc)


# T_DC matrix
Tdc = np.array([[a11, a12], [a21, a22]])

## 2 x 2 T_arm matrix elements
b11 = np.exp(-1J*beta*larm)
b12 = 0
b21 = 0
b22 = np.exp(-1J*beta*larm)

# T_arm matrix
Tarm = np.array([[b11,b12],[b21,b22]])

# T = T_DC * T_arm * T_DC

temp1 = np.dot(Tarm,Tdc) 
T = np.dot(Tdc, temp1) 

#Output model power A_{1}, B_{1}
A = np.zeros(len(phi), dtype=np.complex)
B = np.zeros(len(phi), dtype=np.complex)

ind1 = 0 
for ind2 in phi:
#Single drive T_{arm}
    Tarm = np.array([[1,0],[0,np.exp(1J*ind2)]])
#Push-pull drive T_{arm}
    Tarm = np.array([[np.exp(1J*ind2),0],[0,np.exp(-1J*ind2)]])
#Total system matrix
    temp1 = np.dot(Tarm,Tdc) 
    T = np.dot(Tdc, temp1) 
# A_{1}, B_{1}
    A[ind1] = T[0,0]
    B[ind1]=T[1,0]
    ind1+=1
    
A2 = abs(A)*abs(A)
B2 = abs(B)*abs(B)


argA = A
argB = B

for i in range(A.shape[0]):
    argA[i] = math.atan(A[i])
    argB[i] = math.atan(B[i])


plt.figure()
plt.plot(phi,A2)  
plt.title('Model power of A for MZ interferometer' )
plt.xlabel('phi in radian (0-2 pi) -->');
plt.ylabel('A1 ^ 2 -->');

plt.figure()
plt.plot(phi,B2)
plt.title('Model power of B for MZ interferometer' )
plt.xlabel('phi in radian (0-2 pi) -->');
plt.ylabel('B1 ^ 2 -->');
plt.figure()

plt.plot(phi,argA)
plt.title('arg of A for MZ interferometer' )
plt.xlabel('phi in radian (0-2 pi) -->');
plt.ylabel('arg(A) -->');

plt.figure()
plt.plot(phi,argB)
plt.title('arg of B for MZ interferometer' )
plt.xlabel('phi in radian (0-2 pi) -->');
plt.ylabel('arg(B) -->');







