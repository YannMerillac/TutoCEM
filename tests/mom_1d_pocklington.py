# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def G(z,zp,k,a):
    R = np.sqrt((z-zp)**2+a**2)
    return np.exp(-k*R*1j)/(4.*np.pi*R)

def I_K(m,n,w_z,k,a):
    zm1 = w_z[m]
    zm2 = w_z[m+1]
    z = 0.5*(zm1+zm2)
    zn1 = w_z[n]
    zn2 = w_z[n+1]
    # Integrand function
    def I(zp):
        R = np.sqrt((z-zp)**2+a**2)
        Ke = np.exp(-k*R*1j)/(4.*np.pi*R**5)
        Ke*=((1.+1j*k*R)*(2.*R**2-3.*a**2)+(k*a*R)**2)
        return Ke
    def real_func(zp):
        return np.real(I(zp))
    def imag_func(zp):
        return np.imag(I(zp))
    real_I = quad(real_func, zn1, zn2)[0]
    imag_I = quad(imag_func, zn1, zn2)[0]
    return real_I+1j*imag_I
    
# n disc pts 
npts = 10
N = npts-1
# physical properties
eps_0 = 8.8542e-12
mu_0 = 4.0e-7*np.pi
c = np.sqrt(1/(eps_0*mu_0))
# frequency = 100MHz
#f = 10.0e6
Lambda = 1. #c/f
f = c/Lambda
w = 2*np.pi*f
k = 2.*np.pi/Lambda
print("C = {} m/s, Lambda = {} m".format(c,Lambda))
# wire radius
a = 0.005*Lambda
# Antenna height
h = 0.1*Lambda
# discretization segments
w_z = np.linspace(-0.5*h,0.5*h,npts)

c_z = 0.5*(w_z[1:]+w_z[:-1])

# build impedance matrix
Z = np.zeros((N,N), dtype="complex")
for m in range(N):
    for n in range(N):
        Z[m,n] = I_K(m,n,w_z,k,a)

# build RHS
V = np.zeros(N, dtype="complex")
i_source = int(npts/2-1)
print(i_source)
dz = w_z[i_source+1]-w_z[i_source]
V[i_source] = -1j*w*eps_0/dz

I = np.linalg.solve(Z,V)
print(I)

I_abs = np.absolute(I)
plt.plot(c_z,I_abs)
plt.show()

Iin = I[i_source]
Z = 1./Iin
Zth = 20.*np.pi**2*(h/Lambda)**2
print("Iin = {} A, Zin= {} Ohms".format(Iin,Z))
print(Zth)