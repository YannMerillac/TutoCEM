# -*- coding: utf-8 -*-
from MoM import mom_interface
import numpy as np
import matplotlib.pyplot as plt

# n disc pts 
npts = 50
N = npts-2
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
a = 0.001*Lambda
# Antenna height
h = 0.25*Lambda
# discretization segments
w_z = np.linspace(-0.5*h,0.5*h,npts)

Z = np.zeros((N,N), dtype="complex",order='F')
mom_interface.build_impedance_matrix(w_z,k,a,Z)

i_source = int(npts/2-1)
V = np.zeros(N,dtype="complex")
mom_interface.build_sources(w_z,k,a,i_source,V)
V*=-1j*w*eps_0

I = np.linalg.solve(Z,V)
print(I)

plt.plot(w_z[1:-1],np.absolute(I))
plt.show()

Iin = I[i_source]
Z = 1./Iin
Zth = 20.*np.pi**2*(h/Lambda)**2
print("Iin = {} A, Zin= {} Ohms".format(Iin,Z))
print(Zth)


