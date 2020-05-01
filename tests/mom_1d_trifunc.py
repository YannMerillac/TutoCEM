# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad,dblquad

def G(z,zp,a):
    R = np.sqrt((z-zp)**2+a**2)
    return np.exp(-k*R*1j)/(4.*np.pi*R)

def tri_func(zp,ctrl_z,i):
    z1 = ctrl_z[i-1]
    z2 = ctrl_z[i]
    z3 = ctrl_z[i+1]
    if zp>=z1 and zp<=z2:
        return (zp-z1)/(z2-z1)
    elif zp>=z2 and zp<=z3:
        return (z3-zp)/(z3-z2)
    else:
        return 0.
    
def tri_func_grad(zp,ctrl_z,i):
    z1 = ctrl_z[i-1]
    z2 = ctrl_z[i]
    z3 = ctrl_z[i+1]
    if zp>=z1 and zp<=z2:
        return 1./(z2-z1)
    elif zp>=z2 and zp<=z3:
        return -1./(z3-z2)
    else:
        return 0.
    
def I_A(m,n,w_z,a):
    z_bnd = (w_z[m-1],w_z[m+1])
    zp_bnd = (w_z[n-1],w_z[n+1])
    
    def I_n(z):
        integrand_n = lambda zp: G(z,zp,a)*tri_func(zp,w_z,n)
        integrand_n_R = lambda zp : np.real(integrand_n(zp))
        integrand_n_I = lambda zp : np.imag(integrand_n(zp))
        I_n_R = quad(integrand_n_R, zp_bnd[0], zp_bnd[1])[0]
        I_n_I = quad(integrand_n_I, zp_bnd[0], zp_bnd[1])[0]
        return I_n_R+1j*I_n_I
    
    def real_func(z):
        return np.real(I_n(z)*tri_func(z,w_z,m))
    
    def imag_func(z):
        return np.imag(I_n(z)*tri_func(z,w_z,m))
    
    real_I = quad(real_func, z_bnd[0], z_bnd[1])[0]
    imag_I = quad(imag_func, z_bnd[0], z_bnd[1])[0]
    return real_I+1j*imag_I
    

def I_Phi(m,n,w_z,a):
    z_bnd = (w_z[m-1],w_z[m+1])
    zp_bnd = (w_z[n-1],w_z[n+1])
    
    def I_n(z):
        integrand_n = lambda zp: G(z,zp,a)*tri_func_grad(zp,w_z,n)
        integrand_n_R = lambda zp : np.real(integrand_n(zp))
        integrand_n_I = lambda zp : np.imag(integrand_n(zp))
        I_n_R = quad(integrand_n_R, zp_bnd[0], zp_bnd[1])[0]
        I_n_I = quad(integrand_n_I, zp_bnd[0], zp_bnd[1])[0]
        return I_n_R+1j*I_n_I
    
    def real_func(z):
        return np.real(I_n(z)*tri_func_grad(z,w_z,m))
    
    def imag_func(z):
        return np.imag(I_n(z)*tri_func_grad(z,w_z,m))
    
    real_I = quad(real_func, z_bnd[0], z_bnd[1])[0]
    imag_I = quad(imag_func, z_bnd[0], z_bnd[1])[0]
    return real_I+1j*imag_I

def I_b(m,w_z,i_source):
    z_bnd = (w_z[m-1],w_z[m+1])
    z1_source = w_z[i_source]
    z2_source = w_z[i_source+1]
    print(z1_source,z2_source)
    dz = z2_source-z1_source
    def I(z):
        if z>=z1_source and z<=z2_source:
            Ei = 1./dz
            return Ei*tri_func(z, w_z, m)
        else:
            return 0.
    return quad(I, z_bnd[0], z_bnd[1])[0]
    
# n disc pts 
npts = 20
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
h = 0.1*Lambda
# discretization segments
w_z = np.linspace(-0.5*h,0.5*h,npts)

c_z = 0.5*(w_z[1:]+w_z[:-1])

z_disc = np.linspace(-0.5*h,0.5*h,1000)
f_disc = [tri_func(z,w_z,2) for z in z_disc]
g_disc = [tri_func_grad(z,w_z,2) for z in z_disc]
plt.plot(z_disc,f_disc)
plt.plot(z_disc,g_disc)
plt.show()

i_source = int(npts/2-1)
print(i_source)

Z = np.zeros((N,N), dtype="complex")
V = np.zeros(N,dtype="complex")
for m in range(N):
    for n in range(N):
        Z[m,n] = (k**2)*I_A(m+1,n+1,w_z, a)-I_Phi(m+1,n+1,w_z, a)
    V[m] = -1j*w*eps_0*I_b(m+1,w_z,i_source)
        
#print(Z)
print(V)

I = np.linalg.solve(Z,V)
print(I)

plt.plot(w_z[1:-1],np.absolute(I))
plt.show()

Iin = I[i_source]
Z = 1./Iin
Zth = 20.*np.pi**2*(h/Lambda)**2
print("Iin = {} A, Zin= {} Ohms".format(Iin,Z))
print(Zth)
