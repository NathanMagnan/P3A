##Imports
import math as ma
import sympy as s
import scipy.integrate as sc
import scipy.optimize as so
import numpy as np

## Abaques


def Kp():
    def k(b,m):
        return 0.5*(ma.sin(b))**2*(1-m*(ma.sin(b))**2)**(-1.5)


    m=np.linspace(0,1,100000,endpoint=False)

    def K_p(m):
        if (m>= 1):
            return "Error can't return complexe value"
        return sc.quad(k,0,ma.pi/2,args=(m))[0]

    l=list()
    for i in range (len(m)):
        l.append(K_p(m[i]))


    return l



def Ep():
    def e(b,m):
        return -0.5*(ma.sin(b))**2*(1-m*(ma.sin(b))**2)**(-0.5)

    m=np.linspace(0,1,100000,endpoint=False)

    def E_p(m):
        if (m>= 1):
            return "Error can't return complexe value"
        return sc.quad(e,0,ma.pi/2,args=(m))[0]

    l=list()
    for i in range (len(m)):
        l.append(E_p(m[i]))

    return l

m=np.linspace(0,1,100000,endpoint=False)

def save(m,Kp,Ep):
    np.savetxt("Abacus_m",m)
    np.savetxt("Abacus_Kp",Kp)
    np.savetxt("Abacus_Ep",Ep)
    return

save(m,Kp(),Ep())
