# -*- coding: utf-8 -*-
"""
Created on Thu May 23 18:51:02 2024

@author: ajrp
"""

import matplotlib.pyplot as plt
import numpy as np

a=2
b=1
q=0
d1=5
d2=40

def traza(a,b,q,d1,d2):
    return b-1-a**2-(d1+d2)*q**2

def det(a,b,q,d1,d2):
    return d1*d2*q**4+(d1*a**2+d2*(1-b))*q**2+a**2

def discr(a,b,q,d1,d2):
    return  traza(a, b, q, d1, d2)**2-4*det(a, b, q, d1, d2)


def f1(a,b,q,d1,d2):
    return (1.0/2.0)*(traza(a,b,q,d1,d2)+np.sqrt(traza(a, b, q, d1, d2)**2-4*det(a, b, q, d1, d2)))

def f2(a,b,q,d1,d2):
    return (1.0/2.0)*(traza(a,b,q,d1,d2)-np.sqrt(traza(a, b, q, d1, d2)**2-4*det(a, b, q, d1, d2)))

x1=np.linspace(0, 9)
x2=np.linspace(-6,6)

plt.figure()


plt.axhline(y = 0, color='black', linestyle='-', lw=1)
plt.axvline(x = 0, color='black', linestyle='-', lw=1)

plt.plot(traza(a, x1, q, d1, d2), det(a, x1, q, d1, d2) , label='q=0', zorder=1)

plt.plot(x2, x2**2/4 , color='b', label='$\Delta=0$', zorder=1)


plt.legend()
plt.title("")
plt.xlabel("Número de qubits")
plt.ylabel("Ciclos")

# Mostrar la gráfica
plt.show()
plt.savefig('fourier_cycles')