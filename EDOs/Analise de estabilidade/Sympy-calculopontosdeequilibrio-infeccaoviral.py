import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from sympy.utilities.lambdify import lambdify
from sympy import diff
from sympy import sqrt
import scipy.linalg as slin

#Referências sobre sympy
#https://docs.sympy.org/latest/index.html
#https://pythonforundergradengineers.com/sympy-expressions-and-equations.html

# EDOs do sistema
def dEdt(E,I,V,T,*params):
    return -b*E*V + r*E

def dIdt(E,I,V,T,*params):
    return b*E*V - k*I*T - a*I

def dVdt(E,I,V,T,*params):
    return p*a*I - k2*V 

def dTdt(E,I,V,T,*params):
    return s*I - m*T

#Criar um Symbol para cada variável e parâmetro do sistema
E = sp.Symbol('E') 
I = sp.Symbol('I') 
V = sp.Symbol('V') 
T = sp.Symbol('T') 

b = sp.Symbol('b')
k = sp.Symbol('k')
a = sp.Symbol('a')
p = sp.Symbol('p')
k2 = sp.Symbol('k2')
s = sp.Symbol('s')
m = sp.Symbol('m')
r = sp.Symbol('r')

params = [b,k,a,p,k2,s,m,r]

dEdt_equal0 = sp.Eq(dEdt(E,I,V,T,params),0)
dIdt_equal0 = sp.Eq(dIdt(E,I,V,T,params),0)
dVdt_equal0 = sp.Eq(dVdt(E,I,V,T,params),0)
dTdt_equal0 = sp.Eq(dTdt(E,I,V,T,params),0)

eqpoints = sp.solve((dEdt_equal0, dIdt_equal0, dVdt_equal0, dTdt_equal0),(E,I,V,T))
print("Pontos de equilíbrio: " + str(eqpoints))

nullcline_dEdt = sp.solve(dEdt_equal0,(E,I,V,T))
nullcline_dIdt = sp.solve(dIdt_equal0,(E,I,V,T))
nullcline_dVdt = sp.solve(dVdt_equal0,(E,I,V,T))
nullcline_dTdt = sp.solve(dTdt_equal0,(E,I,V,T))

print("Nullclines de dE/dt: " + str(nullcline_dEdt))
print("Nullclines de dI/dt: " + str(nullcline_dIdt))
print("Nullclines de dV/dt: " + str(nullcline_dVdt))
print("Nullclines de dT/dt: " + str(nullcline_dTdt))

eqMat = sp.Matrix([ dEdt(E,I,V,T,params), dIdt(E,I,V,T,params), dVdt(E,I,V,T,params), dTdt(E,I,V,T,params) ])
Mat = sp.Matrix([ E, I, V, T ])
matrizJacobiana = eqMat.jacobian(Mat)
print('Matriz jacobiana %s' % matrizJacobiana)

trace_jac = matrizJacobiana.trace()
print('Traço\n' + str(matrizJacobiana.trace()))

#mat = np.array(matrizJacobiana)
det_jac = matrizJacobiana.det()
print('Determinante\n' + str(matrizJacobiana.det()))