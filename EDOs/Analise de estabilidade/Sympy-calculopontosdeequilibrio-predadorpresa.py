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
def dHdt(H,p,*params):
    return r*H - a*H*P

def dPdt(H,P,*params):
    return a*H*P - m*P

#Criar um Symbol para cada variável e parâmetro do sistema
H = sp.Symbol('H') 
P = sp.Symbol('P') 

r = sp.Symbol('r')
a = sp.Symbol('a')
m = sp.Symbol('m')

params = [r,a,m]

dHdt_equal0 = sp.Eq(dHdt(H,P,params),0)
dPdt_equal0 = sp.Eq(dPdt(H,P,params),0)

eqpoints = sp.solve((dHdt_equal0, dPdt_equal0),(H,P))
print("Pontos de equilíbrio: " + str(eqpoints))

nullcline_dHdt = sp.solve(dHdt_equal0,(H,P))
nullcline_dPdt = sp.solve(dPdt_equal0,(H,P))

print("Nullclines de dH/dt: " + str(nullcline_dHdt))
print("Nullclines de dP/dt: " + str(nullcline_dPdt))

eqMat = sp.Matrix([ dHdt(H,P,params), dPdt(H,P,params) ])
Mat = sp.Matrix([H, P ])
matrizJacobiana = eqMat.jacobian(Mat)
print('Matriz jacobiana %s' % matrizJacobiana)

trace_jac = matrizJacobiana.trace()
print('Traço\n' + str(matrizJacobiana.trace()))

det_jac = sp.simplify(matrizJacobiana.det())
print('Determinante\n' + str(matrizJacobiana.det()))