import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from sympy import diff
import scipy.linalg as slin

from sympy.solvers.inequalities import solve_poly_inequality
from sympy.solvers.inequalities import reduce_rational_inequalities
from sympy import Poly

n1 = sp.Symbol('n1')
n2 = sp.Symbol('n2')
r1 = sp.Symbol('r1')
r2 = sp.Symbol('r2')
w11 = sp.Symbol('w11')
w12 = sp.Symbol('w12')
w22 = sp.Symbol('w22')
w21 = sp.Symbol('w21')

symbols = [n1,n2,r1,r2,w11,w12,w22,w21]

def dn1dt(u,v):
    return r1*n1*(1 - w11*n1 - w21*n2)

def dn2dt(u,v):
    return r2*n2*(1 - w22*n2 - w12*n1)

dn1dt_equal0 = sp.Eq(dn1dt(n1,n2),0)
dn2dt_equal0 = sp.Eq(dn2dt(n1,n2),0)

eqpoints = sp.solve((dn1dt_equal0, dn2dt_equal0),n1,n2)
print("Eq points: " + str(eqpoints))

nullcline_dn1dt = sp.solve(dn1dt_equal0,n1,n2)
nullcline_dn2dt = sp.solve(dn2dt_equal0,n1,n2)

eqMat = sp.Matrix([ dn1dt(n1,n2), dn2dt(n1,n2) ])
Mat = sp.Matrix([ n1, n2 ])
matrizJacobiana = eqMat.jacobian(Mat)
#print('Matriz jacobiana %s' % matrizJacobiana)

trace_jac = matrizJacobiana.trace()
trace_jac = sp.simplify(trace_jac)
print('Traço: ' + str(trace_jac) + '\n')

#mat = np.array(matrizJacobiana)
det_jac = matrizJacobiana.det()
det_jac = sp.simplify(det_jac)
print('Determinante:' + str(det_jac) + '\n')

alltrs = []
alldets = []
for eqp in eqpoints:
    tr = sp.simplify(trace_jac.subs([(n1,eqp[0]),(n2,eqp[1])]))
    print('trace: '+str(tr))
    #res = sp.solveset(tr < 0,r1,S.Reals)
    #print(res)
    alltrs.append(tr)

    det = sp.simplify(det_jac.subs([(n1,eqp[0]),(n2,eqp[1])]))
    print('determinant: '+str(det))
    alldets.append(det)

print('\n')
print('Condições de estabilidade para o 1º ponto de equilíbrio: ')
detgt0 = sp.solve((alldets[0] >= 0),r1)
print('Det >= 0: ' + str(detgt0))
trlt0 = sp.solve((alltrs[0] < 0),r1)
print('Tr < 0: ' + str(trlt0))
print('\n')

print('Condições de estabilidade para o 2º ponto de equilíbrio: ')
detgt0 = sp.solve((alldets[1] >= 0),w21) #sp.simplify no effect here
print('Det >= 0: ' + str(detgt0))
trlt0 = sp.solve((alltrs[1] < 0),w21)
print('Tr < 0: ' + str(trlt0))
print('\n')

print('Condições de estabilidade para o 3º ponto de equilíbrio: ')
#detgt0 = reduce_rational_inequalities((alldets[2] >= 0),w12) #sp.simplify no effect here
detgt0 = sp.solve((alldets[2] >= 0),w12) #sp.simplify no effect here
print('Det >= 0: ' + str(detgt0))
trlt0 = sp.solve((alltrs[2] < 0),w12)
print('Tr < 0: ' + str(trlt0))
print('\n')

# Not working 
# print('Condições de estabilidade para o 4º ponto de equilíbrio: ')
# #detgt0 = reduce_rational_inequalities((alldets[2] >= 0),w12) #sp.simplify no effect here
# detgt0 = sp.solve((alldets[3] >= 0),w12,w21) #sp.simplify no effect here
# print('Det >= 0: ' + str(detgt0))
# trlt0 = sp.solve((alltrs[3] < 0),symbols)
# print('Tr < 0: ' + str(trlt0))
# print('\n')
