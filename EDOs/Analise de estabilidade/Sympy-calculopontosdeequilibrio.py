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

""" zn = sp.Symbol('zn')
Na_e = sp.Symbol('Na_e')
K_e = sp.Symbol('K_e')
g_na = sp.Symbol('g_na')
g_k = sp.Symbol('g_k')
R = sp.Symbol('R')

Na_int(1) = 10;  % condicoes iniciais do modelo 
K_int(1) = 130;
alpha_int(1) = 0.869565217;
V(1) = -70;
pint(1) = -266.64;

zn   = 1
Na_e = 145
K_e  = 3.5
g_na = 1.22258163*10^(-2)
g_k  = 6.28564508*10^(-2)
R = 8.314462
F = 96400
T = 370.0
Cm = 4.16121841*10^(-3)
gama = 6.4103*10^(-3)
S = 1.0;
pe = 266.64;
alpha_int0 = 0.869565217;
ni_w = 1.7*10^(-5);
I_max = 13;
X_int = 3.669*10^(-6);
w_t = 4.314*10^(-7);
K_k  = 2;
K_na = 7.7;


Na_int = sp.Symbol('Na_int')
K_int = sp.Symbol('K_int')
alpha_int = sp.Symbol('alpha_int')
V = sp.Symbol('V')
pint = sp.Symbol('pint')
 

def dNa_int():
    I_na = g_na*(V(i) - (R*T/F)*log(Na_e/Na_int(i)))
    return -gama*(I_na + 2*I_max*(1 + (K_k/K_e)^(-2))*(1 + (K_na/Na_int(i))^(-3)))
  
def dK_Int():
    I_k  = g_k*(V(i) - (R*T/F)*log(K_e/K_int(i)));
    return -gama*(I_k - 3*I_max*(1 + (K_k/K_e)^(-2))*(1 + (K_na/Na_int(i))^(-3)))

def dalpha_int():
    return gama*ni_w*(pint(i) - pe - R*T*(X_int/(alpha_int(i)*w_t)) - Na_e - K_e)
  
def dV():
    return (-I_na - I_k)/Cm
    
pint(i) = pe + S*(alpha_int(i) - alpha_int0)
 """

# EDOs do sistema
# EDO da 1ª população
def dudt(u,v,params):
    return params[0]*u - params[1]*u*v

# EDO da 2ª população
def dvdt(u,v,params):
    return params[1]*u*v - params[2]*v

#Criar um Symbol para cada variável e parâmetro do sistema

u = sp.Symbol('u') #representa a presa no predador-presa
v = sp.Symbol('v') #representa o predador no predador-presa

r = sp.Symbol('r')
a = sp.Symbol('a')
m = sp.Symbol('m')
params = [r, a, m]

dudt_equal0 = sp.Eq(dudt(u,v,params),0)
dvdt_equal0 = sp.Eq(dvdt(u,v,params),0)

eqpoints = sp.solve((dudt_equal0, dvdt_equal0),(u,v)) #Resolve a equação para as variáveis u e v
print("Pontos de equilíbrio: " + str(eqpoints))

nullcline_dudt = sp.solve(dudt_equal0,(u,v))
nullcline_dvdt = sp.solve(dvdt_equal0,(u,v))
print("Nullcline de du/dt: " + str(nullcline_dudt))
print("Nullcline de dv/dt: " + str(nullcline_dvdt))

eqMat = sp.Matrix([ dudt(u,v,params), dvdt(u,v,params) ])
Mat = sp.Matrix([ u, v ])
matrizJacobiana = eqMat.jacobian(Mat)
print('Matriz jacobiana %s' % matrizJacobiana)

trace_jac = matrizJacobiana.trace()
print('Traço\n' + str(matrizJacobiana.trace()))

#mat = np.array(matrizJacobiana)
det_jac = matrizJacobiana.det()
print('Determinante\n' + str(matrizJacobiana.det()))

# from sympy import Poly
# >>> from sympy.abc import x
# >>> from sympy.solvers.inequalities import solve_poly_inequality
#
# >>> solve_poly_inequality(Poly(x, x, domain='ZZ'), '==')

#expr = 2*x + y
#expr2 = expr.subs(x, 2)

#eq1 = sp.Eq(det_jac, >0.)
