#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 14:24:56 2016

@author: alexandre
"""
import sys
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from sympy.utilities.lambdify import lambdify
import scipy.optimize as opt
from sympy import diff
from sympy import sqrt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn
import decimal as dc
dc.getcontext().prec = 2
import scipy.integrate as scint
#from sympy.solvers.solveset import nonlinsolve
from scipy.optimize import fsolve

#Graph setup
seaborn.set(style='ticks')
fig = plt.figure(figsize=(10,7))
ax = fig.gca()
ax.grid(True, which='both')
seaborn.despine(ax=ax, offset=0)

#Parameters initialization
fileEndingName = str(sys.argv[1])
N_POINTS = 50
x = np.linspace(0.0, 1.0, N_POINTS)
y = np.linspace(0.0, 1.0, N_POINTS)
x_y = [[[x[p],y[k]] for k in range(0,N_POINTS)] for p in range(0,N_POINTS)]

u_key = str(sys.argv[2])
v_key = str(sys.argv[3])
f_key = 'f'
u = sp.Symbol(u_key)
v = sp.Symbol(v_key)
f = sp.Symbol(f_key)
variables = [u,v]
#bacteria parameters
r = 1.
m_b = 0
a = 1
#neutrophil parameters
s = 1
l = 5

#fibrin parameters
k = 3
p = 3

#weights
wbb=1.
wcoab=1.1 #2*0
wnb=1.1
wfibb=3. #1.2

wbcoa=1.1
wcoacoa=1.
wncoa=1.2
wfibcoa=0.

#wbfib=1.2
wcoafib=0
#wnfib=0
#wfibfib=1

wbn=1.1
wcoan=0.5
wnn=1
wfibn=1

wbto=0
wcoato=0
wnto=0
wfibto=0
wtoto=1

lambda_ = 4
halfsat = 0.2

def inhibition(alpha, gamma, population):
    return (alpha + 1)*pow(population, gamma)/(alpha + pow(population, gamma))

# define system in terms of separated differential equations
#def fib(f,b):
    #return r1*f*(1 - f/k1 - b12*b/k1) -m*f
 #   return k*b*f*(1 - a1*b - a2*f)
def dudt(u,v):
    return k*u*(1 - wbcoa*inhibition(0.15,2,v) - wcoacoa*u)
    #return s*u*v/(1 + wbn*u + wnn*v) - a*u*v
    #return p*v*(1 - v)*(1 - wbcoa*v)/(wcoacoa)
    #return k*(1 - v)*(1 - wbcoa*v)/(wcoacoa)
    #return k*u*(1 - wbcoa*v - wcoacoa*u - wfibcoa*u*v)#*(1 - u - v))
    #return k*u*v*(1 - wbcoa*v - wcoacoa*u - wfibcoa*p*v*(1 - v)*u)
    #return k*u*v*(1 - wbcoa*v - wcoacoa*u - wfibcoa*p*(1 - v)*u)
    #return r*v*(1 - wbb*v - wcoab*1 - wfibb*u) + k*v*(1 - wbcoa*v - wcoacoa*1 - wfibcoa*u)
    #return p*v*(1 - wbcoa*v - wfibcoa*u - p/k)/wcoacoa
    #return k*u*v*(1 - wbcoa*v - wcoacoa*u - wfibcoa*u*v)
    #return k*((1 + halfsat)*pow(v,lambda_))/(halfsat + pow(v,lambda_))*(1 - wbcoa*v - wcoacoa*u - wfibcoa*u*v)
    #return p*v*(1 - wbcoa*v)*(1 - wfibfib*u - wbfib*v) considerando muito rapida a dinamica de coa: dcoa/dt = 0

def dvdt(u,v):
    return r*v*(1 - wcoab*inhibition(0.15,2,u) - wbb*v)
    #return r*v*(1 - wbb*v -wcoab*(1 - wbcoa*v)/(wcoacoa) - wfibb*u)
    #return r*v*(1 - wbb*v -wcoab*u - wfibb*p*(1 - v)*u)
    #return r*v*(1 - wbb*v - wcoab*(1 - wbcoa*v)/(wcoacoa) - wfibb*u)
    #return r*v*(1 - wbb*v -wcoab*u - wfibb*u*v)#*(1 - u - v))
    #return r*v*(1 - wbb*v -wcoab*u - wfibb*p*v*(1 - v)*u) #Ultim0
    #return r*v*(1 - wbb*v - wcoab*u -wfibb*u*v)
    #return r*v*(1 - wbb*v - wcoab*1 - wfibb*u)
    #return r*v*(1 - wbb*v - wcoab*u - wfibb*u*v)
    #return r*v*(1 - wcoab*(1 - wbcoa*v) - wfibb*u - wbb*v) considerando muito rapida a dinamica de coa: dcoa/dt = 0

def f(x):
    return np.array([k*x[0]*x[1]*(1 - wbcoa*x[1] - wcoacoa*x[0] - wfibcoa*x[0]*x[1]),
    r*x[1]*(1 - wbb*x[1] - wcoab*x[0] - wfibb*x[0]*x[1])])

def equations(p):
   u,v = p
   eq1 = dudt(u,v)
   eq2 = dvdt(u,v)
   return (eq1, eq2)

    #r*b*(1 - b - f) - (1 - (1/(l*f + 1)))*b
#def neu(n,b):
#    return k*(pow(b,lambda_))/(phi + pow(b,lambda_))*(1 - a1*b - a2*n)
    #return (s*n*b)*(1 - d1*b - d3*n)
    #return (s*n*b + s*b)*(1 - d1*b - d3*n)
    #return (s*n*b + s)*(1 - d1*b - d3*n)
U = dudt(u,v)
V = dvdt(u,v)

eqMat = sp.Matrix([ U, V ])
Mat = sp.Matrix([ u, v ])
jacMat = eqMat.jacobian(Mat)
print('Jacobian Matrix\n %s' % jacMat)
print('---------------------')

UEqual = sp.Eq(U, 0.)
VEqual = sp.Eq(V, 0.)
fp = []
eqp = []

print("Equilibrium points: " + str(eqp))
# for e in eqp:
#     if (u in e and v in e):
#         #print('value: ' + str(e.get(u)))
#         fp.append((e[u],e[v]))
#
# print("Equilibrium points tuple list: " + str(fp))

if not fp:
    fp = eqp

print("fp var: " + str(fp))

fp_file = open('eqpoints_' + fileEndingName + '.txt','w')
nullclinesFile = open('nullclines' + fileEndingName + '.txt','w')

intersects_eixou = []
intersects_eixov = []
"""
v-nullclines
"""
funcoes_lambda_v = []
res_v = sp.solve(dvdt(u,v))
nullclinesFile.write('v-nullclines: %s' %str(res_v))
nullclinesFile.write('\n--------------------------\n')
for i in range(len(res_v)):
    for key in res_v[i]:
        #funcoes_lambda_b.append(lambdify(key, sp.sympify(res_b[i][key]),'numpy'))
        if str(key) == v_key:
            funcoes_lambda_v.append(lambdify(u, sp.sympify(res_v[i][v]),'numpy'))
        elif str(key) == u_key:
            funcoes_lambda_v.append(lambdify(v, sp.sympify(res_v[i][u]),'numpy'))

nullclines_v =  [[(k,k) for k in range(0,N_POINTS)] for p in range(0,len(res_v))]
#calcula as interseções com os eixos primeiro para ajudar a plotar o gráfico
for i in range(0,len(funcoes_lambda_v)):
    for j in range(0,N_POINTS):
        for key in res_v[i]:
            if str(key) == v_key:
                nullclines_v[i][j] = (x[j],funcoes_lambda_v[i](x[j]))
                #intersects_eixov.append(funcoes_lambda_v[i](0))
            elif str(key) == u_key:
                nullclines_v[i][j] = (funcoes_lambda_v[i](x[j]),x[j])
                #intersects_eixou.append(funcoes_lambda_v[i](0))

"""
u-nullclines
"""
funcoes_lambda_u = []
res_u = sp.solve(dudt(u,v))
nullclinesFile.write('u-nullclines: %s' %str(res_u))
nullclinesFile.write('\n--------------------------\n')
for i in range(len(res_u)):
    for key in res_u[i]:
        if str(key) == v_key:
            funcoes_lambda_u.append(lambdify(u, sp.sympify(res_u[i][v]),'numpy'))
        elif str(key) == u_key:
            funcoes_lambda_u.append(lambdify(v, sp.sympify(res_u[i][u]),'numpy'))

nullclines_u =  [[(k,k) for k in range(0,N_POINTS)] for p in range(0,len(res_u))]

for i in range(0,len(funcoes_lambda_u)):
    for j in range(0,N_POINTS):
        for key in res_u[i]:
            if str(key) == v_key:
                nullclines_u[i][j] = (x[j],funcoes_lambda_u[i](x[j]))
                #intersects_eixov.append(funcoes_lambda_u[i](0))
            elif str(key) == u_key:
                nullclines_u[i][j] = (funcoes_lambda_u[i](x[j]),x[j])
                #intersects_eixou.append(funcoes_lambda_u[i](0))

def find_fixed_points():
    fp = []
    for i in range(0,N_POINTS):
        for j in range(0,N_POINTS):
            if (dvdt(x_y[i][j][0],x_y[i][j][1]) == 0) and (dudt(x_y[i][j][0],x_y[i][j][1]) == 0) :
                fp.append((x_y[i][j][0],x_y[i][j][1]))
                #print('The system has a fixed point in %s,%s' % (x_y[i][j][0],x_y[i][j][1]))
            for m in range(0,len(nullclines_u)):
                for n in range(0,len(nullclines_v)):
                    if (nullclines_u[m][i][0] == nullclines_v[n][j][0] and nullclines_u[m][i][1] == nullclines_v[n][j][1]):
                        fp.append(nullclines_u[m][i])
                        #print("Intersection points: " + str(nullclines_u[m][i]))
    return fp

fp.extend(find_fixed_points())

def eigenvalues(uvalue, vvalue):
    lambda_dF_f = lambdify(variables,diff(U,u))
    a11 = lambda_dF_f(uvalue,vvalue)

    lambda_dF_b = lambdify(variables,diff(U,v))
    a12 = lambda_dF_b(uvalue,vvalue)

    lambda_dB_f = lambdify(variables,diff(V,u))
    a21 = lambda_dB_f(uvalue,vvalue)

    lambda_dB_b = lambdify(variables,diff(V,v))
    a22 = lambda_dB_b(uvalue,vvalue)

    tr = a11 + a22
    det = a11*a22 - a12*a21
    #print('Trace = ' + str(tr))
    #print('Determinant = ' + str(det))
    tr2_4d = tr**2 - 4*det;
    #print('Tr² -4Det = ' + str(tr2_4d))
    sqrt_tr2_4d = sqrt(tr2_4d)
    #lambda1 = (tr - sqrt_tr2_4d)/2 #(tr - sqrt(tr**2 - 4*det))/2
    #lambda2 = (tr + sqrt_tr2_4d)/2
    #print('Check the fixed point  %s, %s' % (x,y))
    if isinstance(sqrt_tr2_4d,complex):
        #print('The real part of the first eigenvalue is %s' %lambda1.real)
        #print('The real part of the second eigenvalue is %s' %lambda2.real)
        if(tr < 0):
            #print('The fixed point in %s, %s is a spiral sink. It is stable.' % (x,y))
            #u_str = "{:.2f}".format()
            fp_file.write('The fixed point in %s, %s is a spiral sink. It is stable.' % (uvalue,vvalue) + "\n")
        elif(tr > 0):
            #print('The fixed point in %s, %s is a spiral source. It is unstable.' % (x,y))
            fp_file.write('The fixed point in %s, %s is a spiral source. It is unstable.' % (uvalue,vvalue) + "\n")
        elif(tr == 0):
            #print('The fixed point in %s, %s is a center.' % (x,y))
            fp_file.write('The fixed point in %s, %s is a center.' % (uvalue,vvalue) + "\n")
    else:
        if(det < 0):
            #if(bvalue > 0.):
            #print('The fixed point in %s, %s is a saddle.' % (x,y))
            fp_file.write('The fixed point in %s, %s is a saddle.' % (uvalue,vvalue) + "\n")
        elif (tr < 0):
            #if(bvalue > 0.):
            #print('The fixed point in %s, %s is a sink. It is stable.' % (x,y))
            fp_file.write('The fixed point in %s, %s is a sink. It is stable.' % (uvalue,vvalue) + "\n")
        elif(tr > 0):
            #if(bvalue > 0.):
            #print('The fixed point in %s, %s is a source. It is unstable.' % (x,y))
            fp_file.write('The fixed point in %s, %s is a source. It is unstable.' % (uvalue,vvalue) + "\n")
        elif(tr == 0):
            #print('The fixed point in %s, %s is a center.' % (x,y))
            fp_file.write('The fixed point in %s, %s is a center.' % (uvalue,vvalue) + "\n")
    #print('----------------------------')

#Determine the equilibrium points stability
for x,y in fp:
    eigenvalues(x,y)

max_fp_f = 0.
max_fp_b = 0.
for i in range(len(fp)):
    if(fp[i][0] == 0. and fp[i][1] > max_fp_b):
        max_fp_b = fp[i][1]
    if(fp[i][1] == 0. and fp[i][0] > max_fp_f):
        max_fp_f = fp[i][0]

#print(str(max_fp_b))
#intersection_f = max(intersects_eixof)
#intersection_b = max(intersects_eixob)
#maximo_f = max(intersection_f,float(max_fp_f))
#maximo_b = max(intersection_b,float(max_fp_b))
#print('Interects eixof:'+str(intersects_eixof))
#print('Interects eixob:'+str(intersects_eixob))
#if(maximo_f > 1.):
#    maximo_f = 1.
#if(maximo_b > 1.):
#    maximo_b = 1.

#Create, plot and save the graph
x = np.linspace(0.0, 1., N_POINTS-20)
y = np.linspace(0.0, 1., N_POINTS-20)
x1 , y1  = np.meshgrid(x,y)                    # create a grid
dx = dudt(x1,y1)
dy = dvdt(x1,y1)

#print(str(dx))
i = 0
newdx = []

hipot = np.hypot(dx, dy)
M = (hipot)                        # norm growth rate
M[ M == 0] = 1.                                 # avoid zero division errors
dx /= M                                        # normalize each arrows
dy /= M
plt.quiver(x1, y1, dx, dy, color='silver')
params = {'legend.fontsize': 'large',
          'figure.figsize': (9, 7),
         'axes.labelsize': '17',
         'axes.titlesize':'17',
         'xtick.labelsize':'17',
         'ytick.labelsize':'17'}
plt.rcParams.update(params)
plt.axis([-0.01, 1.+0.025, -0.01, 1.+0.025])
plt.xlabel(str(sys.argv[2]),fontsize=17);
plt.ylabel(str(sys.argv[3]),fontsize=17);

#Plot nullclines
for i in range(0,len(nullclines_v)):
    if i == 0:
        plt.plot(*zip(*nullclines_v[i]), 'r', label='d'+ v_key + '/dt nullcline', zorder=1);
    else:
        plt.plot(*zip(*nullclines_v[i]), 'r', zorder=1);

for i in range(0,len(nullclines_u)):
    if i == 0:
        plt.plot(*zip(*nullclines_u[i]), 'g--', label='d' + u_key + '/dt nullcline', zorder=1);
    else:
        plt.plot(*zip(*nullclines_u[i]), 'g--', zorder=1);

#fp_set = {x for x in fp}

#Plot equilibrium points
plt.scatter(*zip(*fp), marker='o', color='b', s=22, label='Equilibrium points', zorder=3)

nullclinesFile.close()
fp_file.close()

#calcula os valores de B e F para alguns pontos em torno do ponto de equilíbrio
#para os valores de B e F obtidos determine todos os pontos que satisfazem as equações do sistema

x = 0
y = 0
# deltax = 0.05;
# deltay = 0.05;
# #
u_ini = 0.1
v_ini = 0.3
initials = []
initials.append((u_ini,v_ini))
next_u = 0
next_v = 0
points = []
#
# for eq in fp:
#     ite = 0
#     while (next_u != eq[0] and next_v != eq[1] or ite < 20):
#         next_u = dudt(u_ini+x,v_ini+y)
#         next_v = dvdt(u_ini+x,v_ini+y)
#         x+= deltax
#         y+=deltay
#         points.append((next_u,next_v))
#         ite+= 1
# plt.plot(points,color='k')

# for eq in fp:
#    fvalue = dudt(eq[0]+deltax,eq[1]+deltay)
#    bvalue = dvdt(eq[0]+deltax,eq[1]+deltay)
#    eq_set = sp.solve( (sp.Eq(F, fvalue), sp.Eq(B, bvalue)), f,b)
#    print(str(eq_set))
   #plt.plot(eq_set, color='orange',marker='_')

# for y0 in initials:
#     traj_f = np.empty((0,2))
#     traj_b = np.empty((0,2))
#     t = np.linspace(0,10,100)
#     traj = scint.odeint(f,y0,t)
#     #traj = np.vstack((np.flipud(traj_b),traj_f))
#     plt.plot(traj[:,0],traj[:,1],linewidth=1.2)

#funcoes_f_b1 = sp.sympify(res_b[0][b]);#uma nullcline de b (b: ...)
#funcoes_f_b2 = sp.sympify(res_b[1][b]);#uma nullcline de b (f: ...)

#fig.set_title("Phase plane with nullclines");

leg = plt.legend(loc='upper right',frameon = 1)
frame = leg.get_frame()
frame.set_color('white')
#frame.set_linewidth(0)
plt.show()
pp_filename = 'pp_' + fileEndingName
fig.savefig(pp_filename + '.svg', format='svg', bbox_inches='tight')
