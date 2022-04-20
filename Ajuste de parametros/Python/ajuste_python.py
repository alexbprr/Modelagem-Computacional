#!/usr/bin/python
# -*- coding: utf-8 -*-

from xml.parsers.expat import model
import pandas as pd
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp
from math import factorial as fat
from scipy.optimize import fmin,differential_evolution
import math

path = 'data/'

def odeSystem(t, P, r, k):
    N = P[0]
    dN_dt = r*N*(1 - N/k)
 
    return [dN_dt]  

def setInitialCondition():
    N = 10
    P = [N]
    return P

def setTime():
    dt = 0.010000    
    tfinal = 50
    nsteps = int(tfinal/dt)    
    t=np.linspace(0,tfinal,nsteps)
    return (tfinal,t, dt)

def isReferenceTime(times, ct):
    for t in times: 
        if (abs(ct - t) <= 10**(-5)):
            return True 
    return False

def solve(x):
    (tfinal,t,dt) = setTime()
    P = setInitialCondition()
    r = x[0]
    k = x[1]
    params = (r,k)
    #P[0] = temp[0]
    #temp = np.delete(temp,0)
    def solveOde(t, y):
        return odeSystem(t, y, *params)

    results = solve_ivp(solveOde,(0, tfinal), P, t_eval=t, method='Radau')    

    N = results.y[0,:]    
    error = 0
    i = 0
    j = 0
    for t in np.arange(dt,tfinal+dt,dt):
        if isReferenceTime(reference_times,t):
            #print('t ' + str(t))
            #print('data ' + str(data[i]))
            error += (N[j] - data[i])*(N[j] - data[i])            
            i = i + 1
        j = j + 1
                 
    return error 

if __name__ == "__main__":
    global reference_times, data 
    reference_times = np.loadtxt(path+'t.dat', delimiter='\n')
    data = np.loadtxt(path+'data.dat', delimiter='\n')
    bounds = [
        (0.01, 1), (1,200)
    ]
    #chama evolução diferencial, result contém o melhor individuo
    solucao = differential_evolution(solve, bounds, strategy='best1bin', maxiter=50,popsize=50, disp=True, workers=4)
    print(solucao.x)
    print(solucao.success)
    #saving the best offspring...
    np.savetxt('solucao_ajuste_1.txt',solucao.x, fmt='%.2f')        
    best=solucao.x
   
    #saving the samples for UQ
    #np.savetxt('execution_de.txt',execution_de)
   
    #error, N = solve(best)
    error = solve(best)
    print("ERROR ", error)
    
    """ plt.figure();
    plt.title("Virus")
    plt.plot(ts, np.log10(V), label='Viremia model', linewidth=1.5)
    plt.errorbar(dadosViremia['Day']-1, virus_mean, yerr=[virus_min, virus_max],linestyle='None', label='Data', fmt='o', color='red', capsize=4, elinewidth=2)
    plt.legend(loc="best",prop={'size': 13})
    plt.grid(True)
    plt.tight_layout()    
    plt.savefig('output_survivor_virus.pdf',bbox_inches='tight',dpi = 300) """