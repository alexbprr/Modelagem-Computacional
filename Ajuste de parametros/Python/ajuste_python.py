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
    N0 = x[0]
    r = x[1]
    k = x[2]
    params = (r, k)
    P = [N0]
    
    def solveOde(t, y):
        return odeSystem(t, y, *params)

    results = solve_ivp(solveOde,(0, tfinal), P, t_eval=t, method='Radau')    

    N = results.y[0,:]
    error = 0
    sumobs = 0
    i = 0
    j = 0
    for t in np.arange(0,tfinal,dt):
        if isReferenceTime(reference_times,t):
            error += (N[j] - data[i])*(N[j] - data[i]) 
            sumobs += data[i]*data[i]
            i = i + 1
        j = j + 1

    error = math.sqrt(error/sumobs) #Erro norma 2                  
    return error 

if __name__ == "__main__":
    global reference_times, data 
    reference_times = np.loadtxt(path+'t.dat', delimiter='\n')
    data = np.loadtxt(path+'data.dat', delimiter='\n')
    bounds = [
        (1, 200), (0.01, 1), (1,200)
    ]
    #chama evolução diferencial, result contém o melhor individuo
    solucao = differential_evolution(solve, bounds, tol=10**(-4), 
    strategy='best1bin', maxiter=100,popsize=50, disp=True, workers=4)
    print(solucao.x)
    print(solucao.success)
    #saving the best offspring...
    np.savetxt('solucao_ajuste_1.txt',solucao.x, fmt='%.2f')        
    best=solucao.x
    error = solve(best)
    print("ERROR ", error)