# -*- coding: utf-8 -*-
from scipy.integrate import solve_ivp
import numpy as np
import pandas as pd
# Função que define o sistema de EDOs 
def odeSystem(t, P, b, k, s, m, a, k2, p):
    E, V, I, T = P[0], P[1], P[2], P[3]

    dE_dt = -b*E*V + 0.001*(10000 - E)  #EDO de E 
    dV_dt = +p*a*I -k2*V  #EDO de V  
    dI_dt = +b*E*V -k*I*T -a*I  #EDO de I 
    dT_dt = +s*I -m*T  #EDO de T 
 
    return [dE_dt,dV_dt,dI_dt,dT_dt]  

def solveOde(t, y):
    global model_args
    return odeSystem(t, y, *model_args)

def setInitialCondition():
    E = 10000
    V = 5
    I = 0
    T = 10
    P = [E, V, I, T]
    return P

def setParameters():
    b = 0.002
    k = 0.3
    s = 0.2
    m = 0.05
    a = 0.05
    k2 = 0.4
    p = 5
    model_args = (b, k, s, m, a, k2, p)
    return model_args

def setTime():
    dt = 0.0010000    
    tfinal = 500
    nsteps = int(tfinal/dt)
    print('nsteps ' + str(nsteps))
    t=np.linspace(0,tfinal,nsteps)
    return (tfinal,t)

if __name__ == "__main__":
    global model_args
    (tfinal,t) = setTime()
    P = setInitialCondition()
    model_args = setParameters()     
    results = solve_ivp(solveOde,(0, tfinal), P, t_eval=t, method='Radau')
    dir_ = ''
    df = pd.DataFrame(results.y.transpose(), columns = ['E', 'V', 'I', 'T'])
    df.insert(0, 'Time', results.t)
    df.to_csv(dir_+'oderesults.csv', float_format='%.5f', sep=',')
