# -*- coding: utf-8 -*-
from scipy.integrate import solve_ivp
import numpy as np
import pandas as pd

# Função que define o sistema de EDOs 
def odeSystem(t, u, r, a, b, m):
    H, P = u[0], u[1]

    dH_dt = r*H - a*H*P
    dP_dt = b*H*P - m*P
 
    return [dH_dt,dP_dt]  

def solveOde(t, y):
    global model_args
    return odeSystem(t, y, *model_args)

def setInitialCondition():
    H = 100
    P = 5
    u = [H, P]
    return u

def setParameters():
    r = 0.1
    a = 0.03
    b = 0.05
    m = 0.1
    model_args = (r, a, b, m)
    return model_args

def setTime():
    dt = 0.010000    
    tfinal = 100
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
    df = pd.DataFrame(results.y.transpose(), columns = ['Prey', 'Predator'])
    df.insert(0, 'Time', results.t)
    df.to_csv(dir_+'oderesults.csv', float_format='%.5f', sep=',')
