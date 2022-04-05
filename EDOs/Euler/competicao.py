#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 20:39:45 2016

@author: alexandre
"""

import matplotlib.pyplot as plt
import numpy as np

#parâmetros N1
r1 = 0.5
#parâmetros N2
r2 = 0.5
#coeficientes de competição
w11 = 1.
w21 = 1.5
w12 = 1.2
w22 = 1.

x = 0.2
y = 0.2
xresult = [x]
yresult = [y]

#deltaT=0.0008
#days=20

t_final = 10000
t_ini = 0
num_iteracoes = 1000000 #(t_final - t_ini)/deltaT
t = np.linspace(t_ini, t_final, num_iteracoes)
dt = (t_final - t_ini)/num_iteracoes

for t in range(t_final):
    nextx = (r1*x*(1 - w11*x - w21*y))*dt + x
    nexty = (r2*y*(1 - w22*y - w12*x))*dt + y
    xresult.append(nextx)
    yresult.append(nexty)
    x,y = nextx,nexty

fig = plt.figure(figsize=(9,6))
plt.plot(xresult, 'b-', label='N1')
plt.plot(yresult, 'g--', label='N2')
plt.xlabel('time')
plt.ylabel('amount')
plt.legend()
plt.show()

fig = plt.figure(figsize=(9,6))
plt.plot(xresult,yresult)
plt.xlabel('N1')
plt.ylabel('N2')
plt.legend()
plt.show()
