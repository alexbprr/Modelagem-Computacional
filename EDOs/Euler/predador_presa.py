#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 20:39:45 2016

@author: alexandre
"""

import matplotlib.pyplot as plt
import numpy as np

#parâmetros presa
r = 0.3
#k1 = 50
a = 0.005

#parâmetros predador
m = 0.3
c = 0.8
#k2 = 10

t_final = 100
t_ini = 0
num_iteracoes = 500
t = np.linspace(t_ini, t_final, num_iteracoes)
dt = (t_final - t_ini)/num_iteracoes

h = [0. for i in range(num_iteracoes)]
p = [0. for i in range(num_iteracoes)]
#condição inicial 
h[0] = 30
p[0] = 10

for indice in range(0,num_iteracoes-1):
    h[indice+1] = (r*h[indice] - a*h[indice]*p[indice])*dt + h[indice]
    p[indice+1] = (a*h[indice]*p[indice] - m*p[indice])*dt + p[indice]

fig = plt.figure(figsize=(9,6))
plt.plot(t, h, 'b-', label='presa' )
plt.plot(t, p, 'g--', label='predador')
plt.legend()
plt.show()
fig.savefig('predador_presa_classico.pdf')
