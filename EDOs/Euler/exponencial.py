# -*- coding: utf-8 -*-
"""
Spyder Editor

Este é um arquivo de script temporário.
"""

import matplotlib.pyplot as plt
import numpy as np

N = 11 #Passos de tempo
tini = 0
tfinal = 10
t = np.linspace(tini, tfinal, N)
print(t)
b = 0.1
#dt = 1
dt = (tfinal - tini)/N
#População: x 
x = [] #inicialização de x
x.append(1) #condição inicial

y = np.exp(b*t)

for i in range(0,N-1):
    x.append(x[i] + (b*x[i])*dt) #coloca no final da lista    
    print(x[i+1])

plt.plot(t, x, 'r',label="Crescimento exponencial - Derivada")
plt.plot(t, y, 'b',label="Crescimento exponencial - solucao analitica")

plt.legend(loc='upper left')
plt.show()
