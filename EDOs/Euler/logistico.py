#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np

def logistico(a, n, K, i):
    return a*n[i]*(1 - c*n[i]/K)

N = 1000
n = []
n.append(200) #condição inicial
a = 0.1 #taxa de crescimento
k = 100 #capacidade de suporte

c = 1 #competição intra-específica
t = np.linspace(0, 50, N)
dt = (50 - 0)/N
for i in range(0,N-1):
    n.append(logistico(a, n, k, i)*dt + n[i])

"""
Plot your data
"""
ax = plt.axes()
#ax.grid(True)
#ax.set_xticks(list(np.linspace(0, 50, )))
plt.xlim([0, 50])
plt.ylim([0,2*k+1])
plt.xlabel("tempo")
plt.ylabel("N")
plt.plot(t, n, 'b',label="Crescimento logístico")
plt.legend(loc='lower right')
plt.show()
