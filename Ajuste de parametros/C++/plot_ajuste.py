#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
import sys
import os
from pylab import array 

def getMatrixFromFile(filename, u):
    file = open(filename, 'r')
    if (file):
        temp = [content.split(' ') for content in [content.rstrip() for content in file.readlines()]]
        cont = 0
        for j in temp:
            u[cont] = [float(i) for i in j]
            cont+=1
        return u
    else:
        print('Não foi possível abrir o arquivo!')
        exit

def getArrayFromFile(filename):
    a = []
    with open(filename) as afile:
        a = afile.readlines()
    a = [content.strip() for content in a]
    a = [float(i) for i in a]
    return a

def createFigure(title, x, maxDensity):
    fig = plt.figure(figsize=(9,6))
    plt.axis([0, x[size-1]+0.1, 0, maxDensity + 0.01])
    plt.tick_params(labelsize=18)
    plt.title(title, fontsize=19)
    plt.xlabel('x (mm)',fontsize=18)
    plt.ylabel('concentration (amount/mm³)',fontsize=18)
    ax = fig.gca()
    return [fig,ax]
#'graphs/' +
def generateGraphs(fig, filename):    
    fig.savefig(filename + '.svg', format='svg', bbox_inches='tight')
    
def initializeGraphParameters():
    global scalarMap
    cm = plt.get_cmap('jet')
    cNorm  = colors.Normalize(vmin=0, vmax=10)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

def getMaxValue(mat):
    maxvalue = 0
    for i in range(0,T):
        for j in mat[i]:
            if(j > maxvalue):
                maxvalue = j
    return maxvalue

initializeGraphParameters()
color_index = 0
colorVal = scalarMap.to_rgba(color_index)

h_filename = 'presa'
p_filename = 'predador'

path = 'data/'
t = getArrayFromFile(path+'t.dat')

times = [2,15,30,45.045,49,50]
popdata = [22,
        52.8,
        83.4,
        95.759,
        97,
        97.3]

fig = plt.figure(figsize=(10,7))
ax = fig.gca()
x_val = [x for x in times]
y_val = [x for x in popdata]
ax.plot(x_val,y_val,'o',label='data', color='blue')

def logistico(a, n, K, i):
    return a*n[i]*(1 - n[i]/K)

T = 1000
n = []
n.append(19.87) 
r = 0.10307
k = 98.60712 
tfinal = 50
t = np.linspace(0, tfinal, T)
dt = (tfinal - 0)/T

for i in range(0,T-1):
    n.append(logistico(r, n, k, i)*dt + n[i])

#plt.xlim([0, tfinal])
#plt.ylim([0,2*k+1])
plt.xlabel("tempo")
plt.ylabel("N")
ax.plot(t, n, 'r',label="Crescimento logístico")
ax.legend(loc='lower right')
plt.show()