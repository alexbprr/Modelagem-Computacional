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

t = getArrayFromFile('t.dat')
h = getArrayFromFile(h_filename+'.dat') 
p = getArrayFromFile(p_filename+'.dat') 

#Read experimental data from file 
#h_data = getArrayFromFile('prey_data.dat') 
#p_data = getArrayFromFile('predator_data.dat') 

times = [2,15,30,45.045,49,50]
popdata = [22,
        52.8,
        83.4,
        95.759,
        97,
        97.3]

# fig = plt.figure(figsize=(10,7))
# plt.axis([0, max(t), 0, max(max(h),max(p))])
# plt.tick_params(labelsize=16)
# plt.title('Predador x Presa', fontsize=18)
# plt.xlabel('tempo (dias)',fontsize=16)
# plt.ylabel('população',fontsize=16)
# ax = fig.gca()
# ax.plot(t,h,label="Presa", color=colorVal)
# color_index += 8
# colorVal = scalarMap.to_rgba(color_index)
# ax.plot(t,p,label="Predador", color=colorVal)
# x_val = [x[0] for x in p_data]
# y_val = [x[1] for x in p_data]
# ax.plot(x_val,y_val,'o',label='p_data', color='orange')
# ax.legend(loc='upper right', fontsize=15)
# ax.grid()
# plt.show()
# generateGraphs(fig,'Predador_Presa')

fig = plt.figure(figsize=(10,7))
ax = fig.gca()
x_val = [x for x in times]
y_val = [x for x in popdata]
ax.plot(x_val,y_val,'o',label='data', color='blue')

def logistico(a, n, K, i):
    return a*n[i]*(1 - n[i]/K)

N = 1000
n = []
n.append(19.87) 
a = 0.10307
k = 98.60712 
tfinal = 50
t = np.linspace(0, tfinal, N)
dt = (tfinal - 0)/N

for i in range(0,N-1):
    n.append(logistico(a, n, k, i)*dt + n[i])

#plt.xlim([0, tfinal])
#plt.ylim([0,2*k+1])
plt.xlabel("tempo")
plt.ylabel("N")
ax.plot(t, n, 'r',label="Crescimento logístico")
ax.legend(loc='lower right')
plt.show()