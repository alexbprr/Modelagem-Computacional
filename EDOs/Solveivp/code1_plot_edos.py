# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
import pandas as pd

def readFile(filename):
    file = open(filename, 'r')
    if (file):        
        return pd.read_csv(filename)
    else: 
        print('Error: it was not possible to open the file!')
        exit

def getArrayFromFile(filename):
    file = open(filename, 'r')
    if (file):        
        return np.loadtxt(file, delimiter='\n')
    else: 
        print('Error: it was not possible to open the file!')
        exit

def createFigure(title):
    fig = plt.figure(figsize=(10,7))    
    plt.tick_params(labelsize=18)
    plt.title(title, fontsize=19)
    plt.xlabel('time',fontsize=18)
    plt.ylabel('concentration',fontsize=18)
    ax = fig.gca()
    return [fig,ax]
    
def initializeGraphParameters():
    global scalarMap, color_index, colorVal
    cm = plt.get_cmap('jet')
    cNorm  = colors.Normalize(vmin=0, vmax=50)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    color_index = 0
    colorVal = scalarMap.to_rgba(color_index)

def saveFig(fig, filename):
    fig.savefig(dir_ + filename + '.svg', format='svg', bbox_inches='tight')

dir_ = ''
initializeGraphParameters()

odeValues = readFile(dir_+'oderesults.csv')
print(odeValues)
t = odeValues['Time']
E = odeValues['E']
T = odeValues['T']
I = odeValues['I']
V = odeValues['V']

fig,ax = createFigure('E')
color_index += 50/4
ax.plot(t, E, label="E", color=colorVal)
colorVal = scalarMap.to_rgba(color_index)
ax.legend(loc='upper right', fontsize=15)
ax.grid()
saveFig(fig,dir_+'E')
fig.clear()
ax.clear()

fig,ax = createFigure('V')
color_index += 50/4
colorVal = scalarMap.to_rgba(color_index)
ax.plot(t,V,label="V", color=colorVal)
ax.legend(loc='upper right', fontsize=15)
ax.grid()
saveFig(fig,dir_+'V')
fig.clear()
ax.clear()

fig,ax = createFigure('I')
color_index += 50/4
colorVal = scalarMap.to_rgba(color_index)
ax.plot(t,I,label="I", color=colorVal)
ax.legend(loc='upper right', fontsize=15)
ax.grid()
saveFig(fig,dir_+'I')
fig.clear()
ax.clear()

fig,ax = createFigure('T')
color_index += 50/4
colorVal = scalarMap.to_rgba(color_index)
ax.plot(t,T,label="T", color=colorVal)
ax.legend(loc='upper right', fontsize=15)
ax.grid()
saveFig(fig,dir_+'T')
fig.clear()
ax.clear()