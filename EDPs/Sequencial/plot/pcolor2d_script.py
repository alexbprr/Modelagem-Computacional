import matplotlib.pyplot as plt 
import matplotlib.colors as colors
import matplotlib.cm as cmx
import pandas as pd
import numpy as np 

def plot(fname):
    x = np.arange(0,5.1,0.1)
    y = np.arange(0,5.1,0.1)
    m = np.zeros((len(x),len(y)))
    df = pd.read_csv(fname)
    cm = plt.get_cmap('viridis')
    value = df.value 
    cNorm  = colors.Normalize(vmin=min(value), vmax=max(value))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    k = 0
    for j in np.arange(0,len(y)):
        for i in np.arange(0,len(x)):    
            m[i][j] = value[k]
            k += 1

    fig = plt.figure(figsize=(10,7))
    ax = fig.gca()
    plt.title(fname, fontsize=17)
    plt.xlabel('y',fontsize=17)
    plt.ylabel('x',fontsize=17)
    plt.axis([min(x), max(x), min(y), max(y)])
    pcm = plt.pcolormesh(x, y, m, cmap=cm)
    plt.colorbar(pcm)
    fig.savefig(fname+'.svg', format='svg', bbox_inches='tight')
    plt.close()

plot("A2d_t0.csv")
plot("A2d_t1.csv")
plot("A2d_t5.csv")
plot("A2d_t10.csv")
plot("A2d_t20.csv")
plot("A2d_t30.csv")
plot("A2d_t40.csv")
plot("A2d_t50.csv")

plot("C2d_t0.csv")
plot("C2d_t1.csv")
plot("C2d_t5.csv")
plot("C2d_t10.csv")
plot("C2d_t20.csv")
plot("C2d_t30.csv")
plot("C2d_t40.csv")
plot("C2d_t50.csv")

plot("P2d_t0.csv")
plot("P2d_t1.csv")
plot("P2d_t5.csv")
plot("P2d_t10.csv")
plot("P2d_t20.csv")
plot("P2d_t30.csv")
plot("P2d_t40.csv")
plot("P2d_t50.csv")

plot("N2d_t0.csv")
plot("N2d_t1.csv")
plot("N2d_t5.csv")
plot("N2d_t10.csv")
plot("N2d_t20.csv")
plot("N2d_t30.csv")
plot("N2d_t40.csv")
plot("N2d_t50.csv")