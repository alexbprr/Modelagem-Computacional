import pandas as pd 
import numpy as np 
import plotly.express as px
import argparse
from plotly.subplots import make_subplots
import plotly.graph_objects as go

def split_strip(text, sep):
    """ Separar text por sep e remover espaços em volta """
    return [s.strip() for s in text.split(sep)]

def split_relational_operator(text):
    """ Separar text por um operador relacional
        O resultado é da forma (operador, lado esquedo, lado direito)
    """

    for op in ['==', '!=', '<=', '>=', '<', '>']:
        if op in text:
            return (op, *split_strip(text, op))

    return None

def split_operator(text):
    """ Separar text por um operador (exceto os relacionais)
        O resultado é da forma (operador, lado esquedo, lado direito)
        O resultado será ('=', text, None) se nenhum for encontrado
    """
    for op in ['+', '-', '*', '/', '%', '&', '|', '^']:
        if op in text:
            return (op, *split_strip(text, op))

    return ('=', text.strip(), None)

def readPlotConfigFile(filename):    
    try: 
        fin = open(filename, 'r')
        if (fin):
            tokens = []
            
            for line in fin:
                ibuf = line.rstrip('\n')
                i = 0
                print(line)
                pass 
    except: 
        print('Error in file opening.')

# t = 0:
# 0,0,0
# 0,1,0

def readCsv(fname):
    data = pd.read_csv(fname)
    values = data.values
    size = values.size 
    return (data,values,size) #Return dataframe, values and total size 

if __name__ == '__main__':
    global t,X,Y,xsize,ysize 
    parser = argparse.ArgumentParser(description='Python script for plotting')
    parser.add_argument('-c', help="Entry plot config file")
    args = parser.parse_args()
    #readPlotConfigFile(args.c)

    (tdata,t,tsize) = readCsv("t.csv") 
    (xdata,x,xsize) = readCsv("x.csv") 
    (ydata,y,ysize) = readCsv("y.csv") 
    (zdata,z,zsize) = readCsv("z.csv")
    
    ndim = 3
    if (ysize <= 1 and zsize <= 1):
        ndim = 1
    elif (zsize <= 1):
        ndim = 2            
    time_spacedf = pd.concat([tdata,xdata,ydata], axis=1)

