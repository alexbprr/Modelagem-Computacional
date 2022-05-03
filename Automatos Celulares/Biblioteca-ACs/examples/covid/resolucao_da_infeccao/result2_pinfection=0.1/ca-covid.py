try:
    from ca import *

except ModuleNotFoundError:
    from sys import path
    path.insert(0, '..')
    from ca import *

from time import sleep
from random import randint, random, choice
from sys import argv
import random as rd

VOID = 0
VIRUS = 1
APC = 2
APM = 3
I = 4
THN = 5
THE = 6
TKN = 7
TKE = 8
C = 9

class CA_COVID(CA):
    def rule(self, x, y):
        s = self[x][y]
        n = neighbors8(self, x, y, old=True)        

        # probabilidades de "criar" novas células 
        p_apc = 0.3
        p_thn = 0.2
        p_tkn = 0.2 

        p_apm_death = 0.05
        p_c_decay = 0.1
        p_the_death = 0.1 
        p_tke_death = 0.1             

        p_infection = 0.1
        p_replicacao = 0.25

        virus = len([i for i in neighbors8(self, x, y) if i == VIRUS])
        apm = len([i for i in neighbors8(self, x, y) if i == APM])
        c = len([i for i in neighbors8(self, x, y) if i == C])
        i = len([i for i in neighbors8(self, x, y) if i == I])
        the = len([i for i in neighbors8(self, x, y) if i == THE])
        tke = len([i for i in neighbors8(self, x, y) if i == TKE])

        if s == self.getOld(x, y):

            if s == VOID: 
                if c >= 1:
                    if random() < p_apc: 
                        return APC 
                    elif random() < p_apc + p_thn:  
                        return THN 
                    elif random() <= 1: 
                        return TKN  
                if virus >= 1: 
                    if random() < 0.3: 
                        return VIRUS
                
                if i >= 1 or (apm >= 1 or the >= 1 or tke >= 1 and random() < (1 - p_replicacao)):
                    return C 
                
                if the >= 1 and random() < p_replicacao:
                    return THE                 
                if tke >= 1 and random() < p_replicacao:
                    return TKE                 
                         
            if s == VIRUS:
                if tke >= 1 or apm >= 1: 
                    return VOID 

            if s == APC:
                if virus >= 1 and random() < 0.75: 
                    return APM 
                #else: 
                #    return I  

            if s == APM: 
                if random() < p_apm_death:
                    return VOID 
                if virus >= 1 and random() < p_infection: 
                    return I 
                
                n = [i for i in neighbors8(self, x, y, pos=True)
                     if self[i] == VIRUS]
                if len(n) > 0:
                    p = n[randint(0, len(n)-1)]
                    self.move(x, y, *p) # move em direção ao vírus 

            if s == I:
                if tke >= 1:
                    return VOID 

            if s == C:
                if random() < p_c_decay:
                    return VOID 

            if s == THN:
                if apm >= 1:
                    return THE 

            if s == THE: 
                if random() < p_the_death:
                    return VOID 
                
            if s == TKN:
                if apm >= 1 and random() < 0.75: # probabilidade de ativação aumenta com C 
                    return TKE
                #else: 
                #    return I  
            
            if s == TKE: 
                if random() < p_tke_death:
                    return VOID 
                if virus >= 1 and random() < p_infection:
                    return I 
                
                n = [i for i in neighbors8(self, x, y, pos=True)
                     if self[i] == VIRUS]
                if len(n) > 0:
                    p = n[randint(0, len(n)-1)]
                    self.move(x, y, *p) # move em direção ao vírus 

        return s

c = CA_COVID(100, values=[0]*50 + [1]*1 + [2]*2 + [5]*2 + [7]*1)

#number_of_colors = 10
#color = ["#"+''.join([choice('0123456789ABCDEF') for j in range(6)])            for i in range(number_of_colors)]
#def get_cmap(n, name='hsv'):
#    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
#    RGB color; the keyword argument name must be a standard mpl colormap name.'''
#    return plt.cm.get_cmap(name, n)
#cmap = get_cmap(10)

#Void,Virus,Apc,Apm,I,Thn,The,Tkn,Tke,C

plot(c, N=150, out='cacovid_result1.pdf',
    vmax=9, graphic=True,
    colors=['white','red','white', 'darkblue', 'yellow', 'white', 'darkgreen', 'white', 'purple', 'brown'], 
    names=['Empty','Virus','Apc','Apm','Infected cells', 'Thn','The','Tkn','Tke','Proinflammatory cytokines'])
