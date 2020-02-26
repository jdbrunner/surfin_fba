import numpy as np
import scipy as sp
# import cplex as cp
import matplotlib.pyplot as plt
from scipy.integrate import ode
import cobra as cb
# import json
import pandas as pd

import sys
import surfinFBA as surf


import time
start_time = time.time()
from cycler import cycler

from datetime import datetime



#### Microbe A (Brian):
#
#      y1 -1-> A -2-> B -3-
#                          -3/6->C -7->growth
#      y2 -4-> D -5-> E -6-
#                  5->y3
#

Gamma1A = -np.array([[-1,0,0,0,0,0,0],[0,0,0,-1,0,0,0],[0,0,0,0,1,0,0]])
Gamma2A = np.array([[1,-1,0,0,0,0,0],[0,1,-1,0,0,0,0],[0,0,1,0,0,1,-1],[0,0,0,1,-1,0,0],[0,0,0,0,1,-1,0]])



alsA = np.array([0.5,0.6,0.5])


lbds_exA = np.array([0,0,-100])
upbds_intA = np.array([100,100,100,100,100,100,100])
lbds_intA = np.array([0,0,0,0,0,0,0])

lilgammaA = np.array([0,0,0,0,0,0,1.0])


modelA = surf.SurfMod(Gamma1A,Gamma2A,lilgammaA,lbds_intA,upbds_intA,alsA,lbds_exA,Name = "Brian")


#### Microbe B (Dennis):
#
#      y1 -1-> A -2-> Growth
#
#
#      y3 -3-> B -4-> DEATH
#




#
Gamma1B = -np.array([[-1,0,0,0],[0,0,0,0],[0,0,-1,0]])
Gamma2B = np.array([[1,-1,0,0],[0,0,1,-1]])



alsB = np.array([0.7,0.6,2])


lbds_exB = np.array([0,0,10])####Need to figure out how to poison! Need positive lower bound, which will have to move with upper bound.
upbds_intB = np.array([100,100,100,100])
lbds_intB = np.array([0,0,0,0])

lilgammaB = np.array([0,1,0,-1]).astype(float)

modelB = surf.SurfMod(Gamma1B,Gamma2B,lilgammaB,lbds_intB,upbds_intB,alsB,lbds_exB,Name = "Dennis")



#### Microbe C (Carl):
#
#      y1 -1-> A -2-> B -3-> growth
#                  2->y2


Gamma1C = -np.array([[-1,0,0],[0,1,0],[0,0,0]])
Gamma2C = np.array([[1,-1,0],[0,1,-1]])



alsC = np.array([0.5,0.6,0.5])


lbds_exC = np.array([0,-100,0])
upbds_intC = np.array([100,100,100])
lbds_intC = np.array([0,0,0])

lilgammaC = np.array([0,0,1.0])

modelC = surf.SurfMod(Gamma1C,Gamma2C,lilgammaC,lbds_intC,upbds_intC,alsC,lbds_exC,Name = "Carl")


xA_init = 2
xB_init = 2
xC_init = 2

# y0 = [10,0,0]

###USAGE: Surfin_FBA(model_list,x0,y0,dilution_rates,metabolite_inflow,metabolite_dilution,endtime)

x0 = {'Brian':xA_init,'Dennis':xB_init,'Carl':xC_init}
y0 = {'Surfin USA':2,'Pet Sounds':0,'Good Vibrations':0}

x,y,v,t,usage = surf.Surfin_FBA([modelA,modelB,modelC],x0,y0,[1,1,1],[1,0,0],15,metabolite_names = ['Surfin USA','Pet Sounds','Good Vibrations'], detail_activity = 1, report_activity = 1, initres = 0.01,concurrent = False,solver = 'gb')



fig,ax = plt.subplots(5,1,figsize = (10,10),tight_layout = True)
ax[0].set_prop_cycle(cycler(color = ['green', 'red','blue']))


labels1 = []
labels2 = []


for nm,tc in x.items():
    ax[0].plot(t,tc)
    labels1 +=[nm]
ax[0].legend(labels1,prop={'size': 20})
for nm,tc in y.items():
    ax[1].plot(t,tc)
    labels2 +=[nm]
ax[1].legend(labels2,prop={'size': 20})
#
for met in y0:
    ax[2].plot(usage['Brian'][met],label =met)
    ax[3].plot(usage['Dennis'][met],label =met)
    ax[4].plot(usage['Carl'][met],label =met)
    ax[2].legend()
    ax[3].legend()
    ax[4].legend()



svfgr = False
if svfgr:
    fig.savefig('simulations/toy_community_' + datetime.now().strftime("%B%d%H%M"))
    fig.close()
else:
    plt.show()

print("--- %s seconds ---" % (time.time() - start_time))
