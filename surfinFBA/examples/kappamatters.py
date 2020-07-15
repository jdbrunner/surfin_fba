
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import ode
import cobra as cb
from cobra import util

# import json
import pandas as pd


import sys
sys.path.append("../surfin_fba")
import surfinFBA as surf

# import call_me
import copy


import time
start_time = time.time()
from cycler import cycler
from datetime import datetime

import json

import os

import fcntl
import sys


try:
    os.mkdir("kappas_matter2")
except:
    pass


single_model_filenames = pd.read_csv('bigg_model_file_info.txt',dtype = str)
desired_models = ['E.coli','S.cerevisiae','M.tuberculosis']#['Ea','Pa','Pch','Pf','Pci','Pv','Pp','Sm']
endt = 2



def run_trio(cobra_models,endt,save_trajectories = False,allone = False, savenm = 'kappatest'):


    # if allone:
    #    random_numbers = 0.5*np.ones(10000)
    # else:
    #    random_numbers = np.random.rand(10000)


    my_models,metabolite_list,initial_metabolites = surf.prep_cobrapy_models(cobra_models)#returns dict,list,dict

    kappas = {}
    for mod in my_models.keys():
        kappas[mod] = dict([(metabolite_list[j],my_models[mod].uptakes[j]) for j in range(len(metabolite_list))])



    chc = list(my_models.keys())#list(np.array(xnsm)[np.array([0,1,2])])


    # for chc in xnchs:

    modelllist= [my_models[n] for n in chc]
    x0 = dict(zip(chc,np.ones(len(chc))))

    # print("Running simulation")
    ###USAGE: Surfin_FBA(model_list,x0,y0,metabolite_inflow,metabolite_dilution,endtime)
    with open('kappas_matter2/logfl.txt','a') as logfl:
        x,y,v,t,usage = surf.Surfin_FBA(modelllist,x0,initial_metabolites,{},{},endt,metabolite_names = metabolite_list,concurrent = True,report_activity = True, detail_activity = False,solver = 'both',flobj = logfl)
    # print("Simulation complete")


    xj = copy.deepcopy(x)
    for xx in xj:
        xj[xx] = list(xj[xx])

    yj = copy.deepcopy(y)
    for yy in yj:
        yj[yy] = list(yj[yy])

    vj = copy.deepcopy(v)
    for vv in vj:
        vj[vv] = list(vj[vv])
        for i in range(len(vj[vv])):
            vj[vv][i] = list(vj[vv][i])

    usagej = copy.deepcopy(usage)
    for uu in usagej:
        # usagej[uu] = list(usagej[uu])
        for i in usagej[uu]:
            usagej[uu][i] = list(usagej[uu][i])

    fullthing = {'X':xj, 'Y':yj, 'V':vj, 'U':usagej, 'T':list(t)}


    outcome = {'species':{},'metabolites':{}}
    for bug in fullthing['X'].keys():
        outcome['species'][bug] = fullthing['X'][bug][-1]
    for food in fullthing['Y'].keys():
        outcome['metabolites'][food] = fullthing['Y'][food][-1]

    if save_trajectories:

        print("Making Plots")
        fig,ax = plt.subplots(2,1,figsize = (10,10),tight_layout = True)
        # ax[0].set_prop_cycle(cycler(color = ['green', 'red','blue']))


        labels1 = []
        labels2 = []


        for nm,tc in x.items():
            ax[0].plot(t,tc)
            labels1 +=[nm]
        ax[0].legend(labels1,prop={'size': 14})
        for nm,tc in y.items():
            ax[1].plot(t,tc)
            labels2 +=[nm]
        # ax[1].legend(labels2,prop={'size': 20})

        # svfgr = True
        # if svfgr:
        fig.savefig('simulations/'+savenm+'plot' + ''.join([nm[-5:] for nm in chc]).replace(" ",""))# + datetime.now().strftime("%B%d_%H%M"))
        fig.close()
        # else:
        #     plt.show()

        with open('simulations/'+savenm+'data'+''.join([nm[-5:] for nm in chc]).replace(" ","")+'.json','w') as handle:
            json.dump(fullthing,handle)

    return outcome, kappas, fullthing['T'][-1]


if len(sys.argv) > 1:
    numtodo = int(sys.argv[1])
else:
    numtodo = 0



m9file = "modelsfromBiGG/sample_models/m9med.csv"
m9media_DF = pd.read_csv(m9file)






cobra_models = {}

for mod in desired_models:
    if any(single_model_filenames.Species == mod):
        flnm = single_model_filenames.loc[single_model_filenames.Species == mod,'File'].iloc[0]
        cobra_models[mod] = cb.io.load_json_model(flnm)

        exrxns = [rxn.id for rxn in cobra_models[mod].reactions if 'EX' in rxn.id]
        m9dict = dict([(m9media_DF.loc[i,'exchange reaction id'],m9media_DF.loc[i,'initial metabolite value']) for i in m9media_DF.index if m9media_DF.loc[i,'exchange reaction id'] in exrxns])
        cobra_models[mod].name = mod
        cobra_models[mod].medium = m9dict




    else:
        print("Error: No model of species " + mod)

if numtodo == 0:
    result,kappas, lt = run_trio(cobra_models,endt, allone = True,save_trajectories = True,savenm = 'uniform')
    if lt > 0.99*endt:
        randkey = 'uniform'

        with open('kappas_matter2/results.txt', 'a') as outfile:
            # fcntl.flock(outfile, fcntl.LOCK_EX)
            outfile.write(randkey + '\t' + json.dumps(result) + '\n')
            # fcntl.flock(outfile, fcntl.LOCK_UN)

        with open('kappas_matter2/kappas.txt', 'a') as outfile:
            # fcntl.flock(outfile, fcntl.LOCK_EX)
            outfile.write(randkey + '\t' + json.dumps(kappas) + '\n')
            # fcntl.flock(outfile, fcntl.LOCK_UN)


else:
    randkey = ''.join([str(np.random.choice(list('abcdefg123456789'))) for n in range(10)])

    result,kappas,lt = run_trio(cobra_models,endt,save_trajectories = True, savenm = randkey)


    with open('kappas_matter2/results.txt', 'a') as outfile:
        # fcntl.flock(outfile, fcntl.LOCK_EX)
        outfile.write(randkey + '\t' + json.dumps(result) + '\n')
        # fcntl.flock(outfile, fcntl.LOCK_UN)

    with open('kappas_matter2/kappas.txt', 'a') as outfile:
        # fcntl.flock(outfile, fcntl.LOCK_EX)
        outfile.write(randkey + '\t' + json.dumps(kappas) + '\n')
        # fcntl.flock(outfile, fcntl.LOCK_UN)


    for j in range(1,numtodo):
        result,kappas,lt = run_trio(cobra_models,endt)
        if lt > 0.99*endt:
            randkey = ''.join([str(np.random.choice(list('abcdefg123456789'))) for n in range(10)])

            with open('kappas_matter2/results.txt', 'a') as outfile:
                # fcntl.flock(outfile, fcntl.LOCK_EX)
                outfile.write(randkey + '\t' + json.dumps(result) + '\n')
                # fcntl.flock(outfile, fcntl.LOCK_UN)

            with open('kappas_matter2/kappas.txt', 'a') as outfile:
                # fcntl.flock(outfile, fcntl.LOCK_EX)
                outfile.write(randkey + '\t' + json.dumps(kappas) + '\n')
                # fcntl.flock(outfile, fcntl.LOCK_UN)
