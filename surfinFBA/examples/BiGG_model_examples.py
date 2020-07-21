import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import ode
import cobra as cb
# import json
import pandas as pd
import sys
sys.path.append("../../")
import surfinFBA as surf

import copy
import json
# import call_me


import time
start_time = time.time()
from cycler import cycler
from datetime import datetime



endt = 1.0




single_model_filenames = pd.read_csv('bigg_model_file_info.txt',dtype = str)



desired_models = ['E.coli']#S.cerevisiae,C.difficile,H.pylori,P.putida,M.tuberculosis
cobra_models = {}

for mod in desired_models:
    if any(single_model_filenames.Species == mod):
        flnm = single_model_filenames.loc[single_model_filenames.Species == mod,'File'].iloc[0]
        cobra_models[mod] = cb.io.load_json_model(flnm)
        if not cobra_models[mod].name:
            cobra_models[mod].name = single_model_filenames.loc[single_model_filenames.Species == mod,'Species'].iloc[0] + "_" + single_model_filenames.loc[single_model_filenames.Species == mod,'ID'].iloc[0]
    else:
        print("Error: No model of species " + mod)


print("Loaded " + str(len(cobra_models)) + " models successfully")




my_models,metabolite_list,initial_metabolites = surf.prep_cobrapy_models(cobra_models)

print("Prepped all models successfully")

# xnsm = list(my_models.keys())
# modelllist= [my_models[n] for n in desired_models]
x0 = dict(zip(desired_models,np.ones(len(desired_models))))

death_rates = dict([(ky,1+0.2*np.random.rand()) for ky in desired_models])
for mod in desired_models:
    my_models[mod].deathrate = death_rates[mod]

met_in = dict([(ky,1) for ky in metabolite_list])
met_out = dict([(ky,1) for ky in metabolite_list])

print("Running simulation")
###USAGE: Surfin_FBA(model_list,x0,y0,met_in,met_out,endtime,model_names = [],metabolite_names = [],ptimes = True, report_activity = True, detail_activity = True, initres = 0.001,enoughalready = 10)
with open("bigg_model_log.txt",'w') as logfl:
    x,y,v,t,usage = surf.Surfin_FBA(my_models,x0,initial_metabolites,met_in,met_out,endt,metabolite_names = metabolite_list,concurrent = False,solver = 'both', flobj = logfl,report_activity = True, detail_activity = True)
print("Simulation complete")


print("Making Plots")
fig,ax = plt.subplots(2,1,figsize = (10,10),tight_layout = True)
ax[0].set_prop_cycle(cycler(color = ['green', 'red','blue']))


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

svfgr = False
if svfgr:
    fig.savefig('simulations/real_model_'+''.join(desired_models) + datetime.now().strftime("%B%d_%H%M"))
    fig.close()
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

    with open('simulations/real_model_data'+ ''.join(desired_models) + datetime.now().strftime("%B%d_%H%M")+'.json','w') as handle:
            json.dump(fullthing,handle)

else:
    plt.show()





tottime = divmod(time.time()-start_time,60)
print("--- %s minutes, %s seconds ---" % tottime)

# phone = "Real Model Examples Done, --- %s minutes, %s seconds ---" % tottime
#
# call_me.send(phone)
