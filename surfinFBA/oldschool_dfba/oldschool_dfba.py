from os.path import expanduser, join
from os import makedirs
import pkg_resources

import pandas as pd
import numpy as np
# import cobra as cb
#import cplex as cp

import matplotlib.pyplot as plt
import seaborn as sb

from cycler import cycler
import time








def dfba_oldschool(x0,y0,T,dt,model_li,alphas_li,V,infl,ddict):
    '''imliments dynamic fba for a community using the algorithm from Varma and Palsson
    give x0 and model_li as lists, ddict as dictionary of dilutions,infl as inflow dict
    '''
    t1 = time.time()



    mess = "Complete"
    num_sp = len(x0)
    if len(model_li) != num_sp:
        print('missing a model or two...')
        return None
    t = np.array([0])
    x = np.array([x0])
    x_dict = dict(zip([mod.name for mod in model_li],x0))
    xkys = list(x_dict.keys())
    ykys = list(y0.keys())
    y_out = np.array([[y0[yk] for yk in ykys]])
    y = y0.copy()
    media_keys = dict()
    rev_med_keys = dict()
    for mod in model_li:
        media_keys[mod.name] =  dict([(ky,mod.reactions.get_by_id(ky).name) for ky in list(mod.medium.keys())])
        rev_med_keys[mod.name] = dict([(mod.reactions.get_by_id(ky).name,ky) for ky in list(mod.medium.keys())])

    initime = time.time() - t1
    reoptime = 0
    optcount = 0
    #Reconcile the medium files with the exchanged metabolites - for each model we need a dict of
    #metabolites:reaction. Maybe pass that in with the initial conditions?
    while t[-1]<T:
        #update the media and optimize
        optimal_growths = dict()
        optimal_fluxes = dict()
        if any(np.array(list(y.values())) < 0):
            print('dfba_oldschool: overdepleated resource')
            mess = "Failed"
            break

        t3 = time.time()
        for mod in model_li:
            #set the media based on available metabolites
            tmp_med = mod.medium
            kydict = media_keys[mod.name]
            al_dict = alphas_li[mod.name]#need alphas_li to be a dict of dicts with keys model name and rxn_id
            for ky in tmp_med:
                tmp_med[ky] = al_dict[ky]*y[kydict[ky]]/V
            mod.medium = tmp_med
            #optimize
            modsol = mod.optimize()
            optimal_growths[mod.name] = modsol.objective_value
            optimal_fluxes[mod.name] = modsol.fluxes

        reoptime += time.time() - t3
        optcount += 1

        tmp_xd = {}
        for sp in x_dict:
            tmp_xd[sp] = x_dict[sp]*np.exp(dt*(optimal_growths[sp] - ddict[sp]))
        y_new = y
        for yi in y_new:
            yn = y_new[yi]
            for j in model_li:
                jnm = j.name
                if yi in rev_med_keys[jnm].keys():
                    if optimal_growths[jnm]-ddict[jnm] != 0:
                        yn += (-optimal_fluxes[jnm].loc[rev_med_keys[jnm][yi]]/(optimal_growths[jnm]-ddict[jnm]))*x_dict[jnm]*(1-np.exp(dt*(optimal_growths[sp]-ddict[sp])))
                    else:
                        yn += -dt*optimal_fluxes[jnm].loc[rev_med_keys[jnm][yi]]*x_dict[jnm]
            y_new[yi] = yn + infl[yi]*dt
        x_dict = tmp_xd
        t = np.append(t,[t[-1]+dt])
        x = np.append(x,[[x_dict[xk] for xk in xkys]],axis=0)
        y_out = np.append(y_out,[[y_new[yk] for yk in ykys]],axis=0)
        y = y_new


    t2 = time.time() - t1
    minuts,sec = divmod(t2,60)
    print("End t = ",t[-1])
    print("-----")
    print("-----")
    print("dfba_oldschool: ", mess, "  in ",int(minuts)," minutes, ",sec," seconds.")
    print("dfba_oldschool: Initialization was ", 100*initime/t2, "% of computational time.")
    print("dfba_oldschool: Reinitialization was ", 100*reoptime/t2, "% of computational time.")
    print("dfba_oldschool: Required ",optcount," reinitializations.")

    x_asdict = dict([(xkys[i],x[:,i]) for i in range(len(xkys))])
    y_asdict = dict([(ykys[i],y_out[:,i]) for i in range(len(ykys))])

    return x_asdict,y_asdict,t


def optimize_noncobra(toy,available):
    G1,G2,obj,inner_bounds,als,lwex = toy
    intk = als*available
    for i in range(len(lwex)):
        if lwex[i] > intk[i]:
            lwex[i] = intk[i]
    # print(lwex,intk)
    MatrixA = np.concatenate([G1,G1,G2],axis = 0)
    growth = cp.Cplex()
    growth.set_results_stream(None)
    growth.set_warning_stream(None)
    sparms = ['w'+str(j) for j in range(MatrixA.shape[1])]
    s_lbs = np.array(inner_bounds[0]).astype(float)
    s_ubs = np.array(inner_bounds[1]).astype(float)
    growth.variables.add(obj = obj, lb = s_lbs, ub = s_ubs, names = sparms)
    growth.objective.set_sense(growth.objective.sense.maximize)
    g1p2 = [list(g.astype(float)) for g in MatrixA]
    g1p2_wi = [[sparms, g] for g in g1p2]
    bds_vec = np.concatenate([lwex,intk,np.zeros(len(G2))])
    bdtypes = 'G'*len(G1)+'L'*len(G1)+'E'*len(G2)
    growth.linear_constraints.add(lin_expr = g1p2_wi, senses = bdtypes,  rhs = bds_vec)
    growth.solve()

    grate = growth.solution.get_objective_value()

    wi = growth.solution.get_values()

    ydots = np.dot(G1,wi)

    return [grate,ydots]



def dfba_oldschool_notcobra(x0,y0,T,dt,model_li,infl,dlist,V = 1):
    '''impliments dynamic fba for a community using the algorithm from Varma and Palsson
    give x0 and model_li as lists
    model_li should be lists of lists [G1,G2,obj,inner_bounds,als,exchng_lower_bounds]
    '''
    num_sp = len(x0)
    if len(model_li) != num_sp:
        print('missing a model or two...')
        return None
    t = np.array([0])
    x = np.array([x0])
    y = np.array([y0])
    while t[-1]<T:
        optimal_growths = np.empty(len(x[-1]))
        optimal_fluxes = np.empty((len(x[-1]),len(y[-1])))
        if any(y[-1] < 0):
            print('overdepleated resource')
            break
        j = 0
        for mod in model_li:
            # print(mod[5])
            gr,fl = optimize_noncobra(mod,y[-1])
            optimal_growths[j] = gr
            optimal_fluxes[j] = fl
            j+=1
        x_new = x[-1]*np.exp(dt*(optimal_growths-np.array(dlist)))
        mult_vect = dt*x[-1]
        for i in range(len(x[-1])):
            if optimal_growths[i]-dlist[i] != 0:
                mult_vect[i] = x[-1][i]/(optimal_growths[i]-dlist[i])*(1-np.exp(dt*(optimal_growths[i]-dlist[i])))
        y_new = y[-1] + np.dot(optimal_fluxes.T,mult_vect) + infl*dt
        t = np.append(t,[t[-1]+dt])
        x = np.append(x,[x_new],axis=0)
        y = np.append(y,[y_new],axis=0)
    return x,y,t
