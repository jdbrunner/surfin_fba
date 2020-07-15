import numpy as np
import pandas as pd
import cobra as cb
from scipy.integrate import ode
import matplotlib.pyplot as plt
import json
from datetime import date
import os

sv = True
today = date.today().strftime("%B%d%Y")



try:
    direct = "direct_meth/" + today
    os.mkdir(direct)
except:
    for i in range(100):
        direct = "direct_meth/" + today + "_" + str(i)
        try:
            os.mkdir(direct)
            break
        except:
            pass



#xdot = x*growth rate
#ydot = -x*flux
def dfba_RHS(t,y,params):
    model = params[0]
    exchg_order = params[1] #list of rxn ids in media file
    kappas = params[2]
    deathrate = params[3]
    inflow = params[4]
    x = y[0]
    yvec = y[1:]
    #reset the media file with kappa*y
    tmpmed = model.medium
    for i in range(len(yvec)):
        if exchg_order[i] in tmpmed.keys():
            tmpmed[exchg_order[i]] = kappas[i]*yvec[i]
    #optimize with cobra
    model.medium = tmpmed
    try:
        optm = model.optimize()
        xd = x*optm.objective_value - x*deathrate
        yd = list(np.array([x*optm.fluxes[exrn] if exrn in tmpmed.keys() else 0 for exrn in exchg_order]) +  np.array(yvec)*np.array(inflow))
        dfba_RHS.optcount += 1
    except:
        xd = 0
        yd = [0 for i in exchg_order]
        dfba_RHS.failtimes += [t]
    return np.array([xd] + yd)



def com_dfba_RHS(t,y,params):
    models = params[0]
    exchg_order = params[1] #list of rxn ids in media file
    kappas = params[2]
    deathrates = params[3]
    inflow = params[4]
    xs = y[:len(models)]
    yvec = y[len(models):]
    xdts = []
    ydots = np.zeros(len(yvec))
    for i in range(len(models)):
        xx = xs[i]
        state_vec = np.array([xx] + list(yvec))
        rhs = dfba_RHS(t,state_vec,[models[i],exchg_order,kappas[i],deathrates[i],inflow])
        xdts += [rhs[0]]
        ydots = ydots + np.array(rhs[1:])
    return np.array([xdts + list(ydots)])





combos=[['E.coli','M.tuberculosis'],['E.coli','S.cerevisiae'],['S.cerevisiae','M.tuberculosis'],['E.coli','P.putida'],['E.coli','S.cerevisiae','M.tuberculosis'],['E.coli','S.cerevisiae','P.putida','M.tuberculosis']]
for desired_models in combos:
    model_file_info = pd.read_csv('bigg_model_file_info.txt',dtype = str)
    endt = 15

    cobra_models = {}

    for mod in desired_models:
        if any(model_file_info.Species == mod):
            flnm = model_file_info.loc[model_file_info.Species == mod,'File'].iloc[0]
            cobra_models[mod] = cb.io.load_json_model(flnm)
            if not cobra_models[mod].name:
                cobra_models[mod].name = model_file_info.loc[model_file_info.Species == mod,'Species'].iloc[0] + "_" + model_file_info.loc[model_file_info.Species == mod,'ID'].iloc[0]
        else:
            print("Error: No model of species " + mod)


    print("Loaded " + str(len(cobra_models)) + " models successfully")


    m9file = "modelsfromBiGG/sample_models/m9med.csv"
    kappa_fl = "modelsfromBiGG/sample_models/model_exchange.json"

    m9media_DF = pd.read_csv(m9file)
    with open(kappa_fl) as fl:
        kappa_vals = json.load(fl)



    #compute dfba
    # modl1 = cobra_models[desired_models[0]]
    # modl2 = cobra_models[desired_models[1]]
    ## Get exchange reactions from minimal medium (minimized by components)
    # initial_growth1 = modl1.slim_optimize()
    # min_med1 = cb.medium.minimal_medium(modl1,initial_growth1,minimize_components = True)
    # initial_growth2 = modl2.slim_optimize()
    # min_med2 = cb.medium.minimal_medium(modl2,initial_growth2,minimize_components = True)
    cust_media = {}
    modids = []
    kappas = []
    exchg_ids = []

    for model in desired_models:
        flnm = model_file_info.loc[model_file_info.Species == model,'File'].iloc[0]
        cobmod  = cb.io.load_json_model(flnm)
        exrxns = [rxn.id for rxn in cobmod.reactions if 'EX' in rxn.id]
        modid = model_file_info.loc[model_file_info.Species == model,'ID'].iloc[0]
        # kappas_used[model] = dict([(m9media_DF.loc[i,'exchange reaction id'],1.0) for i in m9media_DF.index if m9media_DF.loc[i,'exchange reaction id'] in kappa_vals[modid].keys()])
        # m9dict = dict([(m9media_DF.loc[i,'exchange reaction id'],(1.0*m9media_DF.loc[i,'initial metabolite value'])) for i in m9media_DF.index if m9media_DF.loc[i,'exchange reaction id'] in kappa_vals[modid].keys()])

        # kappas_used[model] = dict([(m9media_DF.loc[i,'exchange reaction id'],kappa_vals[modid][m9media_DF.loc[i,'exchange reaction id']]) for i in m9media_DF.index if m9media_DF.loc[i,'exchange reaction id'] in kappa_vals[modid].keys()])
        m9dict = dict([(m9media_DF.loc[i,'exchange reaction id'],(m9media_DF.loc[i,'initial metabolite value'])*kappa_vals[modid][m9media_DF.loc[i,'exchange reaction id']]) for i in m9media_DF.index if m9media_DF.loc[i,'exchange reaction id'] in kappa_vals[modid].keys()])
        cust_media[model] = m9dict
        cobra_models[model].medium = m9dict
        # modids += [modid]
        exchg_ids += list(m9dict.keys())

    exchg_ids = list(np.unique(exchg_ids))

    for m in desired_models:
        kappas += [[kappa_vals[modid][exch] if exch in kappa_vals[modid].keys() else 0 for exch in exchg_ids]]


    y_init = [m9media_DF.loc[m9media_DF.loc[:,'exchange reaction id'] == er,'initial metabolite value'].values[0] for er in exchg_ids]
    flow = []
    for i in y_init:
        if i == 'D-Glucose':
            flow += [0.2]
        elif i == 'O2 O2':
            flow += [0.2]
        else:
            flow += [0]

    # exchg_ids = list(np.unique(list(cust_media[desired_models[0]].keys()) + list(cust_media[desired_models[1]].keys())))
    # modl1.medium = cust_media[desired_models[0]]
    # modl2.medium = cust_media[desired_models[1]]

    # modid1 = model_file_info.loc[model_file_info.Species == desired_models[0],'ID'].iloc[0]
    # modid2 = model_file_info.loc[model_file_info.Species == desired_models[1],'ID'].iloc[0]


    # kappa1 = [kappa_vals[modid1][exch] if exch in kappa_vals[modid1].keys() else 0 for exch in exchg_ids]
    # kappa2 = [kappa_vals[modid2][exch] if exch in kappa_vals[modid2].keys() else 0 for exch in exchg_ids]

    # print('EX_ca2_e' in cust_media[desired_models[0]].keys())
    # print('EX_ca2_e' in cust_media[desired_models[1]].keys())
    # print('EX_ca2_e' in kappa_vals[modid1].keys())
    # print('EX_ca2_e' in kappa_vals[modid2].keys())


    # y_init = [cust_media[desired_models[0]][exchg_ids[i]]/kappa_vals[modid1][exchg_ids[i]] if exchg_ids[i] in cust_media[desired_models[0]].keys() else cust_media[desired_models[1]][exchg_ids[i]]/kappa_vals[modid2][exchg_ids[i]] for i in range(len(exchg_ids))]


    with open(direct + "/" + "_".join(desired_models) + "_log.txt",'w') as logfl:

        logfl.write(" & ".join(desired_models))



        x_init = [0.3]*len(desired_models)
        for meth in ["vode","zvode","lsoda","dopri5","dop853"]:

            logfl.write("Both:" + meth + "\n")

            #### intitialize ODE
            dfba2 = ode(com_dfba_RHS).set_integrator(meth)
            init_v = x_init + y_init
            dfba2.set_initial_value(init_v,0)
            dfba2.set_f_params([[cobra_models[desmod] for desmod in desired_models],exchg_ids,kappas,[0.2]*len(desired_models),flow])


            x = [np.array(x_init)]
            yv = [np.array(y_init)]
            t = [0]

            dfba_RHS.optcount = 0
            dfba_RHS.failtimes = []

            # x = dfba2.integrate(endt)

            while dfba2.t < endt and dfba2.successful():
                sl = dfba2.integrate(dfba2.t + 0.01)
                x += [sl[:len(desired_models)]]
                yv += [np.array(sl[len(desired_models):])]
                t += [dfba2.t]

            yv = np.array(yv)

            y2 = dict([(exchg_ids[i],yv.T[i]) for i in range(len(yv.T))])

            xd = dict([(desired_models[i],np.array(x).T[i]) for i in range(len(desired_models))])

            logfl.write(str(dfba_RHS.optcount) + '\n\n')
            # logfl.write(dfba_RHS.failtimes)

            if meth == "vode":
                fig,ax = plt.subplots(2,1,figsize = (10,10),tight_layout = True)
                if len(x) < 9:
                    ax[0].set_prop_cycle(cycler(color = ['green', 'red','blue','purple','cyan','deeppink','goldenrod','slategray']))



                labels1 = []
                labels2 = []


                for nm,tc in xd.items():
                    ax[0].plot(t,tc)
                    labels1 +=[nm]
                ax[0].legend(labels1,prop={'size': 14})
                for nm,tc in y2.items():
                    if nm == 'EX_glc__D_e' or nm == 'EX_o2_e':
                        ax[1].plot(t,tc)
                        labels2 +=[nm]
                ax[1].legend(labels2,prop={'size': 14})

                plt.savefig(direct + "/" + "_".join([st.replace('.','') for st in desired_models]))
