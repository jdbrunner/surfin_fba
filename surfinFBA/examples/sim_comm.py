import os
import pandas as pd
import cobra as cb
import sys
sys.path.append("../../")
import surfinFBA as surf
from datetime import date
import json

try:
    mkdir("brunner_meth")
except:
    pass


sv = True
today = date.today().strftime("%B%d%Y")
try:
    direct = "brunner_meth/" + today
    os.mkdir(direct)
except:
    for i in range(100):
        direct = "brunner_meth/" + today + "_" + str(i)
        try:
            os.mkdir(direct)
            break
        except:
            pass

model_file_info = pd.read_csv('bigg_model_file_info.txt',dtype = str)
endt = 15

combos=[['E.coli'],['M.tuberculosis'],['S.cerevisiae'],['P.putida'],['E.coli','M.tuberculosis'],['E.coli','S.cerevisiae'],['S.cerevisiae','M.tuberculosis'],['E.coli','P.putida'],['E.coli','S.cerevisiae','M.tuberculosis'],['E.coli','S.cerevisiae','P.putida','M.tuberculosis']]#


#[['E.coli','S.cerevisiae']]#
m9file = "modelsfromBiGG/sample_models/m9med.csv"
kappa_fl = "modelsfromBiGG/sample_models/model_exchange.json"

m9media_DF = pd.read_csv(m9file)

# metin = {'D-Glucose':0.2,'O2 O2':0.2}


with open(kappa_fl) as fl:
    kappa_vals = json.load(fl)


for comm in combos:
    save_name = '-'.join([model_file_info.loc[model_file_info.Species == sp,'ID'].iloc[0] for sp in comm])


    desired_models = comm#['E.coli']#S.cerevisiae,C.difficile,H.pylori,P.putida,M.tuberculosis
    x_init = dict([(model,0.3) for model in desired_models])
    # death_rates = dict([(model,0.2) for model in desired_models])

    minmed = dict([(mod,'minimal') for mod in desired_models])
    cust_media = {}
    kappas_used = {}
    for model in desired_models:
        flnm = model_file_info.loc[model_file_info.Species == model,'File'].iloc[0]
        cobmod  = cb.io.load_json_model(flnm)
        exrxns = [rxn.id for rxn in cobmod.reactions if 'EX' in rxn.id]
        modid = model_file_info.loc[model_file_info.Species == model,'ID'].iloc[0]
        # kappas_used[model] = dict([(m9media_DF.loc[i,'exchange reaction id'],1.0) for i in m9media_DF.index if m9media_DF.loc[i,'exchange reaction id'] in kappa_vals[modid].keys()])
        # m9dict = dict([(m9media_DF.loc[i,'exchange reaction id'],(1.0*m9media_DF.loc[i,'initial metabolite value'])) for i in m9media_DF.index if m9media_DF.loc[i,'exchange reaction id'] in kappa_vals[modid].keys()])

        kappas_used[model] = dict([(m9media_DF.loc[i,'exchange reaction id'],kappa_vals[modid][m9media_DF.loc[i,'exchange reaction id']]) for i in m9media_DF.index if m9media_DF.loc[i,'exchange reaction id'] in kappa_vals[modid].keys()])
        m9dict = dict([(m9media_DF.loc[i,'exchange reaction id'],(m9media_DF.loc[i,'initial metabolite value'])*kappa_vals[modid][m9media_DF.loc[i,'exchange reaction id']]) for i in m9media_DF.index if m9media_DF.loc[i,'exchange reaction id'] in kappa_vals[modid].keys()])
        cust_media[model] = m9dict

        with open(direct+'/'+save_name+'m9media_kappas.json','w') as fl:
            json.dump(kappas_used,fl)




    surf.sim_cobraPY_comm(desired_models,model_file_info,endt,media = cust_media,x_init = x_init,y_init = {},death_rates = {},uptake_dicts = kappas_used,allinflow = 0,alloutflow = 0,met_inflow = {},met_outflow = {}, extracell = 'e', random_kappas = "ones", save = sv,save_fl = direct+'/'+save_name+'m9media',concurrent = False, solver = 'both',met_plots = ['D-Glucose','O2 O2'])

    # surf.sim_cobraPY_comm(desired_models,model_file_info,endt,media = cust_media,x_init = x_init,y_init = {},death_rates = {},uptake_dicts = {},allinflow = 0,alloutflow = 0,met_inflow = {},met_outflow = {}, extracell = 'e', random_kappas = "ones", save = sv,save_fl = direct+'/'+save_name+'m9media',concurrent = False, solver = 'both')
