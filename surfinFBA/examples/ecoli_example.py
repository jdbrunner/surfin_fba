import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import ode
import cobra as cb
# import json
import pandas as pd
import sys
sys.path.append("../../../surfin_fba")
# import Surfin_FBA as surf
import surfinFBA as surf


import copy
import json
# import call_me


import time
start_time = time.time()
from cycler import cycler



endt = 10.0




single_model_filenames = pd.read_csv('model_file_info.txt',dtype = str)
model_file_dict = dict([(single_model_filenames.loc[i,"Species"],single_model_filenames.loc[i,"File"]) for i in single_model_filenames.index])


desired_models = ["Ec"]


#### By default, the model uses the model.medium dictionary (averaging in the case of communities) to determine initial external metabolite concentration.
##### With this helper function, those files cannot be manipulated. However, we can pass a dictionary y_init = {metabolite:concentration} that
#### includes any initial concentration that we wish to change

### Any metaboite not in the media files that is marked as with extracellular by the model is assumed not initially present.

### Collect the metabolites in the media file for E. coli:

#Ecmod = cb.io.load_json_model(model_file_dict["Ec"])
#media_metabolites = surf.getMediaMetabolites(Ecmod)
#y0 = dict(zip(media_metabolites,100*np.ones(len(media_metabolites))))

y0 = {}
y0['O2'] = 0.24
y0['D-Glucose'] = 15.5
y0['D-Xylose'] = 8


#### Any metabolite marked "extracellular" by the model is exchanged.
Tup = surf.sim_cobraPY_comm(desired_models, model_file_dict,endt,y_init = y0,x_init = {'Ec':1})





