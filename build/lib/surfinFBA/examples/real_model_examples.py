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



endt = 1.0




single_model_filenames = pd.read_csv('model_file_info.txt',dtype = str)
model_file_dict = dict([(single_model_filenames.loc[i,"Species"],single_model_filenames.loc[i,"File"]) for i in single_model_filenames.index])
print(model_file_dict)



desired_models = ["Ea"]#['Ea','Pa','Pch']#,'Pf','Pci','Pv','Pp','Sm']#['Pv','Pp']#['Ea']#["Cr","Sc"]#['Pa','Pf','Pci']#


#### By default, the model uses the model.medium dictionary (averaging in the case of communities) to determine initial external metabolite concentration.
##### With this helper function, those files cannot be manipulated. However, we can pass a dictionary y_init = {metabolite:concentration} that
#### includes any initial concentration that we wish to change

### Any metaboite not in the media files that is marked as with extracellular by the model is assumed not initially present.

### Collect the metabolites in the media file for Ea:

Eamod = cb.io.load_json_model(model_file_dict["Ea"])
media_metabolites = surf.getMediaMetabolites(Eamod)
y0 = dict(zip(media_metabolites,np.zeros(len(media_metabolites))))

#### Any metabolite marked "extracellular" by the model is exchanged.
Tup = surf.sim_cobraPY_comm(desired_models, model_file_dict,endt,y_init = y0)
