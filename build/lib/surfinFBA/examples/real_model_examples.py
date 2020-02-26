import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import ode
import cobra as cb
# import json
import pandas as pd
import sys
# sys.path.append("../../surfin_fba")
# import Surfin_FBA as surf
import surfinFBA as surf


import copy
import json
# import call_me


import time
start_time = time.time()
from cycler import cycler
from datetime import datetime



endt = 1.0




single_model_filenames = pd.read_csv('model_file_info.txt',dtype = str)
model_file_dict = dict([(single_model_filenames.loc[i,"Species"],single_model_filenames.loc[i,"File"]) for i in single_model_filenames.index])
print(model_file_dict)



desired_models = ["Ea"]#['Ea','Pa','Pch']#,'Pf','Pci','Pv','Pp','Sm']#['Pv','Pp']#['Ea']#["Cr","Sc"]#['Pa','Pf','Pci']#


# phone = "Real Model Examples Done, --- %s minutes, %s seconds ---" % tottime
#
# call_me.send(phone)


Tup = surf.sim_cobraPY_comm(desired_models, model_file_dict)
