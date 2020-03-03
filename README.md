Dynamic FBA for use with COBRAPy metabolic models.

To simulate a community of microbes given as .json CORBA models, you may use:

sim_cobraPY_comm(desired_models,model_info)

required (positional) paramters:

* desired models = list - a list of keys for the models in the community to be simulated
* model_info = dict - dictionary of model .json file paths, given as strings.

keyword parameters:

* x_init = dict - initial microbe biomasses, keyed by model keys. Any model not in dict will default to initial biomass 1
* y_init = dict - By default, the model uses the model.medium dictionary (averaging in the case of communities) to determine initial external metabolite concentration. With this helper function, those files cannot be manipulated. However, we can pass a dictionary y_init = {metabolite:concentration} that includes any initial concentration that we wish to change
* death_rates = dict - death/dilution rates of microbes. Defaults to 0
* uptake_dicts = {} - dict of dicts keyed by model key (from cobra_models dict) and metabolite. If empty, random parameters are generated
* allinflow = float - default metabolite inflow rate
* alloutflow = float - detault metabolite outflow rate
* met_inflow = dict - pass a dict to change only certain metabolite inflows
* met_outflow = dict - pass a dict to change only certain metabolite outflows.
* save = bool - if true, will save simulation plots and trajectory. Trajectory is saved as a .json file (load to python dictionary) with keys "X" for microbe biomass, "Y" for metabolite biomass, "V" for internal fluxes, "U" for metabolite usage, and "T" for time.
* save_fl = name of save files. Appended with '_fig_'+ ''.join(desired_models) + ".png" and '_data_' + ''.join(desired_models)  + ".json"
* extracell = string - name of extracellular compartment in COBRAPy model
* random_kappas = string - If this is a file containing random numbers (and so repeatable for debugging) these will be loaded and used for the uptake values kij. If this is "ones" then all uptake parameters will be set to 1. Otherwise, random numbers will be generated.

returns the tuple (biomasses,metabolite_bioms,internal_flux,t,ydotconts):
* biomasses = dict - microbe biomasses at each time point, keys are mod.Name for given models
* metabolite_bioms = dict - metabolite biomasses at each time point, keys are metabolite names
* internal_flux = dict - internal fluxes for each model, keys are mod.Name for given models
* t = list - time points of simulation
* ydotconts = dict of dicts - usage of each metabolite by species. Keys are mod.Name, second level keys are metabolite names
'''




The main function is

Surfin_FBA(model_list,x0,y0,met_in,met_out,endtime,metabolite_names = [], report_activity = False, detail_activity = False, initres = 0.001,concurrent = True, solver = 'both',enoughalready = 10,flobj = None)

Surfin_FBA parameters:
* model_list (positional) - list or dictionary of SurfMod objects (as returned by prep_cobrapy_models)
* x0 (positional) - list or dictionary: initial microbial biomass (keys should be mod.Name for each model in model_list)
* y0 (positional) - list or dictionary: initial metabolite biomass
* met_in (positional) - list or dict: metabolite inflow rates
* met_out (positional) - list or dict: metabolite outflow rates
* endtime (positional) - int: simulation length
* metabolite_names (keyword) - list: metabolite names, which are keys for metabolite related dicts.
* report_activity (keyword) - Bool: whether or not to log all simulation steps
* detail_activity (keyword) - Bool: whether or not to log all simulation substeps, especially in basis finding
* initres (keyword) - float: initial time step resolution
* concurrent (keyword) - Bool: option to attempt to parallelize initialization
* solver (keyword) - string: LP solver to use, 'cp' to use CPLEX, 'gb' to use Gurobi, 'both' to use a combination (fastest option)
* enoughalready (keyword) - int: -log of minimum step size
* flobj (keyword) - file-like: file onto which to print all simulation output (excluding some error messages), if None will print to stdout

Which returns the tuple (biomasses,metabolite_bioms,internal_flux,t,ydotconts):
* biomasses = dict - microbe biomasses at each time point, keys are mod.Name for given models
* metabolite_bioms = dict - metabolite biomasses at each time point, keys are metabolite names
* internal_flux = dict - internal fluxes for each model, keys are mod.Name for given models
* t = list - time points of simulation
* ydotconts = dict of dicts - usage of each metabolite by species. Keys are mod.Name, second level keys are metabolite names


A Genome scale metabolic model can be written as a stoichiometric matrix S which can be partitioned

[I G1]  
[0 G2]

Surfin_FBA defines the class SurfMod for simulation. This is a GEM of a single microbe re-formatted for dFBA simulation.

Instance attributes required are:
* Gamma1 = numpy array - part of stoichiometric matrix (G1)
* Gamma2 = numpy array - part of stoichiometric matrix (G2)
* objective = numpy array - growth objective for FBA
* intLB = numpy array - lower bounds on internal flux vectors (ordered)
* intUB = numpy array - upper bounds on internal flux vectors (ordered)
* uptakes = numpy array - uptake efficiencies [k1,...,kn] upper bound of exchange reaction is kij*yj where yj is the metabolite exchanged with microbe xi (this model).
* exchgLB = numpy array - lower bound of exchange reactions
* Name = string - model name
* deathrate = float - death or dilution rate of

The following attributes are created at initialization of an instance:
* MatrixA =  np.concatenate([G1,-G1,np.eye(G1.shape[1]),-np.eye(G1.shape[1])],axis = 0)
* statbds = np.concatenate([-elbs,iubs,-ilbs])#np.empty(0)

The method self.prep_indv_model(self,initial_N,secondobj = [],report_activity = True, solver = 'gb',flobj = None) further preps the model for use by Surfin_FBA. This is called by Surfin_FBA.

To use dFBA with COBRAPy model, first format the desired models using prep_cobrapy_models, which takes as input
a COBRApy model and outputs a model as a SurfMod object.

prep_cobrapy_models(cobra_models, uptake_dicts = {}, extracell = 'e', random_kappas = "new")

parameters:
* cobra_models = dict - dictionary of SurfMod objects
* uptake_dicts = {} - dict of dicts keyed by model key (from cobra_models dict) and metabolite. If empty, random parameters are generated
* extracell = string - name of extracellular compartment in COBRAPy model
* random_kappas = string - If this is a file containing random numbers (and so repeatable for debugging) these will be loaded and used for the uptake values kij. If this is "ones" then all uptake parameters will be set to 1. Otherwise, random numbers will be generated.

returns:
* dict of SurfMod objects
* list of exchanged metabolites (in an "agreed upon" order)
* dict of initial metabolite biomass
