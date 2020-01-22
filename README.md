Dynamic FBA for use with COBRApy metabolic models.

To use dFBA, first format the desired models using prep_cobrapy_models, which takes as input
a COBRApy model and outputs a model in the format used by surfin_fba.

prep_cobrapy_models(cobra_models: dictionary of models,, uptake_dicts = {}: dict of dicts, random_nums=random_numbers: list of saved random numbers (for debugging) for random uptake rates)#returns dict,list,dict

#can provide metabolite uptake dictionary as dict of dicts {model_key1:{metabolitename1:val,metabolitename2:val}}

next, run surfin_fba

Example usage:


my_models,metabolite_list,initial_metabolites,_ = surfin_fba.prep_cobrapy_models(cobra_models)

xnsm = list(my_models.keys())
modelllist= [my_models[n] for n in xnsm]
x0 = dict(zip(xnsm,np.ones(len(xnsm))))

death_rates = dict([(ky,0) for ky in xnsm])
met_in = dict([(ky,0) for ky in metabolite_list])
met_out = dict([(ky,0) for ky in metabolite_list])

###USAGE: Surfin_FBA(model_list,x0,y0,death,met_in,met_out,endtime,model_names = [],metabolite_names = [], report_activity = False, detail_activity = False, initres = 0.001,concurrent = True, solver = 'gb',enoughalready = 10,flobj = None)
with open("real_model_log.txt",'w') as logfl:
    x,y,v,t = surf.Surfin_FBA(modelllist,x0,initial_metabolites,death_rates,met_in,met_out,endt,model_names = xnsm,metabolite_names = metabolite_list,concurrent = False,solver = 'both', flobj = logfl,report_activity = True, detail_activity = True)

Surfin_FBA parameters:
model_list (positional) - list or dictionary of models as returned by prep_cobrapy_models

x0 (positional) - list or dictionary: initial microbial biomass

y0 (positional) - list or dictionary: initial metabolite biomass

death (positional) - list or dict: microbial dilution or death rates

met_in (positional) - list or dict: metabolite inflow rates

met_out (positional) - list or dict: metabolite outflow rates

endtime (positional) - int: simulation length

model_names (keyword) - list: model (microbe) names, which are keys for microbe related dicts. Establishes order

metabolite_names (keyword) - list: metabolite names, which are keys for metabolite related dicts. Establishes order

report_activity (keyword) - Bool: whether or not to log all simulation steps

detail_activity (keyword) - Bool: whether or not to log all simulation substeps, especially in basis finding

initres (keyword) - float: initial time step resolution

concurrent (keyword) - Bool: option to attempt to parallelize initialization

solver (keyword) - string: LP solver to use, 'cp' to use CPLEX, 'gb' to use Gurobi, 'both' to use a combination (fastest option)

enoughalready (keyword) - int: -log of minimum step size

flobj (keyword) - file-like: file onto which to print all simulation output (excluding some error messages), if None will print to stdout
