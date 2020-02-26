import numpy as np

try:
    import cplex as cp
except ImportError:
    pass

try:
    import gurobipy as gb
except ImportError:
    pass

from scipy.integrate import ode

import pandas as pd
import time
from joblib import Parallel, delayed, cpu_count



# [Gamma1ar,Gamma2ar,lilgamma,internal_lower_bounds,internal_upper_bounds,kappas,exchng_lower_bounds]

class SurfMod:
    def __init__(self,G1,G2,lilg,ilbs,iubs,kaps,elbs,Name = None,deathrate = 0):
        self.Gamma1 = np.array(G1)
        self.Gamma2 = np.array(G2)
        self.objective = np.array(lilg)
        self.intLB = np.array(ilbs)
        self.intUB = np.array(iubs)
        self.uptakes = np.array(kaps)
        self.exchgLB = np.array(elbs)
        if Name == None:
            self.Name = ''.join([str(np.random.choice(list('abcdefg123456789'))) for n in range(5)])
        else:
            self.Name = Name
        self.MatrixA =  np.concatenate([np.array(G1),-np.array(G1),np.eye(np.array(G1).shape[1]),-np.eye(np.array(G1).shape[1])],axis = 0)
        self.statbds = np.concatenate([-np.array(elbs),np.array(iubs),-np.array(ilbs)])#np.empty(0)
        self.deathrate = deathrate

    def prep_indv_model(self,initial_N,secondobj = [],report_activity = True, solver = 'gb',flobj = None):

        Gamma1 = self.Gamma1
        Gamma2 = self.Gamma2
        obje = self.objective
        low_int = self.intLB
        up_int = self.intUB
        alphas = self.uptakes
        low_exch = self.exchgLB
        MatrixA = self.MatrixA

        t1 = time.time()
        Gamma1 = Gamma1.astype(float)
        Gamma2 = Gamma2.astype(float)
        obje = obje.astype(float)
        low_int = low_int.astype(float)
        up_int = up_int.astype(float)####Should check to make sure all LB <= UB
        alphas = alphas.astype(float)
        low_exch = np.minimum(low_exch,alphas*initial_N)

        if solver == 'gb':




            # MatrixA = np.concatenate([Gamma1,-Gamma1,np.eye(Gamma1.shape[1]),-np.eye(Gamma1.shape[1])],axis = 0)
            upbds_exch = initial_N*alphas

            if report_activity:
                try:
                    flobj.write("prep_indv_model: initializing LP\n")
                except:
                    print("prep_indv_model: initializing LP")
            growth = gb.Model("growth")
            growth.setParam( 'OutputFlag', False )


            sparms = [growth.addVar(lb = - gb.GRB.INFINITY,ub = gb.GRB.INFINITY, name = "s" + str(i)) for i in range(MatrixA.shape[1])]
            growth.update()
            objv = gb.quicksum([a[0]*a[1] for a in zip(obje,sparms)])
            growth.setObjective(objv,gb.GRB.MAXIMIZE)

            bds_vec = np.concatenate([upbds_exch,-low_exch,up_int,-low_int])
            if report_activity:
                try:
                    flobj.write("prep_indv_model: Adding constraints\n")
                except:
                    print("prep_indv_model: Adding constraints")

            growth.addConstrs((gb.quicksum([MatrixA[i][l]*sparms[l] for l in range(len(sparms))]) <= bds_vec[i] for i in range(len(MatrixA))), name = 'LE')
            growth.addConstrs((gb.quicksum([Gamma2[i][l]*sparms[l] for l in range(len(sparms))]) == 0 for i in range(len(Gamma2))), name = 'Kernal')
            growth.update()

            if report_activity:
                try:
                    flobj.write("prep_indv_model: optimizing LP\n")
                    flobj.write("prep_indv_model: optimizing with " + str(len(growth.getConstrs())) + " constraints\n" )
                except:
                    print("prep_indv_model: optimizing LP")
                    print("prep_indv_model: optimizing with ",len(growth.getConstrs()) ," constraints" )
            growth.optimize()


            status = growth.status
            # if status == 2:
            #
            # print(status)

            statusdic = {1:"LOADED",2:"OPTIMAL",3:"INFEASIBLE",4:"INF_OR_UNBD",5:"UNBOUNDED"}
            if status in statusdic.keys():
                if report_activity:
                    try:
                        flobj.write("find_waves: LP Status: " +  statusdic[status] + '\n')
                    except:
                        print("find_waves: LP Status: ", statusdic[status])
            else:
                if report_activity:
                    try:
                        flobj.write("find_waves: LP Status: Other\n")
                    except:
                        print("find_waves: LP Status: Other")

            if status == 2:


                # wi = np.array([v.x for v in growth.getVars()])#growth.solution.get_values()
                val = growth.objVal

                if len(secondobj) != len(sparms):#if not given a valid second objective, minimize total flux
                    secondobj = -np.ones(len(sparms))


                growth.addConstr(objv == val)
                growth.update()
                newobj = gb.quicksum([a[0]*a[1] for a in zip(secondobj,sparms)])
                growth.setObjective(newobj,gb.GRB.MAXIMIZE)
                growth.update()
                growth.optimize()

                wi = np.array([v.x for v in growth.getVars()])



                # static2 = np.concatenate([-low_exch,up_int,-low_int])
                if report_activity:
                    minuts,sec = divmod(time.time() - t1, 60)
                    try:
                        flobj.write("prep_indv_model: Done in " + str(int(minuts)) + " minutes, " + str(sec) + " seconds.\n")
                    except:
                        print("prep_indv_model: Done in ",int(minuts)," minutes, ",sec," seconds.")


                # self.statbds = static2
                return wi#,(MatrixA,static2,alphas,Gamma1,Gamma2,obje,death)
            else:
                return "failed to prep"

        elif solver == 'cp':
    #
            MatrixA = np.concatenate([Gamma1,-Gamma1,np.eye(Gamma1.shape[1]),-np.eye(Gamma1.shape[1])],axis = 0).astype(float)
            Gamma2 = Gamma2.astype(float)
            upbds_exch = initial_N*alphas

            if report_activity:
                try:
                    flobj.write("prep_indv_model: initializing LP\n")
                except:
                    print("prep_indv_model: initializing LP")




            growth = cp.Cplex()


            sparms = ["s" + str(i) for i in range(MatrixA.shape[1])]
            s_lbs = [-cp.infinity]*MatrixA.shape[1]
            s_ubs = [cp.infinity]*MatrixA.shape[1]

            growth.variables.add(obj = obje, lb = s_lbs, ub = s_ubs, names = sparms)
            growth.objective.set_sense(growth.objective.sense.maximize)

            growth.set_results_stream(None)
            growth.set_warning_stream(None)

            if report_activity:
                try:
                    flobj.write("prep_indv_model: Adding constraints\n")
                except:
                    print("prep_indv_model: Adding constraints")

            bdtypes = np.array(['L']*len(MatrixA) + ['E']*len(Gamma2))
            bds_vec = np.concatenate([upbds_exch,-low_exch,up_int,-low_int,np.zeros(len(Gamma2))])


            g1p2 = [list(g) for g in MatrixA] + [list(g) for g in Gamma2]
            g1p2_wi = [[sparms, g] for g in g1p2]
            growth.linear_constraints.add(lin_expr = g1p2_wi, senses = bdtypes,  rhs = bds_vec)
            if report_activity:
                try:
                    flobj.write("prep_indv_model: optimizing LP\n")
                except:
                    print("prep_indv_model: optimizing LP")
            growth.solve()

            status = growth.solution.get_status()
            statusdic = {1:"OPTIMAL",3:"INFEASIBLE",4:"INF_OR_UNBD",2:"UNBOUNDED"}


            if status in statusdic.keys():
                if report_activity:
                    try:
                        flobj.write("find_waves: LP Status: " +  statusdic[status] + '\n')
                    except:
                        print("find_waves: LP Status: ", statusdic[status])
            else:
                if report_activity:
                    try:
                        flobj.write("find_waves: LP Status: Other\n")
                    except:
                        print("find_waves: LP Status: Other")



            if status == 1:

                # wi = np.array([growth.solution.get_values("s"+str(i)) for i in range(MatrixA.shape[1])])
                # wi2 = np.array(growth.solution.get_values()) This should be the same but you never freaking know.

                val = growth.solution.get_objective_value()

                if len(secondobj) != len(sparms):#if not given a valid second objective, minimize total flux
                    secondobj = -np.ones(len(sparms)).astype(float)
                else:
                    secondobj = np.array(secondobj).astype(float)

                new_const = [sparms,list(obje)]
                growth.linear_constraints.add(lin_expr = [new_const], senses = ['E'], rhs = [val])
                newobj = [(sparms[i],secondobj[i]) for i in range(len(sparms))]
                growth.objective.set_linear(newobj)
                growth.solve()

                wi = np.array([growth.solution.get_values("s"+str(i)) for i in range(MatrixA.shape[1])])




                # static2 = np.concatenate([-low_exch,up_int,-low_int])
                if report_activity:
                    minuts,sec = divmod(time.time() - t1, 60)
                    try:
                        flobj.write("prep_indv_model: Done in " + str(int(minuts)) + " minutes " + str(sec) + " seconds.\n")
                    except:
                        print("prep_indv_model: Done in ",int(minuts)," minutes, ",sec," seconds.")

                # self.statbds
                return wi#,(MatrixA,static2,alphas,Gamma1,Gamma2,obje,death)


            else:
                return "failed to prep"


        else:
            print("Please select solver Gurobi: 'gb' or CPlex: 'cp'")
            return "failed to prep"



def get_expr_coos(expr, var_indices):
    rw = np.zeros(len(var_indices))
    wherewhat = [(var_indices[expr.getVar(i)],expr.getCoeff(i)) for i in range(expr.size())]
    for tp in wherewhat:
        rw[tp[0]] = tp[1]
    return rw



def prep_cobrapy_models(models,uptake_dicts = {},extracell = 'e', random_kappas = "new"):

    #can provide metabolite uptake dictionary as dict of dicts {model_key1:{metabolite1:val,metabolite2:val}}

    from cobra import util

    if not isinstance(models,dict):
        modeldict = {}
        for mod in models:
            modeldict[mod.name] = mod

    metaabs = {}
    y0s = {}
    exrn = {}
    metabids = {}
    nametoid = {}
    nametorxnid = {}
    urts = {}
    if len(uptake_dicts) == 0:
        try:
            random_nums = np.load(random_kappas)
            loadedrand = True
        except:
            random_nums = np.empty(0)
            loadedrand = False
            if random_kappas == "ones":
                print("Will use uniform uptake parameters = 1")
            else:
                print("Will create random uptake")
        rand_str_loc = 0

    for modelkey in models.keys():
        model = models[modelkey]

        # exchng_metabolite_names = [met.name for met in model.metabolites if met.compartment == extracell]
        exchng_metabolite_ids =  [met.id for met in model.metabolites if met.compartment == extracell]
        exchng_reactions = []
        exchng_metabolite_names = []
        for met in exchng_metabolite_ids:
            exchng_metabolite_names += [model.metabolites.get_by_id(met).name]
            exchng_reactions += [rxn.id for rxn in model.metabolites.get_by_id(met).reactions if 'EX_' in rxn.id]
        nutrient_concentrations = {}
        if len(uptake_dicts) ==0:
            if loadedrand:
                if rand_str_loc < len(random_nums):
                    uptake_rate = random_nums[rand_str_loc:(rand_str_loc + len(exchng_reactions))]
                    rand_str_loc = rand_str_loc + len(exchng_reactions)
                    uptkdict = dict(zip(exchng_metabolite_names,uptake_rate))
                else:
                    random_nums = np.concatenate([random_nums,np.random.rand(len(exchng_reactions))])
                    uptake_rate = random_nums[rand_str_loc:(rand_str_loc + len(exchng_reactions))]
                    rand_str_loc = rand_str_loc + len(exchng_reactions)
                    uptkdict = dict(zip(exchng_metabolite_names,uptake_rate))
            else:
                if random_kappas == "ones":
                    random_nums = np.concatenate([random_nums,np.ones(len(exchng_reactions))])
                    uptake_rate = random_nums[rand_str_loc:(rand_str_loc + len(exchng_reactions))]
                    rand_str_loc = rand_str_loc + len(exchng_reactions)
                    uptkdict = dict(zip(exchng_metabolite_names,uptake_rate))
                else:
                    random_nums = np.concatenate([random_nums,np.random.rand(len(exchng_reactions))])
                    uptake_rate = random_nums[rand_str_loc:(rand_str_loc + len(exchng_reactions))]
                    rand_str_loc = rand_str_loc + len(exchng_reactions)
                    uptkdict = dict(zip(exchng_metabolite_names,uptake_rate))

        else:
            uptkdict = uptake_dicts[modelkey]
            uptake_rate = [uptkdict[met] for met in exchng_metabolite_names]


        i = 0


        for er in exchng_reactions:
            al = uptake_rate[i]
            i += 1
            if er in model.medium.keys():
                nutrient_concentrations[er] = model.medium[er]/(al*500)
            else:
                nutrient_concentrations[er] = 0
            # uptake_rate+= [al]



        i = 0
        nmid = dict(zip(exchng_metabolite_ids,exchng_metabolite_names))
        idnm = dict(zip(exchng_metabolite_names,exchng_metabolite_ids))
        rxnid = dict(zip(exchng_reactions,exchng_metabolite_ids))
        nmtorxn = dict(zip(exchng_metabolite_names,exchng_reactions))
        y_init = dict([(nmid[rxnid[ky]],nutrient_concentrations[ky]) for ky in nutrient_concentrations.keys()])

        metaabs[model.name] = exchng_metabolite_names
        y0s[model.name] = y_init
        exrn[model.name] = exchng_reactions
        metabids[model.name] = exchng_metabolite_ids
        nametoid[model.name] = idnm
        nametorxnid[model.name] = nmtorxn
        urts[model.name] = uptkdict



    ##### NOW: we have to reconcile the exchanged metabolites. Swapping order means swapping rows of Gamma1! So
    ### we must agree on an order.
    masterlist = []
    for li in metaabs.values():
        masterlist += li
    masterlist = np.unique(masterlist)



    #### Initial y is not as difficult. Average them out.
    mastery0 = {}
    for nm in masterlist:
        yyy0 = 0
        ctt = 0
        for mod in y0s.values():
            if nm in mod.keys():
                yyy0 += mod[nm]
                ctt += 1
        if ctt:
            mastery0[nm] = yyy0/ctt
        else:
            mastery0[nm] = 0





    real_model = {}
    namemap ={}
    for modelkey in models.keys():
        model = models[modelkey]

        #Get the stoichiometric matrix and break it apart
        ###Index is metabolite ID, columns are rxn ID
        Gamma = util.array.create_stoichiometric_matrix(model, array_type = 'DataFrame')


        for meta in masterlist:
            if meta not in metaabs[model.name]:
                blnk = pd.DataFrame([np.zeros(len(Gamma.columns))],columns = Gamma.columns, index = [meta])
                Gamma = Gamma.append(blnk)
                # print(meta)
            elif nametoid[model.name][meta] not in Gamma.index:
                blnk = pd.DataFrame([np.zeros(len(Gamma.columns))],columns = Gamma.columns, index = [nametoid[model.name][meta]])
                Gamma = Gamma.append(blnk)


        mastertoids = [nametoid[model.name][nm] if nm in nametoid[model.name].keys() else nm for nm in masterlist]



        internal_reactions = np.array(Gamma.columns)[[rxn not in exrn[model.name] for rxn in Gamma.columns]]


        internal_metabs = np.array(Gamma.index)[[((met not in metabids[model.name]) and (met not in masterlist)) for met in Gamma.index]]




        Gamma1 = Gamma.loc[np.array(mastertoids),internal_reactions]

        # =============================================================================


        Gamma2 = Gamma.loc[internal_metabs,internal_reactions]
        Gamma1ar = Gamma1.values
        Gamma2ar = Gamma2.values



        #Next we need the objective function that identifies growth - the flux that COBRA optimizes
        growth_col = pd.Series(np.zeros(len(internal_reactions)),index = Gamma2.columns)
        biom_rxns = [rxn.id for rxn in util.solver.linear_reaction_coefficients(model).keys()]
        growth_col.loc[biom_rxns] = 1
        lilgamma = growth_col.values


        real_reactions = [ky for ky in nametorxnid[model.name].keys() if model.reactions.has_id(nametorxnid[model.name][ky])]


        exchng_lower_bounds = np.array([-model.reactions.get_by_id(nametorxnid[model.name][nm]).bounds[1] if nm in real_reactions else 0 for nm in masterlist])

        internal_upper_bounds = np.array([rxn.bounds[1] for rxn in model.reactions if rxn.id not in exrn[model.name]])
        internal_lower_bounds = np.array([rxn.bounds[0] for rxn in model.reactions if rxn.id not in exrn[model.name]])


        kappas = np.array([urts[model.name][nm] if nm in urts[model.name].keys() else 0 for nm in masterlist])

        real_model[modelkey] = SurfMod(Gamma1ar,Gamma2ar,lilgamma,internal_lower_bounds,internal_upper_bounds,kappas,exchng_lower_bounds,Name = model.name)
        #[Gamma1ar,Gamma2ar,lilgamma,internal_lower_bounds,internal_upper_bounds,kappas,exchng_lower_bounds]
        # namemap[modelkey] = model.name


    return real_model,masterlist,mastery0




def unpack_model(y,num_orgs,num_mets):
    x = y[:num_orgs]
    N = y[num_orgs:num_mets+num_orgs]
    return np.array(x),np.array(N)#,w

def pack_model(x,N):
    y = list(x) + list(N)
    return np.array(y)



def find_waves(y,v,bddts,surfmodel,headsup = [], model = None,report_activity = True, solver = 'gb',flobj = None):
    '''
    bddts must be 0 (for internal bound) or kappa_j ydot_j

    MatA is
    [Gamma1 ]
    [-Gamma1]
    [   I   ]
    [  -I   ]
    where Gamma1 is GammaStar and Gamma2 is GammaDagger

    statbds needs to be (negative) lower exchange then upper internal then (negative) lower internal


    Returns a basis matrix B and ind so vdot is solution to solve(B,bddts[ind])
    '''

    t = time.time()
    precision = 8
    minprecision = 0
    MatA = surfmodel.MatrixA
    statbds = surfmodel.statbds
    kappas = surfmodel.uptakes
    Gamma1 = surfmodel.Gamma1
    Gamma2 = surfmodel.Gamma2
    obje = surfmodel.objective
    death = surfmodel.deathrate
    ###Determine which constraints we are at
    bds = np.concatenate([y*kappas,statbds])

    dists = np.dot(MatA,v) - bds
    if len(headsup):
        dists[headsup] = 0.0
    atthem = ([],)#

    LP_not_ready = True

    if report_activity:
        try:
            flobj.write("find_waves: Initializing LP\n")
        except:
            print("find_waves: Initializing LP")

    if solver == 'gb':

        if model == None:

            findbasis = gb.Model("findbasis")
            findbasis.setParam( 'OutputFlag', False )

            if report_activity:
                try:
                    flobj.write("find_waves: Adding Variables\n")
                except:
                    print("find_waves: Adding Variables")
                taddk = time.time()
            vdots = [findbasis.addVar(lb = - gb.GRB.INFINITY,ub = gb.GRB.INFINITY, name = "s" + str(i)) for i in range(MatA.shape[1])]
            findbasis.update()
            objv = gb.quicksum([a[0]*a[1] for a in zip(obje,vdots)])
            findbasis.setObjective(objv,gb.GRB.MAXIMIZE)
            findbasis.addConstrs((gb.quicksum([Gamma2[i][l]*vdots[l] for l in range(len(vdots))]) == 0 for i in range(len(Gamma2))), name = 'Kernal')
            if report_activity:
                try:
                    flobj.write("find_waves: Added Kernal constraints in " + str(time.time() - taddk) + " seconds.\n")
                except:
                    print("find_waves: Added Kernal constraints in ",time.time() - taddk, " seconds.")
            inthere = []

        else:
            if report_activity:
                try:
                    flobj.write("find_waves: Removing Constraints\n")
                except:
                    print("find_waves: Removing Constraints")
            findbasis = model
            findbasis.reset(1)
            vdots = findbasis.getVars()
            inthere = [con.ConstrName for con in [findbasis.getConstrByName('LE[' + str(i) + ']') for i in range(len(bddts))] if con != None]

    elif solver == 'cp':
        if model == None:

            findbasis = cp.Cplex()
            if report_activity:
                try:
                    flobj.write("find_waves: Adding Variables\n")
                except:
                    print("find_waves: Adding Variables")
                taddk = time.time()
            vdots = ["s" + str(i) for i in range(MatA.shape[1])]
            lbs = [-cp.infinity]*len(vdots)
            ups = [cp.infinity]*len(vdots)
            findbasis.variables.add(obj = obje,lb = lbs,ub = ups,names = vdots)
            findbasis.objective.set_sense(findbasis.objective.sense.maximize)

            findbasis.set_results_stream(None)
            findbasis.set_warning_stream(None)

            constrmat = [[vdots,list(g2)] for g2 in Gamma2.astype(float)]
            findbasis.linear_constraints.add(lin_expr = constrmat, senses = ['E']*len(Gamma2), rhs = [0.0]*len(Gamma2), names = ['Kernal['+str(i)+']' for i in range(len(Gamma2))])

            inthere = []

        else:
            if report_activity:
                try:
                    flobj.write("find_waves: Removing Constraints\n")
                except:
                    print("find_waves: Removing Constraints")
            findbasis = model
            vdots = findbasis.variables.get_names()
            inthere = [nm for nm in findbasis.linear_constraints.get_names() if 'LE' in nm]

    else:
        print("Please select solver Gurobi: 'gb' or CPlex: 'cp'")
        return None


    while (LP_not_ready) and (precision >= minprecision):

        atthem = np.where(dists.round(precision) >= 0)


        the_rows_of_A = MatA[atthem]


        A_res_and_G2 = np.concatenate([the_rows_of_A,Gamma2])



        if np.linalg.matrix_rank(A_res_and_G2) < MatA.shape[1] and precision > minprecision:
            precision = precision - 1

        elif np.linalg.matrix_rank(A_res_and_G2) < MatA.shape[1]:
            print('find_waves: Precision = ', precision)
            print('find_waves: Returning BAD LP: Could not get enough bounds to construct LP')
            return 'BAD LP'

        else:




            ##Now we need to find a basis for the rows of this includes rows Gamma2
            ### (actually a basis for the row space of Gamma2) and such that
            ### if we solve Bv = bddts[B], then we get A_jv \leq bddts[j]
            #### We can get that by solving a LP



            if report_activity:
                try:
                    flobj.write("find_waves: Removing/Adding Constraints\n")
                except:
                    print("find_waves: Removing/Adding Constraints")
                taddrmv = time.time()

            constrnames = {} #Name of constraint to index in my problem
            for ii in range(MatA.shape[0]):
                constrnames['LE[' + str(ii) + ']'] = ii
            for jj in range(Gamma2.shape[0]):
                constrnames['Kernal[' + str(jj) + ']'] = jj + MatA.shape[0]

            dontrmv = ['LE[' + str(i) + ']' for i in atthem[0]]

            if solver == 'gb':

                for nm in inthere:
                    findbasis.getConstrByName(nm).RHS = bddts[constrnames[nm]]

                findbasis.update()



                findbasis.remove([con for con in [findbasis.getConstrByName(conNm) for conNm in inthere if conNm not in dontrmv]])
                findbasis.update()



                findbasis.addConstrs((gb.quicksum([MatA[i][l]*vdots[l] for l in range(len(vdots))]) <= bddts[i] for i in atthem[0] if (not ('LE[' + str(i) + ']' in inthere))), name = 'LE')

                findbasis.update()

                inthere = [con.ConstrName for con in [findbasis.getConstrByName('LE[' + str(i) + ']') for i in range(len(bddts))] if con != None]


                if report_activity:
                    try:
                        flobj.write("find_waves: Added new constraints in " + str(time.time()-taddrmv) + " seconds.\n")
                        flobj.write("find_waves: optimizing with " + str(findbasis.NumConstrs) + " constraints\n" )
                    except:
                        print("find_waves: Added new constraints in ",time.time()-taddrmv, " seconds.")
                        print("find_waves: optimizing with ",findbasis.NumConstrs," constraints" )


            else:
                for nm in inthere:
                    findbasis.linear_constraints.set_rhs(nm,bddts[constrnames[nm]])
                findbasis.linear_constraints.delete([con for con in inthere if con not in dontrmv])

                needtoadd = [at for at in atthem[0] if (not ('LE[' + str(at) + ']' in inthere))]

                newconstrs = [[vdots,list(MatA.astype(float)[rw])] for rw in needtoadd]
                new_rhs = [bddts[i] for i in needtoadd]
                new_names = ['LE[' + str(i) + ']' for i in needtoadd]

                findbasis.linear_constraints.add(lin_expr = newconstrs, senses = ['L']*len(newconstrs), rhs = new_rhs, names = new_names)

                inthere =  [nm for nm in findbasis.linear_constraints.get_names() if 'LE' in nm]




                if report_activity:
                    try:
                        flobj.write("find_waves: Added new constraints in " + str(time.time()-taddrmv) + " seconds.\n")
                        flobj.write("find_waves: optimizing with " + str(findbasis.linear_constraints.get_num()) + " constraints\n" )
                    except:
                        print("find_waves: Added new constraints in ",time.time()-taddrmv, " seconds.")
                        print("find_waves: optimizing with ",findbasis.linear_constraints.get_num()," constraints" )



            if solver == 'gb':
                findbasis.optimize()

                findbasis.update()


                status = findbasis.status
                # if status == 2:
                #
                # print(status)

                statusdic = {1:"LOADED",2:"OPTIMAL",3:"INFEASIBLE",4:"INF_OR_UNBD",5:"UNBOUNDED"}

                if report_activity:
                    if status in statusdic.keys():
                        try:
                            flobj.write("find_waves: LP Status: " + statusdic[status] + '\n')
                        except:
                            print("find_waves: LP Status: ", statusdic[status])
                    else:
                        try:
                            flobj.write("find_waves: LP Status: Other")
                        except:
                            print("find_waves: LP Status: Other")




                if statusdic[status] == "INFEASIBLE":
                    if report_activity:
                        try:
                            flobj.write('infeasible, will need to lower step size\n')
                        except:
                            print('infeasible, will need to lower step size')
                    return 'infeasible'


                if statusdic[status] == "INF_OR_UNBD":
                    # print('STAT4')
                    findbasis.reset()
                    findbasis.Params.DualReductions = 0
                    findbasis.optimize()

                    status = findbasis.status
                    # print(status)

                    if report_activity:
                        if status in statusdic.keys():
                            try:
                                flobj.write("find_waves: LP Status: " + statusdic[status] + '\n')
                            except:
                                print("find_waves: LP Status: ", statusdic[status])
                        else:
                            try:
                                flobj.write("find_waves: LP Status: Other")
                            except:
                                print("find_waves: LP Status: Other")


                    if statusdic[status] == "INFEASIBLE":
                        return 'infeasible'


                if statusdic[status] == "UNBOUNDED":
                    if precision > minprecision:
                        findbasis.printStats()
                        precision = precision - 1
                    else:
                        print('find_waves: Precision = ', precision)
                        findbasis.printStats()
                        print('find_waves: Returning BAD LP: Could not get enough basis vectors')
                        return 'BAD LP'


                if (statusdic[status] != "OPTIMAL") and (statusdic[status] != "UNBOUNDED"):
                    # print('STATBAD')
                    print('find_waves: Precision = ', precision)
                    findbasis.printStats()
                    print('find_waves: Returning BAD LP - ', statusdic[status])
                    return 'BAD LP'

            else:
                findbasis.solve()

                statusdic = {1:"OPTIMAL",3:"INFEASIBLE",4:"INF_OR_UNBD",2:"UNBOUNDED"}
                status = findbasis.solution.get_status()

                if report_activity:
                    if status in statusdic.keys():
                        try:
                            flobj.write("find_waves: LP Status: " + statusdic[status] + '\n')
                        except:
                            print("find_waves: LP Status: ", statusdic[status])
                    else:
                        try:
                            flobj.write("find_waves: LP Status: Other")
                        except:
                            print("find_waves: LP Status: Other")

                if statusdic[status] == "INFEASIBLE":
                    if report_activity:
                        print('infeasible, will need to lower step size')
                    return 'infeasible'

                if statusdic[status] == "INF_OR_UNBD":
                    checkobj = [[vdots,list(obje)]]
                    findbasis.linear_constraints.add([checkobj],senses = ['L'],rhs = [1], names = ['check'])
                    findbasis.solve()
                    chckstatus = findbasis.solution.get_status()
                    if chckstatus == 1:
                        status = 2
                    else:
                        status = 3

                    if report_activity:
                        if status in statusdic.keys():
                            try:
                                flobj.write("find_waves: LP Status: " + statusdic[status] + '\n')
                            except:
                                print("find_waves: LP Status: ", statusdic[status])
                        else:
                            try:
                                flobj.write("find_waves: LP Status: Other")
                            except:
                                print("find_waves: LP Status: Other")


                    if statusdic[status] == "INFEASIBLE":
                        return 'infeasible'

                if statusdic[status] == "UNBOUNDED":
                    if precision > minprecision:
                        findbasis.printStats()
                        precision = precision - 1
                    else:
                        print('find_waves: Precision = ', precision)
                        print('find_waves: Returning BAD LP: Could not get enough basis vectors')
                        return 'BAD LP'

                if (statusdic[status] != "OPTIMAL") and (statusdic[status] != "UNBOUNDED"):
                    print('find_waves: Precision = ', precision)
                    print('find_waves: Returning BAD LP - ', statusdic[status])
                    return 'BAD LP'




            if statusdic[status] == "OPTIMAL":
                LP_not_ready = False
                allbddts = np.concatenate([bddts,np.zeros(len(Gamma2))])

                allem = np.concatenate([atthem[0],np.arange(len(MatA),len(MatA)+len(Gamma2))])
                #

                sparmsnames = {'s' + str(i):i for i in range(MatA.shape[1])}#variable name to variable index in MatA
                if solver == 'gb':
                    vard = {v.getAttr(gb.GRB.Attr.VarName):v for v in findbasis.getVars()}#variable name to variable object
                    vardict = {vard[v]:sparmsnames[v] for v in sparmsnames.keys()}#variable object to variable index in MatA

                    rws = [(c.getAttr(gb.GRB.Attr.ConstrName),get_expr_coos(findbasis.getRow(c),vardict),c.getAttr(gb.GRB.Attr.Slack),c.getAttr(gb.GRB.Attr.RHS)) for c in findbasis.getConstrs() if c.getAttr(gb.GRB.Attr.CBasis) == -1]#tuples of row name, row, of basic rows (i.e. non-basic slacks)

                    the_basis_i_need = np.array([r[1] for r in rws])#non-basic constraints.



                    the_index_of_basis = [constrnames[r[0]] for r in rws]

                else:


                    # Returns a pair (head, x), where head is a list of variable indices and x is a list of floats indicating the values of those variables. Indices of basic slacks are specified by -rowindex - 1.
                    basic_head,basic_vals = findbasis.solution.basis.get_header()



                    allnms = findbasis.linear_constraints.get_names()

                    var_ord = findbasis.variables.get_names()
                    fullmat = np.array([findbasis.linear_constraints.get_coefficients([(cnm,s) for s in vdots]) for cnm in allnms])

                    nonbasvars = np.array([vdots[i] for i in range(fullmat.shape[1]) if i not in basic_head])#list of non-basic variables. need to put them in



                    all_slack_vals = np.zeros(len(allnms))
                    bslcks = []
                    for j in range(len(all_slack_vals)):
                        if -j-1 in basic_head:
                            all_slack_vals[j] = np.array(basic_vals)[np.where(np.array(basic_head) == -j-1)][0]
                            bslcks += [j]


                    #we may need to pivot in non-basic variables.
                    for ii in nonbasvars:
                        findbasis.advanced.pivot(ii,findbasis.advanced.no_variable,0)


                    basic_head,basic_vals = findbasis.solution.basis.get_header()
                    allnms = findbasis.linear_constraints.get_names()
                    #
                    var_ord = findbasis.variables.get_names()
                    #
                    nonbasvars = np.array([i for i in range(fullmat.shape[1]) if i not in basic_head])#list of non-basic variables. need to put them in
                    #
                    #
                    #
                    all_slack_vals = np.zeros(len(allnms))
                    bslcks = []
                    for j in range(len(all_slack_vals)):
                        if -j-1 in basic_head:
                            all_slack_vals[j] = np.array(basic_vals)[np.where(np.array(basic_head) == -j-1)][0]
                            bslcks += [j]


                    rws = [[allnms[i],findbasis.linear_constraints.get_coefficients([(allnms[i],s) for s in vdots])] for i in range(len(allnms)) if i not in bslcks]


                    the_basis_i_need = np.array([r[1] for r in rws])#non-basic constraints.



                    the_index_of_basis = [constrnames[r[0]] for r in rws]




                all_bddts = np.concatenate([bddts, np.zeros(len(Gamma2))])





                if np.linalg.matrix_rank(the_basis_i_need) < MatA.shape[1]:
                    if solver == 'gb':
                        if precision > minprecision and np.linalg.matrix_rank(np.concatenate([the_rows_of_A,Gamma2])) < MatA.shape[1]:
                            precision = precision - 1
                            LP_not_ready = True
                        else:
                            print('find_waves: Precision = ', precision)
                            print("find_waves: Returning BAD LP, not enough basis vectors")
                            return 'BAD LP'
                    else:
                        print("Missing dimensions?")
                        return 'BAD LP'



    if len(the_index_of_basis) > MatA.shape[1]:
        print('find_waves: Precision = ', precision)
        print("find_waves: Returning BAD LP, too many basis vectors.")
        return 'BAD LP'

    # print('find_waves: Done')
    if report_activity:
        try:
            flobj.write("find_waves: Found Basis, Computing QR Decomposition\n")
        except:
            print("find_waves: Found Basis, Computing QR Decomposition")
        for n in headsup:
            if n in the_index_of_basis:
                try:
                    flobj.write("find_waves: Constraint" + str(n) + " in basis.\n")
                except:
                    print("find_waves: Constraint", n," in basis.")
            else:
                if solver == 'gb':
                    modelvars = [findbasis.getVarByName('s' + str(i)) for i in range(MatA.shape[1])]
                    dvdtdotn = np.dot(MatA[n],np.array([v.x for v in modelvars]))
                else:
                    dvdtdotn = np.dot(MatA[n],np.array([findbasis.solution.get_values("s"+str(i)) for i in range(MatA.shape[1])]))
                try:
                    flobj.write("find_waves: Constraint " + str(n) + " - decreasing as " + str(bddts[n]) + " - dv/dt dot constraint n " + str(dvdtdotn) + "\n")
                except:
                    print("find_waves: Constraint ", n," - decreasing as ", bddts[n]," - dv/dt dot constraint n  ",dvdtdotn )


    QB,RB = np.linalg.qr(the_basis_i_need)
    if report_activity:
        try:
            flobj.write("find_waves: Done  in " + str(time.time() - t) +  " seconds.\n")
        except:
            print("find_waves: Done  in ",time.time() - t, " seconds.")

    return QB,RB,the_index_of_basis, findbasis


def compute_v(y,kappa,internal,G2len,basis):#Q,R,basisindex):
    Q,R,basisindex = basis
    ybds = kappa*y

    lbds_ex_min = np.minimum(ybds,-internal[:len(y)])
    upatlow_li = np.where(-internal[:len(y)] > lbds_ex_min)

    static = np.concatenate([internal,np.zeros(G2len)])

    static[upatlow_li] = -ybds[upatlow_li]

    bds = np.concatenate([ybds,static])


    QTb = np.dot(Q.T,bds[basisindex])
    v =np.linalg.solve(R,QTb)

    return v



def evolve_sys2(t,y,surfmod):##Computes xdot for an organism, and that organism's contribution to ydot
    '''
    MatA is
    [Gamma1 ]
    [-Gamma1]
    [   I   ]
    [  -I   ]
    where Gamma1 is GammaStar

    Gamma2 is GammaDagger

    statbds needs to be (negative) lower exchange then upper internal then (negative) lower internal
    '''
    # _,_,alphas,Gamma1,_,lilg,death,_,_,_ = params

    alphas = surfmod.uptakes#params[2]
    Gamma1 = surfmod.Gamma1#params[3]
    lilg = surfmod.objective#params[5]
    death = surfmod.deathrate#params[6]


    #Unpack state vector
    x = y[0]
    N = y[1:len(alphas) + 1]
    v = y[len(alphas)+1:]
    #compute xdot and ydot
    xd = np.array([x*(np.dot(v,lilg) - death)])
    Nd = -x*np.dot(Gamma1,v)

    return [xd,Nd]

def evolve_comm(t,y,all_params):

    '''Computes
        dx/dt, dy/dt, dv/dt
    '''

    model_list,bases,num_mets,v_sizes,metab_infl,metab_dil = all_params
    num_orgs = len(v_sizes)
    ###State vector is given as [x1,...,xn,y1,...,ym,v11,...,v1l,v21,...,v2o,...,vnp]
    unpk = unpack_model(y,num_orgs,num_mets)
    x = unpk[0]
    N = unpk[1]


    v = [compute_v(N,model_list[i].uptakes,model_list[i].statbds,len(model_list[i].Gamma2),bases[i]) for i in range(len(model_list))]



    xds = []
    Nds = np.zeros(len(N))
    for i in range(num_orgs):
        yds = evolve_sys2(t,np.concatenate([[x[i]],N,v[i]]),model_list[i])
        xds += list(yds[0])
        Nds = Nds + yds[1]
    Nds = Nds + metab_infl - metab_dil*N

    return np.concatenate([xds,Nds])

def Surfin_FBA(model_list,x0,y0,met_in,met_out,endtime,metabolite_names = [], report_activity = False, detail_activity = False, initres = 0.001,concurrent = True, solver = 'both',enoughalready = 10,flobj = None):
    '''
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

    '''

    if solver =='both':
        try:
            dir(cp)
            dir(gb)
            solver1 = 'cp'
            solver2 = 'gb'
        except:
            try:
                dir(gb)
                solver1 = 'gb'
                solver2 = 'gb'
                print('Cplex not found, using Gurobi')
            except:
                solver1 = 'cp'
                solver2 = 'cp'
                print('Gurobi not found, using Cplex')
    elif solver == 'gb':
        try:
            dir(gb)
            solver1 = 'gb'
            solver2 = 'gb'
        except:
            print('Gurobi not found, using Cplex')
            solver1 = 'cp'
            solver2 = 'cp'
    elif solver == 'cp':
        try:
            dir(cp)
            solver1 = 'cp'
            solver2 = 'cp'
        except:
            print('Cplex not found, using Gurobi')
            solver1 = 'gb'
            solver2 = 'gb'


    mess = 'Complete'

    t1 = time.time()

    ####Model names constrols the order of models!

    fwvtime = 0
    preptime = 0


    np.set_printoptions(linewidth = 200)

    if isinstance(model_list,dict):
        model_names = [mod.Name for mod in model_list.values()]#establish an order for the models
        model_keys = dict([(ky,mod.Name) for ky,mod in model_list.items()])
        model_keys_inverse = dict([(mod.Name,ky) for ky,mod in model_list.items()])

        model_list0 = [model_list[model_keys_inverse[nm]] for nm in model_names]#
        model_list = model_list0
        model_list0 = 0
    else:
        model_names = [mod.Name for mod in model_list]#or just get the names.
        model_keys = dict([(mod,mod) for mod in model_names])



    # if len(model_names) != len(x0):
    #     model_names = ['x' + str(i) for i in range(len(x0))]
    if len(metabolite_names) != len(y0):
        metabolite_names = ['y' + str(i) for i in range(len(y0))]


    ####Can pass in metabolite in/out as dicts.
    if isinstance(met_in,dict):
        metin = np.zeros(len(y0))
        for nm in met_in.keys():
            metin[np.where(np.array(metabolite_names) == nm)] = met_in[nm]
        met_in = metin
        metin = 0
    #
    if isinstance(met_out,dict):
        metout = np.zeros(len(y0))
        for nm in met_out.keys():
            metout[np.where(np.array(metabolite_names) == nm)] = met_out[nm]
        met_out = metout
        metout = 0

    ###Can pass in x0/y0 as dicts too I guess. Why not.
    if isinstance(x0,dict):
        xx0 = np.zeros(len(x0))
        for nm in x0.keys():
            xx0[np.where(np.array(model_names) == model_keys[nm])] = x0[nm]
        x0 = xx0
        xx0 = 0

    if isinstance(y0,dict):
        yy0 = np.zeros(len(y0))
        for nm in y0.keys():
            yy0[np.where(np.array(metabolite_names) == nm)] = y0[nm]
        y0 = yy0
        yy0 = 0




    # if isinstance(death,dict):
    #     dth = np.zeros(len(x0))
    #     for nm in death.keys():
    #         dth[np.where(np.array(model_names) == nm)] = death[nm]
    #     death = dth
    #     dth = 0

    #The models can be passed in a dict!




    x0 = np.array(x0)
    y0 = np.array(y0)
    death = np.array([mod.deathrate for mod in model_list])
    met_in = np.array(met_in)
    met_out = np.array(met_out)


    yd0 = np.zeros(len(y0))
    if report_activity:
        try:
            flobj.write('Surfin_FBA: Initializing with models\n'+', '.join(model_names) + '\n')
        except:
            print('Surfin_FBA: Initializing with models\n',model_names)

    lbds_li = [modl.exchgLB.astype(float) for modl in model_list]
    lbds_m_li = [np.minimum(modl.exchgLB,modl.uptakes*y0) for modl in model_list]

    static2_nox = []

    if report_activity:
        try:
            flobj.write('Surfin_FBA: Prepping models\n')
        except:
            print('Surfin_FBA: Prepping models')

    #Gamma1,Gamma2,lilgamma,lbds_int,upbds_int,kaps,lbds_ex = model_list[i]
    t_re = time.time()


    numjobs = min(len(model_list),cpu_count())

    if concurrent:
        preps = dict(Parallel(n_jobs=numjobs)(delayed(lambda i: (i,model_list[i].prep_indv_model(y0,report_activity = 0,solver = solver1,flobj = None)))(i) for i in range(len(model_list))))#preps contains the initial internal fluxes
    else:
        preps = dict([(i,model_list[i].prep_indv_model(y0,report_activity = 1,solver = solver1,flobj = flobj)) for i in range(len(model_list))])



    if "failed to prep" in [list(pva) for pva in preps.values()]:
        if report_activity:
            t2 = time.time() - t1
            minuts,sec = divmod(t2,60)
            try:
                flobj.write("Surfin_FBA: Failed to prep models  in " + str(int(minuts)) + " minutes, " + str(sec) + " seconds.\n")
            except:
                print("Surfin_FBA: Failed to prep models  in ",int(minuts)," minutes, ",sec," seconds.")
        return None,None,None,None,None


    yd0 = sum([-x0[i]*np.dot(model_list[i].Gamma1,preps[i]) for i in range(len(model_list))])
    # paramlist1 = [preps[i][1] for i in range(len(model_list))] #
    initial_vs = [preps[i] for i in range(len(model_list))]

    static2_nox = [np.concatenate([-modl.exchgLB,modl.intUB,-modl.intLB]) for modl in model_list]


    preptime = preptime + (time.time()-t_re)

    yd0 = yd0 + met_in - y0*met_out
    bases = []
    gb_modellist = []

    upperatlower = []
    for i in range(len(model_list)):


        t_re = time.time()

        lbds_ex = model_list[i].exchgLB.astype(float)#lbds_li[i]
        lbds_ex_min = np.minimum(model_list[i].exchgLB,model_list[i].uptakes*y0)#lbds_m_li[i]


        exup_mov= model_list[i].uptakes*yd0
        stati = np.zeros(len(model_list[i].MatrixA)-len(yd0))



        upperatlower += [np.where([lbds_ex[j] > lbds_ex_min[j] for j in range(len(lbds_ex))])]




        stati[upperatlower[i]] = -exup_mov[upperatlower[i]]

        bddts0 = np.concatenate([exup_mov,stati])

        if report_activity:
            try:
                flobj.write('Surfin_FBA: finding basis for ' + model_names[i] + '\n')
            except:
                print('Surfin_FBA: finding basis for ',model_names[i])

        fwreturn = find_waves(y0,initial_vs[i],bddts0,model_list[i],report_activity = detail_activity,solver = solver2,flobj = flobj)

        if ((fwreturn != 'BAD LP') and (fwreturn !='infeasible')):
            gb_modellist += [fwreturn[-1]]
            wavs = fwreturn[:-1]
            bases +=  [list(wavs)]
        else:
            break

        preptime = preptime + (time.time()-t_re)

    if ((fwreturn == 'BAD LP') or (fwreturn == 'infeasible')):
        if report_activity:
            t2 = time.time() - t1
            minuts,sec = divmod(t2,60)
            try:
                flobj.write("Surfin_FBA: Failed  in " + str(int(minuts)) + " minutes, " + str(sec) + " seconds.\n")
            except:
                print("Surfin_FBA: Failed  in ",int(minuts)," minutes, ",sec," seconds.")
        return None,None,None,None,None

    if report_activity:
        t2 = time.time() - t1
        minuts,sec = divmod(t2,60)
        try:
            flobj.write("Surfin_FBA: Models prepped and initial basis found in " + str(minuts) + " minutes, " + str(sec) + " seconds.\n")
        except:
            print("Surfin_FBA: Models prepped and initial basis found in ",minuts," minutes, ",sec," seconds.")




    v_sizes = np.array([len(vv) for vv in initial_vs])

    parameters = [model_list,bases,len(y0),v_sizes,met_in,met_out]

    ics = pack_model(x0,y0)

    dfba = ode(evolve_comm).set_integrator('vode')
    dfba.set_f_params(parameters)
    dfba.set_initial_value(ics,0)

    chk_round = 8
    # initres = 0.01: now a passed in option
    resolution = initres
    # enoughalready = 10 #how small to let resolution get before giving up. (-log) now a passed in option

    x = [x0]###list of arrays of xs
    y = [y0]###list of arrays of ys
    v = [initial_vs]### List of list of arrays of vs. So v[i][j][k] is v[k] for organism x[j] and time t[i]
    t = [0]

    ydparts = [[evolve_sys2(t[-1],np.concatenate([[x[-1][i]],y[-1],v[-1][i]]),model_list[i])[1] for i in range(len(model_list))]]### List of list of arrays of vs. So ydpart[i][j][k] is organism x[j]'s contribution to ydot[k] at time t[i]

    # tn = 0
    ok = [True]

    brokeit = [(np.array([]).astype(int),) for pl in model_list]


    notmoving = False
    reinits = 0
    check = True
    if report_activity:
        try:
            flobj.write('Surfin_FBA: Simulating Forward\n')
        except:
            print('Surfin_FBA: Simulating Forward')

    while dfba.t < endtime:

        full_state = dfba.integrate(dfba.t + resolution)
        unpacked = unpack_model(full_state,len(x0),len(y0))

        x_t = unpacked[0]
        y_t = unpacked[1]
        v_t = [compute_v(unpacked[1],model_list[i].uptakes,model_list[i].statbds,len(model_list[i].Gamma2),bases[i]) for i in range(len(model_list))]


###################################### Make sure the solution is still feasible
        ok = np.empty(len(model_list), dtype = bool)
        breakers = []
        upperatlower_new = []
        for j in range(len(model_list)):



            #Here, we can handle lower bounds greater than 0 by having them match the upper bounds. However, that means
            # they aren't changing smoothly, which isn't good.
            lbds_ex = model_list[j].exchgLB
            lbds_ex_min = np.minimum(lbds_ex,model_list[j].uptakes*y_t)
            upperatlower_new += [np.where(lbds_ex > lbds_ex_min)]
            no_longer = [num for num in upperatlower[j][0] if num not in upperatlower_new[j][0]]
            model_list[j].statbds[upperatlower_new[j]] = -(model_list[j].uptakes*y_t)[upperatlower_new[j]]
            model_list[j].statbds[(no_longer,)] = -lbds_ex[(no_longer,)]

            breakers += [np.dot(model_list[j].MatrixA,v_t[j]).round(chk_round) <= np.concatenate([model_list[j].uptakes*y_t,model_list[j].statbds]).round(chk_round)]

            ok[j] = np.all(breakers[j]) and np.all(np.dot(model_list[j].Gamma2,v_t[j]).round(chk_round) ==0)

        check = all(ok)

        ####Check tells us whether or not to reinitialize - if False then we need to. If true then just go on.

        if check:
            x += [x_t]
            y += [y_t]
            v += [v_t]
            t += [dfba.t]
            ydparts += [[evolve_sys2(t[-1],np.concatenate([[x[-1][i]],y[-1],v[-1][i]]),model_list[i])[1] for i in range(len(model_list))]]
            notmoving = False
            if resolution < initres:
                resolution = resolution*5
            upperatlower = upperatlower_new

        ############################## OTHERWISE find a new basis.
        else:
            reinits += 1
            ###We have to reinitialize
            t_re = time.time()

            ###Haven't managed to go anywhere so give up
            if (resolution < 10**(-enoughalready)):
                if report_activity:
                    try:
                        flobj.write('Surfin_FBA: Minimum step size reached at time t = ' + str(t[-1]) + '\n')
                        flobj.write('Surfin_FBA: Having trouble the following constraints: \n' + str(brokeit) + '\n')

                    except:
                        print('Surfin_FBA: Minimum step size reached at time t = ',t[-1])
                        print('Surfin_FBA: Having trouble the following constraints: \n',brokeit)
                mess = 'Failed:Min Step Size.'
                break

            if report_activity:
                try:
                    flobj.write("+++++++++++++++++++++++++++++++++++++++++++++\n")
                    flobj.write("+++++++++++++++++++++++++++++++++++++++++++++\n")
                    flobj.write('Surfin_FBA: Re-Initializing at system time {0}\n'.format(t[-1]))

                except:
                    print("+++++++++++++++++++++++++++++++++++++++++++++")
                    print("+++++++++++++++++++++++++++++++++++++++++++++")
                    print('Surfin_FBA: Re-Initializing at system time {0}'.format(t[-1]))

            ###What went wrong? Answer is stored in "breakers".


            #### Need to Re-Kajigger with previous timepoint
            yd0 = np.zeros(len(y0))
            for j in range(len(model_list)):
                ### Gamma^star_j is model_list[j].Gamma1
                yd0 += -x[-1][j]*np.dot(model_list[j].Gamma1,v[-1][j])
            yd0 = yd0 + met_in - met_out*y[-1]




            ###The actual basis finding step for each microbe:

            for j in range(len(model_list)):

                if not ok[j]:

                    exup_mov= model_list[j].uptakes*yd0
                    stati = np.zeros(len(model_list[j].MatrixA)-len(yd0))


                    stati[upperatlower[j]] = -exup_mov[upperatlower[j]]


                    bddts0 = np.concatenate([exup_mov,stati])



                    brokeit[j] = np.where(np.invert(breakers[j]))


                    fwreturn = find_waves(y[-1],v[-1][j],bddts0,model_list[j],headsup = brokeit[j][0], model = gb_modellist[j],report_activity = detail_activity,solver = solver2,flobj = flobj)



                    if fwreturn == 'infeasible':
                        notmoving = True

                    elif fwreturn != 'BAD LP':
                        wavs = fwreturn[:-1]
                        gb_modellist[j] = fwreturn[-1]
                        bases[j] =  list(wavs)


                    else:
                        ### Just stop
                        if report_activity:
                            try:
                                flobj.write('Surfin_FBA: Failed due to linear program unboundedness at system time t = ' + str(t[-1]) + '\n')
                            except:
                                print('Surfin_FBA: Failed due to linear program unboundedness at system time t = ',t[-1])
                        mess = 'Failed'
                        break


            if fwreturn == 'BAD LP':
                break

            #### Re-set the ODE and you're on your way!

            if notmoving:
                resolution = resolution*0.1



            parameters = [model_list,bases,len(y0),v_sizes,met_in,met_out]#[paramlist2,len(y0),v_sizes,met_in,met_out,[ml[-1] for ml in model_list]]

            ############################################################################################################
            ############################################################################################################
            ############################################################################################################


            ics = pack_model(x[-1],y[-1])
            dfba = ode(evolve_comm)
            dfba.set_f_params(parameters)
            dfba.set_initial_value(ics,t[-1])

            fwvtime = fwvtime + (time.time()-t_re)
            if report_activity:
                try:
                    flobj.write('Surfin_FBA: Simulating Forward\n')
                    flobj.write('Surfin_FBA: Time resultion:' + str(resolution) + '\n')
                except:
                    print('Surfin_FBA: Simulating Forward')
                    print('Surfin_FBA: Time resultion:',resolution)


            ##Count that we had to reinitialize - gets reset to 0 with a successful step forward
            notmoving = True





    x = np.array(x)
    biomasses = {}
    internal_flux = {}
    ydotconts = {}
    ydparts = np.array(ydparts)
    for i in range(len(x0)):
        biomasses[model_names[i]] = x[:,i]
        internal_flux[model_names[i]] = [vv[i] for vv in v]
        ydotconts[model_names[i]] = {}
        for j in range(len(y0)):
            ydotconts[model_names[i]][metabolite_names[j]] = ydparts[:,i,j]
    y = np.array(y)
    metabolite_bioms = {}
    for i in range(len(y0)):
        metabolite_bioms[metabolite_names[i]] = y[:,i]


    t2 = time.time() - t1
    minuts,sec = divmod(t2,60)
    try:
        flobj.write("Surfin_FBA: " +  mess+"  in " +str(int(minuts)) + " minutes, " + str(sec)+" seconds.\n")
        flobj.write("Surfin_FBA: Initialization was " + str(100*preptime/t2) + "% of computational time.\n")
        flobj.write("Surfin_FBA: Reinitialization was "+ str(100*fwvtime/t2)+ "% of computational time.\n")
        flobj.write("Surfin_FBA: Required " + str(reinits) + " reinitializations.\n")
    except:
        print("Surfin_FBA: ", mess, "  in ",int(minuts)," minutes, ",sec," seconds.")
        print("Surfin_FBA: Initialization was ", 100*preptime/t2, "% of computational time.")
        print("Surfin_FBA: Reinitialization was ", 100*fwvtime/t2, "% of computational time.")
        print("Surfin_FBA: Required ",reinits," reinitializations.")


    return biomasses,metabolite_bioms,internal_flux,t,ydotconts


def sim_cobraPY_comm(desired_models,model_info,x_init = {},death_rates = {},uptake_dicts = {},allinflow = 0,alloutflow = 0,met_inflow = {},met_outflow = {}, extracell = 'e', random_kappas = "new", save = False,save_fl = ''):
    '''
    paramters:

    * desired models = list - a list of keys for the models in the community to be simulated
    * model_info = dict - dictionary of model .json file paths.
    * x_init = dict - initial microbe biomasses, keyed by model keys. Any model not in dict will default to initial biomass 1
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



    cobra_models = {}

    for mod in desired_models:
        if mod in model_info.keys():
            flnm = model_info[mod]
            cobra_models[mod] = cb.io.load_json_model(flnm)
            if not cobra_models[mod].name:
                cobra_models[mod].name = mod
        else:
            print("Error: No model of species " + mod)


    print("Loaded " + str(len(cobra_models)) + " models successfully")




    my_models,metabolite_list,initial_metabolites = prep_cobrapy_models(cobra_models,uptake_dicts = uptake_dicts, extracell = extracell, random_kappas = random_kappas)

    print("Prepped all models successfully")

    x0 = dict(zip(desired_models,np.ones(len(desired_models))))
    for mod in desired_models:
        if mod in x_init.keys():
            x0[mod] = x_init[mod]

    for mod in desired_models:
        if mod in death_rates.keys():
            my_models[mod].deathrate = death_rates[mod]
        else:
            my_models[mod].deathrate = 0


    met_in = dict([(ky,allinflow) for ky in metabolite_list])
    for met in met_inflow.keys():
        met_in[met] = met_inflow[met]


    met_out = dict([(ky,alloutflow) for ky in metabolite_list])
    for met in met_outflow.keys():
        met_out[met] = met_outflow[met]


    print("Running simulation")
    ###USAGE: Surfin_FBA(model_list,x0,y0,met_in,met_out,endtime,model_names = [],metabolite_names = [],ptimes = True, report_activity = True, detail_activity = True, initres = 0.001,enoughalready = 10)
    with open("real_model_log.txt",'w') as logfl:
        x,y,v,t,usage = Surfin_FBA(my_models,x0,initial_metabolites,met_in,met_out,endt,metabolite_names = metabolite_list,concurrent = False,solver = 'both', flobj = logfl,report_activity = True, detail_activity = True)
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
    if save:
        fig.savefig(save_fl +'_fig_'+ ''.join(desired_models))
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

        with open(save_fl + '_data_' + ''.join(desired_models) + '.json','w') as handle:
                json.dump(fullthing,handle)

    else:
        plt.show()





    tottime = divmod(time.time()-start_time,60)
    print("--- %s minutes, %s seconds ---" % tottime)

    return x,y,v,t,usage
