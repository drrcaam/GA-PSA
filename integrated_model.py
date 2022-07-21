# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 10:02:11 2022
Integrated model system.
@author: Owner
"""
import pickle
import gen_transit_param
import gen_ga_param
import fixed_transit_GA
import service_area_GA as sa_ga


### STEP 0. Generate parameters
transit_param = gen_transit_param.gen()
route_GA_param, sa_GA_param = gen_ga_param.GA_param()

### STEP 1. Generate fixed route designs
#Generate the set of routes
fixed_route_designs, design_id = fixed_transit_GA.routes_GA(route_GA_param, transit_param)


# ### STEP 2. Explore additional service area configurations
num_trials = 10
for trial in range(0,num_trials):

    #Gather the minimum service areas generated for the generated fixed route 
    #network designs
    best_ft_designs = {}
    all_upgraded_designs = {}
    for design_archive in fixed_route_designs.values():
        all_designs = design_archive['all_designs']
        parent_ids  = design_archive["parent_ids"] 
        
        counter = 0
        for id_ft_design in parent_ids:
            
            counter += 1
            best_ft_designs[id_ft_design] = all_designs[id_ft_design]
    
            #Search for better paratransit coverage using 
            #service area genetic algoritm
            all_upgraded_designs, design_id = sa_ga.sa_ga_post_gen(best_ft_designs[id_ft_design],id_ft_design, design_id,  transit_param, sa_GA_param, all_upgraded_designs)
    
    
            filehandler = open("upgraded_designs_tr_" +str(trial) + "_prt_" + str(counter) + ".obj", 'wb') 
            pickle.dump([all_upgraded_designs, design_id] , filehandler)
            filehandler.close()

