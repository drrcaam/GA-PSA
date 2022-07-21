# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 15:35:15 2022

@author: Owner
"""
import numpy as np

def non_dom_sorting(F):
    """
    Implements the Fast Nondominated Sorting Approach proposed by Deb et al.
    "A Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II"

    PARAMETERS:
        + F      -- 2D array with objective function values. Columns contain different
                    objective function values. Rows represent different designs.
    
    OUTPUT:
        + p_rank -- 1D array that contains the Pareto rank of each point.
    
    """
    num_sol, num_objectives = F.shape
    num_dominator = np.zeros((num_sol,))
    S_p = {}
    front = {1: []} #dictionary that stores the 
    p_rank = np.zeros((num_sol,), dtype='int64') #array to store rank of points
    
    for p in range(num_sol):
        #Find set of points dominated by point p
        mask = (np.sum(F[p,:] <= F ,axis = 1) == num_objectives) &\
               (np.sum(F[p,:] <  F ,axis = 1) >  0) 
        idx = mask.nonzero()[0]
        S_p[p] = idx
        
        #Find set of points that dominate p
        mask = (np.sum(F <= F[p,:] ,axis = 1) == num_objectives) &\
               (np.sum(F <  F[p,:] ,axis = 1) >  0) 
        num_dominator[p] = np.count_nonzero(mask)
        
        if num_dominator[p] == 0:
            p_rank[p] = 1
            front[1].append(p)
                   
    i = 1
    while len(front[i]) > 0:
        next_front = []
        for p in front[i]:
            idx = S_p[p]
            num_dominator[idx] = num_dominator[idx] - 1
            mask = num_dominator[idx] == 0
            p_rank[idx[mask]] = i + 1
            next_front += list(idx[mask])
        i += 1
        front[i] = next_front
        
    return p_rank
    
def crowding_distance(F):
    """
    Implements the crowding distance measure proposed by Deb et al.
    "A Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II"

    PARAMETERS:
        + F      -- 2D array with objective function values. Columns contain different
                    objective function values. Rows represent different designs.
    
    OUTPUT:
        + p_rank -- 1D array that contains the Pareto rank of each point.
    
    """    
    num_sol, num_objectives = F.shape
    distance = np.zeros((num_sol,))
    
    for m in range(num_objectives):
       
        idx = np.argsort(F[:,m])
        f_min = F[idx[0],m]
        f_max = F[idx[-1],m]
        
        if f_min == f_max:
            f_min = 0
            f_max = 1
        
        distance[idx[0]] = np.inf
        distance[idx[-1]] = np.inf
        
        for i in range(1, num_sol - 1):
            distance[idx[i]] += (F[idx[i+1], m] - F[idx[i-1], m])/(f_max - f_min)
        
    return distance


def nsga_rules(all_designs, pop_ids, population, OF, pop_size, budget_threshold):
    """
    Integrates the non_dom_sorting and crowding_distance functions.
    Adjust ranking of solutions that do not satisfy the budget constraint. 
    
    PARATEMERS:
        + all_designs      -- dictionary with all the designs
        + pop_ids          -- list with the ids of the designs being considered
        + population       -- 2D array with the route sets
        + OF               -- 2D array with objective function values. Columns contain different
                               objective function values. Rows represent different designs.
        + pop_size         -- number of designs that must be selected.
        + budget_threshold --
        
    OUTPUT:
        + selected_ids     -- ids of the designs that were selected
        + new_population   -- 2D array; each row contains the selected route sets
        + new_OF           -- 2D array; contains the objective function values of selected
                            designs
    """
    
    pop_ids = np.array(pop_ids)
    
    #Collect the remaining budget of all designs and determine if they
    #contain the original coverage:
    all_rem_budget = []
    # all_coverage = []
    for network_id in pop_ids:
        all_rem_budget.append(all_designs[network_id]["rem_budget"])
        # all_coverage.append(all_designs[network_id]["sa"].contains(og_coverage))
        
    all_rem_budget = np.array(all_rem_budget)
    # all_coverage = np.array(all_coverage)
    
    #Identifying Feasible solutions
    mask = (all_rem_budget >= budget_threshold) #& all_coverage #this zero 
    print("number of feasible solutions: ", np.sum(mask))
    
    #Creating array that will store pareto ranks
    p_rank = np.zeros((OF.shape[0],))
    
    #Finding and setting pareto ranks of feasible solutions
    feasible_rank = non_dom_sorting(OF[mask,:])
    
    p_rank[mask] = feasible_rank
    
    #Setting the pareto rank of infeasible solutions
    max_rank = max(feasible_rank)
    p_rank[mask == False] = max_rank + 1

    #Elitist selection
    for rank_limit in range(1,max_rank + 2):
        mask = p_rank <= rank_limit
        if np.sum(mask) == pop_size:
            selected_ids = pop_ids[mask]
            new_population = population[mask, :]
            new_OF = OF[mask,:]
        elif np.sum(p_rank <= rank_limit) > pop_size:
            selection_mask = p_rank < rank_limit #accepted solutions
            
            remainder = pop_size - np.sum(selection_mask) #need to find extra
            maybe_idx = (p_rank == rank_limit).nonzero()[0] #candidate extra solutions
            distances = crowding_distance(OF[maybe_idx,:]) #distances of extra
            saved = np.argsort(distances)[::-1][:remainder] #getting indeces of solutions that were not pruned
            selection_mask[maybe_idx[saved]] = True
            
            selected_ids = pop_ids[selection_mask]
            new_population = population[selection_mask, :]
            new_OF = OF[selection_mask,:]
            break

            
    return selected_ids, new_population, new_OF 

def nsga_rules_sa(all_designs, pop_ids, population_array, population_geo, OF, pop_size, budget_threshold):
    """
    Integrates the non_dom_sorting and crowding_distance functions.
    Adjust ranking of solutions that do not satisfy the budget constraint. 
    
    PARATEMERS:
        + all_designs          -- dictionary with all the designs
        + pop_ids              -- list with the ids of the designs being considered
        + population_array     --
        + population_geo       --
        + population           -- list of arrays with vertices that represent service array
        + OF                   -- 2D array with objective function values. Columns contain different
                                  objective function values. Rows represent different designs.
        + pop_size             -- number of designs that must be selected.
        + budget_threshold     --
        
    OUTPUT:
        + selected_ids         -- ids of the designs that were selected
        + new_population_array -- 
        + new_population_geo   --
        + new_population       -- 2D array; each row contains the selected route sets
        + new_OF               -- 2D array; contains the objective function values of selected
                                  designs
    """
    
    pop_ids = np.array(pop_ids)
    
    #Collect the remaining budget of all designs:
    all_rem_budget = []
    for network_id in pop_ids:
        all_rem_budget.append(all_designs[network_id]["rem_budget"])
    all_rem_budget = np.array(all_rem_budget)
    
    #Identifying Feasible solutions
    mask = all_rem_budget >= budget_threshold #this zero 
    
    #Creating array that will store pareto ranks
    p_rank = np.zeros((OF.shape[0],))
    
    #Finding and setting pareto ranks of feasible solutions
    feasible_rank = non_dom_sorting(OF[mask,:])
    
    p_rank[mask] = feasible_rank
    
    #Setting the pareto rank of infeasible solutions
    max_rank = max(feasible_rank)
    p_rank[mask == False] = max_rank + 1

    #Elitist selection
    new_population_array = []
    new_population_geo   = []
    for rank_limit in range(1,max_rank + 2):
        mask = p_rank <= rank_limit
        if np.sum(mask) == pop_size:
            selected_ids = pop_ids[mask]
            new_OF = OF[mask,:]
            
            for idx in mask.nonzero()[0]:
                new_population_array.append(population_array[idx])
                new_population_geo.append(population_geo[idx])
                
        elif np.sum(p_rank <= rank_limit) > pop_size:
            selection_mask = p_rank < rank_limit #accepted solutions
            
            remainder = pop_size - np.sum(selection_mask) #need to find extra
            maybe_idx = (p_rank == rank_limit).nonzero()[0] #candidate extra solutions
            distances = crowding_distance(OF[maybe_idx,:]) #distances of extra
            saved = np.argsort(distances)[::-1][:remainder] #getting indeces of solutions that were not pruned
            selection_mask[maybe_idx[saved]] = True
            
            selected_ids = pop_ids[selection_mask]
            new_OF = OF[selection_mask,:]
            
            for idx in selection_mask.nonzero()[0]:
                new_population_array.append(population_array[idx])
                new_population_geo.append(population_geo[idx])
                
            break

            
    return selected_ids, new_population_array, new_population_geo, new_OF 

def main():

    F = np.array([[1,2,3], [3,2,1], [1,1,1], [0,0,0], [4,4,4]])
    pareto_ranking = non_dom_sorting(F)
    print(pareto_ranking)
    
    F = np.array([[1,4],[3,3], [3,2], [8,1]])
    distance = crowding_distance(F)
    print(distance)
    
    #Testing nsga_rules
    all_designs = {}
    all_designs[1] = {"rem_budget":1}
    all_designs[2] = {"rem_budget":1}
    all_designs[3] = {"rem_budget":1}
    all_designs[4] = {"rem_budget":-41}
    all_designs[5] = {"rem_budget":1}
    all_designs[6] = {"rem_budget":1}
    all_designs[7] = {"rem_budget":1}
    all_designs[8] = {"rem_budget":1}

    
    OF = [[1,10],\
          [4,3],\
          [10,1],\
          [4.1,3.2],\
          [5,8],\
          [5.1,7.9],\
          [5.2,7.8],\
          [10,10]]
    OF = np.array(OF)
            
    pop_size = 5
    budget_threshold = 0
    pop_ids = list(all_designs.keys())
    
    population = []
    for k in all_designs.keys():
        population.append(3*[k])
    population = np.array(population)
    
    selected_ids, new_population, new_OF = nsga_rules(all_designs, pop_ids, population, OF, pop_size, budget_threshold)
    print("ids: ",selected_ids)
    print("new_population: ",new_population)
    print("new_OF: ", new_OF)
if __name__ == "__main__":
    main()
