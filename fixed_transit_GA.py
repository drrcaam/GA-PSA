# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 14:53:13 2022

@author: Owner
"""
import numpy as np
import random
import networkx as nx
import objective_functions as of
import moo_selector


def parent_crossover(x1, x2, p_crossover, num_set_route):
  """
  Cross over function.

  PARAMETERS:
    + x1            -- one dimensional numpy array; parent, vector of decision variables
    + x2            -- one dimensional numpy array; another parent, vector of decision variables
    + p_crossover   -- probability of cross over operation
    + num_set_route -- number of routes

  OUTPUT:
    + y1            -- one dimensional numpy array; first offspring, vector of decision variables
    + y2            -- one dimensional numpy array; second offspring, vector of decision variables

  """
  #Copy information to offspring
  y1 = x1.copy()
  y2 = x2.copy()

  #Find differences in route set
  not_in_1 = list(set(x2) - set(x1))
  not_in_2 = list(set(x1) - set(x2))
  num_diff = len(not_in_1)

  #Swap information (integers => routes):
  counter = 0
  for k in range(num_diff):
    if np.random.random() < p_crossover:
      
      idx = np.random.randint(num_set_route-1) #location to change
      y1[idx] = not_in_1[k]
      y2[idx] = not_in_2[k]
      counter += 1

  #Ensure that something is transferred
  if counter == 0: #if nothing change
    idx = np.random.randint(num_set_route-1) #location to change
    y1[idx] = not_in_1[0]
    y2[idx] = not_in_2[0]

  return y1, y2


def crossover(population, p_crossover, num_set_route, index_list, pop_size):
    """
    Coordinates crossover operation.
    
    PARAMETERS:
      + population    -- 2D array with route designs
      + p_crossover   -- probability of cross over operation
      + num_set_route -- number of routes
      + index_list    -- list with indices
      + pop_size      -- number of solutions in population
    
    OUTPUT:
      + offspring     -- 2D array with new route designs
    """
    
    #Shuffle the index used to assign partners for reproduction
    random.shuffle(index_list)
    
    #List with the children solutions
    offspring = np.zeros((pop_size, num_set_route))
    
    #For each parent:
    for k in range(0,pop_size-1, 2):
    
      #Perform crossover
      y1, y2 = parent_crossover(population[index_list[k],:], population[index_list[k+1],:], p_crossover, num_set_route)
      #Store:
      offspring[k,:] = y1.copy()
      offspring[k+1,:] = y2.copy()
    
    return offspring


def real_coded_mutation(population, p_mutation, all_routes, pop_size, num_set_route):
  """
  Functions assumes that if a mutation occurs, only one route is changed.

  PARAMETERS:
    + population    -- 2D array with route designs
    + p_mutation    -- probability of mutation
    + all_routes    -- list with all route ids
    + pop_size      -- size of population
    + num_set_route -- number of routes

  OUTPUTS: 
    + population    -- 2D array with route designs
  """
  for k in range(pop_size):
    #Mutations:
    if np.random.random() < p_mutation:
      not_in_set = list(set(all_routes) - set(population[k,:]))
      idx = np.random.randint(num_set_route)
      population[k,idx] = np.random.choice(not_in_set)

  return population


def gen_initial_pop(all_routes, num_set_route, init_pop_size, random_gen):
      """
      Generates an array of size (init_pop_size, num_set_route) that contains sets
      of randomly chosen integers that represent different route set designs. 
    
      PARAMETERS:
      + all_routes    -- integer list that represents all possible routes
      + num_set_route -- number of routes that must be in a route network
      + init_pop_size -- number of sets that will be generated in this initial population.
      + random_gen    -- numpy object that generates random numbers
    
      OUTPUT:
      + population    -- 2D array that contains candidate designs.
      """
    
      population = np.zeros((init_pop_size, num_set_route), dtype = "int64")
    
      #First design
      population[0,:] = np.sort(random_gen.choice(all_routes, size = (num_set_route,), replace = False))
    
      #Generating the rest of the route set desings:
      counter = 1
      while counter < init_pop_size:
    
        #Generate candidates.
        candidate = np.sort(random_gen.choice(all_routes, size = (num_set_route,), replace = False))
    
        #Check if it has been previously generated
        repeat = np.sum(np.sum(population[:counter,:] == candidate, axis = 1) == num_set_route)
        if repeat == 0:
          population[counter,:] = candidate
          counter += 1
    
      return population


def gen_universal_route_set(route_GA_param, transit_param):
    """
    Function generate the universe of possible routes.  It configures all 
    candidate routes for the current transit network using Dijkstra's shortest
    path and Yen's k-shortest path algorithms.

    Requires the user to define the minimum and maximum route lengths.

    Based on Wei Fan's Initial Candidate Route Set Generation Procedure (ICRSGP).
    
    PARAMETERS:
        + route_GA_param     -- dictionary with parameters that control GA 
        + transit_param      -- dictionary with parameters, including the network
                                and the length thresholds
                           
    OUTPUT:
        + route_set_universe -- dictionary with the set of all routes
        + all_routes         -- list with the ids of all the routes in the universe.
    """
    #Output:
    route_set_universe = {}
    
    #Get parameters:
    
    graph = transit_param["network_graph"]
    minimum_length = transit_param["minimum_route_length"]
    maximum_length = transit_param["maximum_route_length"]
    num_zones = transit_param["num_zones"]
    zone_2_node = transit_param["zone_2_node"]
    
    K_paths = route_GA_param['K_paths']
    
    route_id = 0
    
    #For each origin and destination:
    for origin in range(num_zones):
        for destination in range(num_zones):
            if origin != destination:
                #Create path generator
                origin_node = zone_2_node[origin]
                destination_node = zone_2_node[destination]
                path_generator = nx.shortest_simple_paths(graph, origin_node, destination_node, weight = "weight")
        
                for i in range(K_paths):# generate K paths
                    
                    try: 
                        path = next(path_generator) #get path
                        route_length = nx.path_weight(graph, path, weight = "weight")
                            
                        # Route constraints
                        if (route_length >= minimum_length) and (route_length <= maximum_length):
                            #Add route
                            route_set_universe[route_id] = {'path': path, 'length': route_length} 
                            route_id += 1
                    except Exception:
                        pass
           
    all_routes = list(range(route_id))
    return route_set_universe, all_routes

def design_evaluation(designs, route_set_universe, pop_size, max_eval, num_objectives, transit_param, all_designs, design_id ):
    """
    Coordinates the evaluation of the designs and the record keeping of solutions.
    
    PARAMETERS:
        + designs            -- 2D array that contains candidate designs
        + route_set_universe -- dictionary with the set of all routes
        + network_design     -- dictionary with the set of all routes
        + pop_size           -- number of designs; size of GA parent population
        + max_eval           -- maximum number solutions to consider 
        + num_objectives     -- number of design objectives
        + transit_param      -- dictionary with transit related parameters
        + all_designs        -- dictionary; archive of all evaluated solutions
        + design_id          -- counter of design; unique id
    
    OUTPUT:
        + OF                 -- 2D array of size (pop_size, num_objectives) with objective
                                function values of all designs
        + all_designs        -- as before, but updated with new designs
        + idx                -- list with the design ids of solutions evaluated in function
        + design_id          -- as before, but updated given new designs
        + designs            -- 2D array that contains candidate designs
    
    """
    OF = np.zeros((max_eval, num_objectives))
    ids = []
    counter = 0
    survivor_idx = []
    for k in range(pop_size):
        design = of.fixed_transit_OF(designs[k,:],route_set_universe, transit_param)
        if design != -1:
            all_designs[design_id] = design
            OF[counter, :] = design['OF']
            ids.append(design_id)
            design_id += 1
            counter +=1 
            survivor_idx.append(k)
            if counter >= max_eval:
                break


    return OF[:counter,:], all_designs, ids, design_id, designs[survivor_idx,:]
        
def routes_GA(route_GA_param, transit_param):

    """
    Genetic algorithm that searches for bus network designs

    PARAMETERS:
     + route_GA_param      -- dictionary with GA parameters
     + transit_param       -- dictionary with operations and demand parameters

    OUTPUT:
      + fixed_route_design -- dictionary with generated designs
      + design_id          -- number; counter of desings; ID
    """

    #Get parameters
    init_pop_size   = route_GA_param["init_pop_size"]
    pop_size        = route_GA_param["pop_size"]
    max_eval        = route_GA_param["max_eval"]
    random_gen      = route_GA_param["random_gen"]
    p_crossover     = route_GA_param["p_crossover"]
    p_mutation      = route_GA_param["p_mutation"]
    num_objectives  = route_GA_param["num_objectives"]
    max_iter        = route_GA_param["max_iter"]
    set_num_routes  = transit_param["set_num_routes"]
    budget_threshold = transit_param["budget_threshold"]

    
    #Generate universe of routes and their ids:
    route_set_universe, all_routes = gen_universal_route_set(route_GA_param, transit_param)    
    
    #Dictionary to store best designs for each maximum number of routes
    fixed_route_designs = {}
    
    #Counter for design ids
    design_id = 0
    
    #For each number of route 
    for num_routes in set_num_routes:
          
        #Create index list
        index_list = list(range(pop_size))


        #Generate initial population.

        population = gen_initial_pop(all_routes, num_routes, init_pop_size, random_gen)
      
        #Evaluate initial population
        all_designs = {}
        OF, all_designs, parent_ids, design_id, population = \
            design_evaluation(population, route_set_universe, init_pop_size, max_eval, num_objectives, transit_param, all_designs, design_id )
      
  
      
        #Perform NSGA-2 based selection of solutions
        parent_ids, population, OF =\
        moo_selector.nsga_rules(all_designs, parent_ids, population, OF, pop_size, budget_threshold)
          
        for k in range(max_iter):
            
            #Crossover.
            offspring = crossover(population, p_crossover, num_routes, index_list, pop_size)
        
            #Mutation.
            offspring = real_coded_mutation(offspring, p_mutation, all_routes, pop_size, num_routes)
        
            #Evaluate offspring.
            child_OF, all_designs, child_ids, design_id, offspring = \
                design_evaluation(offspring,route_set_universe, pop_size, max_eval, num_objectives, transit_param, all_designs, design_id )
             
            #Merge
            population = np.vstack((population, offspring))
            OF         = np.vstack((OF, child_OF)) 
            merged_ids = list(parent_ids) + child_ids

            #Perform selection based on nsga rules
            parent_ids, population, OF = \
                moo_selector.nsga_rules(all_designs, merged_ids, population, OF, pop_size, budget_threshold)
        

            
        fixed_route_designs[num_routes] =  {'all_designs':all_designs, "parent_ids": parent_ids}                       
    return fixed_route_designs, design_id