#Base libraries
import numpy as np
import random
from shapely.ops import orient
import polygon_misc
from grid_mutation import grid_mutation_coordinator_v2
from crossover import sa_crossover
import objective_functions as of
import moo_selector
import pickle



def initial_pop(N, base_buffer, min_buffer_size, max_buffer_size):
  """
  Based on the fixed route network, this function generates an initial set of 
  service area boundaries using buffer functions. 

  #PARAMETERS:
    + N               -- size of the population of design solutions
    + base_buffer     --
    + network         -- object that represents sets of fixed transit routes
    + min_buffer_size -- minimum size of service area buffer 
    + max_buffer_size -- maximum size of service area buffer 

  #OUTPUT:
    + population      -- list of 2D numpy arrays; each array contains coordinates of a service area (i.e., polygon)
    + population_array --
    + population_geo   --

  """
  #Storage of arrays
  population_array = []
  population_geo = []
  all_sizes =  np.linspace(min_buffer_size, max_buffer_size, N)

  for size in all_sizes:
    #Generate buffer
    new_buffer = base_buffer.buffer(size)
    
    # Combining all buffers to form one single buffer polygon around routes
    new_buffer = new_buffer.unary_union
    
    # Orienting vertices to counter_clockwise
    new_buffer = orient(new_buffer, sign=1.0)
    
    # Extracting vertices from Shapely polygon
    buffer_vertices = new_buffer.exterior.coords
    
    # Buffer vertices as a 2D numpy array
    buffer = np.array(buffer_vertices)

    #Add to population list
    population_array.append(buffer)
    population_geo.append(polygon_misc.polygon_from_vertices(buffer))
    
  return population_array, population_geo

def design_evaluation(service_area, service_area_geo, num_objectives, \
                      transit_param, ft_design, id_ft_design, all_designs, design_id ):
    
    """
    Coordinates the evaluation of the designs and the record keeping of solutions.
    
    PARAMETERS:
        + service_area     -- list with 2D arrays representing service area
        + service_area_geo -- list with shapely polygons representing service area
        + num_objectives   -- number of design objectives
        + transit_param    -- dictionary with transit related parameters
        + ft_design        -- components underlying fixed transit (ft) route network
        + id_ft_design     -- id of the ft design
        + all_designs      -- dictionary; archive of all evaluated solutions
        + design_id        -- counter of design; unique id
    
    OUTPUT:
        + OF               -- 2D array of size (pop_size, num_objectives) with objective
                              function values of all designs
        + all_designs      -- as before, but updated with new designs
        + idx              -- list with the design ids of solutions evaluated in function
        + design_id        -- as before, but updated given new designs
    
    """
    pop_size = len(service_area)
    OF = np.zeros((pop_size, num_objectives))
    ids = []
    
    for k in range(pop_size):
        design = of.paratransit_OF(service_area[k], service_area_geo[k], ft_design, id_ft_design, transit_param)
        all_designs[design_id] = design
        OF[k, :] = design['OF']
        ids.append(design_id)
        design_id += 1

    return OF, all_designs, ids, design_id

def offspring_feasibility_check(offspring_array, min_buffer, project_crs):
    """
    Checks if offsprings are valid and, if so, collects the geoseries 
    representation.
    
    PARAMETERS
        + offspring_array -- list with 2D arrays representing SA
        + min_buffer      -- shapely polygon representing minimum coverage
        + project_crs     --
        
    OUTPUT
        + valid_sa_array  -- list with 2D arrays representing SA
        + valid_sa_geo    -- list with geoseries representing SA
    """
    #Containers
    valid_sa_array = []
    valid_sa_geo = []
    
    for polygon_vertices in offspring_array:
        valid, sa_geo = polygon_misc.service_area_geo_validity(polygon_vertices, min_buffer, project_crs)
        if valid:
            valid_sa_array.append(polygon_vertices)
            valid_sa_geo.append(sa_geo)

    return valid_sa_array, valid_sa_geo


def sa_ga_post_gen(ft_design ,id_ft_design, design_id,  transit_param, ga_param, all_upgraded_designs):
    """
    [Description]

    PARAMETERS:
        + ft_design            --
        + id_ft_design         --
        + design_id            --
        + transit_param        --
        + ga_param             --
        + all_upgraded_designs -- 

    OUTPUT:
        + all_upgraded_designs --
        + design_id            --
    """
    
    #Get parameters
    num_generations = ga_param["p_num_generations"]
    p_crossover     = ga_param["p_crossover"]
    p_mutation      = ga_param["p_mutation"]
    pop_size        = ga_param["p_population_size"]
    index           = ga_param["index"]
    grid_parameters = ga_param["grid_parameters"]
    mut_parameters  = ga_param["mut_parameters"]
    min_buffer_expand = ga_param["min_buffer_expand"]
    max_buffer_expand = ga_param["max_buffer_expand"]
    num_objectives  = ga_param["num_objectives"]
    
    #Read solution feautes 
    min_buffer = ft_design['sa']
    
    #Set initial population values
    population_array, population_geo = initial_pop(pop_size, min_buffer, min_buffer_expand, max_buffer_expand)
    
    #Evaluate initial population
    OF, all_upgraded_designs, parent_ids, design_id = design_evaluation(population_array, population_geo, num_objectives, \
                          transit_param, ft_design, id_ft_design, all_upgraded_designs, design_id )
    
    #Commence evolution:
    for g in range(num_generations): 
        print("sa generation :", g)
        
        #Shuffle the index used to assign partners for reproduction
        random.shuffle(index)
      
        #List with the children solutions
        offspring_array = []
      
        #For each parent:
        for k in range(0,pop_size):
            
            #Cross over operation
            if np.random.random() < p_crossover:
        
              #Get parents   
              parent_1 = population_array[index[k]]
              parent_2 = population_array[index[k-1]]
              
              #Perform cross over 
              child_1, child_2 = sa_crossover(parent_1, parent_2, ga_param)
        
              #Consider mutation operation
              if np.random.random() < p_mutation:
                  child_1 = grid_mutation_coordinator_v2(child_1, mut_parameters, grid_parameters)
                  # child_2 = grid_mutation_coordinator_v2(child_2, mut_parameters, grid_parameters)
        
              #Add to offspring list
              offspring_array.append(child_1)
              # offspring_array.append(child_2)    
        
            else:
              #Mutation operation (no crossover)
              parent_1 = population_array[index[k]]
              child_1  = grid_mutation_coordinator_v2(parent_1, mut_parameters, grid_parameters)
              
              #Add to offspring list
              offspring_array.append(child_1)
          
        #Feasibility check: 
        offspring_array, offspring_geo = \
            offspring_feasibility_check(offspring_array, min_buffer, transit_param['project_crs'])
        
        #Evaluate the fitness of offspring. 
        child_OF, all_upgraded_designs, child_ids, design_id = design_evaluation(offspring_array, offspring_geo, num_objectives, \
                              transit_param, ft_design, id_ft_design, all_upgraded_designs, design_id )
      
        #Merge
        population_array = population_array + offspring_array
        population_geo = population_geo + offspring_geo
        OF         = np.vstack((OF, child_OF)) 
        merged_ids = list(parent_ids) + child_ids
    
        #Perform selection based on nsga rules
        parent_ids, population_array, population_geo, OF  = \
            moo_selector.nsga_rules_sa(all_upgraded_designs, merged_ids, \
           population_array, population_geo, OF, pop_size, transit_param['budget_threshold'])


    return all_upgraded_designs, design_id
    
