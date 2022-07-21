# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 11:02:50 2022

@author: Owner
"""
import numpy as np
def GA_param():

    """
    [Description]

    OUTPUT:
        + route_GA_param -- 
        + sa_ga_param    --
    """

    route_GA_param = {}    
    route_GA_param["init_pop_size"] = 2000
    route_GA_param["max_eval"] = 100
    route_GA_param["pop_size"] = 64
    route_GA_param["random_gen"] =  np.random.default_rng()
    route_GA_param["p_crossover"] = 0.5
    route_GA_param["p_mutation"] = 0.2
    route_GA_param["num_objectives"] = 3
    route_GA_param["max_iter"] = 30
    route_GA_param['K_paths'] = 2 #to generate routes in univeral set of routes
    
    #Service area genetic algorithm
    sa_ga_param ={}
    sa_ga_param["p_num_generations"]  = 30
    sa_ga_param["p_crossover"] = 0.75
    sa_ga_param["p_mutation"] = 0.25
    sa_ga_param["p_population_size"] = 64
    sa_ga_param["index"] = list(range(sa_ga_param["p_population_size"]))
    sa_ga_param["random_gen"] =  np.random.default_rng()
    sa_ga_param['num_objectives'] = 3
    sa_ga_param["min_buffer_expand"] = 100
    sa_ga_param["max_buffer_expand"] = 25000
    
    mut_parameters = {}
    mut_parameters["min_num_mut"] = 1 # OLD: 2
    mut_parameters["max_num_mut"] = 4 # OLD: 5
    mut_parameters["p_expansion"] = 0.8 
    mut_parameters["min_selection_size"] = 0.025 #COULD CHANGE!!!!
    mut_parameters["max_selection_size"] = 0.05#COULD CHANGE!!!!
    mut_parameters["angle_scope"] = np.pi/6
    
    sa_ga_param["mut_parameters"] = mut_parameters
    
    #Parameters for grid operations
    grid_parameters = {}
    grid_parameters["box_size_buffer"] = 3960
    grid_parameters["grid_num_pts"] = 30
    grid_parameters["grid_point_2D_space"] =  np.zeros((grid_parameters["grid_num_pts"],2))
    grid_parameters["p_col_in"] = 0.5
    grid_parameters["window_length"] = 3
    
    sa_ga_param["grid_parameters"] = grid_parameters
    
    #Parameters for cross-over operation
    sa_ga_param["min_angle_slice"] = np.pi/6 
    sa_ga_param["angle_slice_span"] = np.pi/3
    
    
    
    return route_GA_param, sa_ga_param