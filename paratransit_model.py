# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 10:30:01 2022
Paratransit demand model
@author: Owner
"""

import numpy as np
import matplotlib.path as mpltPath
import polygon_misc


def od_demand_model(sa_array, pzones_coor, od_matrix):
    """
    Function implements the paratransit demand estimation model.
    
    PARAMETERS:
        + sa_array           --  2D array with the sa polygon vertices
        + p_zones_coor       -- coordinates for the paratransit zones
        + od_matrix          -- DataFrame with trips for each origin-destination pair.
        
    OUTPUTS:
        + total_met_demand   -- accumulator with total trips that can be serviced
        + total_unmet_demand -- accumulator with total trips that cant be serviced
    """
    
    #Outputs
    total_met_demand = 0
    total_unmet_demand = 0
    
    #Determine the pzones that are withing the service area.
    path = mpltPath.Path(sa_array)
    inside_polygon = path.contains_points(pzones_coor)
    
    #Compute demand served and not served
    for idx, row in od_matrix.iterrows():
        if inside_polygon[row["origin"]] and inside_polygon[row["destination"]]:
            total_met_demand += row["trips"]
        else:
            total_unmet_demand += row["trips"]
            
    return total_met_demand, total_unmet_demand

def od_demand_model_pd(sa_array, pzones_data, od_matrix):
    """
    Function implements the paratransit demand estimation model.
    
    PARAMETERS:
        + sa_array           --  2D array with the sa polygon vertices
        + pzones_data        -- dataframe with coordinates for the paratransit zones
        + od_matrix          -- DataFrame with trips for each origin-destination pair.
        
    OUTPUTS:
        + total_met_demand   -- accumulator with total trips that can be serviced
        + total_unmet_demand -- accumulator with total trips that cant be serviced
    """
    
    
    #Determine the pzones that are withing the service area.
    path = mpltPath.Path(sa_array)
    inside_polygon = path.contains_points(pzones_data[['x', 'y']])
    inside_ids = pzones_data[inside_polygon]['id'].values
    
    #Compute demand served and not served
    mask = od_matrix['origin'].isin(inside_ids) & od_matrix['destination'].isin(inside_ids)
    total_met_demand =  od_matrix[mask]['trips'].sum()
    # total_unmet_demand =  #od_matrix[mask==False]['trips'].sum()  
    
    return total_met_demand #, total_unmet_demand #, inside_ids


def min_paratransit_service(transit_routes, transit_param):
    """
    Computes the paratransit objective function and service cost, as well
    as the equity objective for the integrated service. 
    
    PARAMETERS:
        + network_routes       -- list of list. Each list represents a different
                                  route. All the lists are the transit network.
        + transit_param        -- dictionary with parameters

    OUTPUT:
        + total_met_demand     --
        + min_paratransit_cost --
        + network_buffer       --
        + min_sa_geo           --
    
    """
    
    #Get parameters
    pzones_data = transit_param['pzones_data'] #data frame with pzone coordinates
    paratransit_od_matrix = transit_param['paratransit_od_matrix']
    paratransit_trip_cost = transit_param['paratransit_trip_cost']
    
    #Generate buffer.
    network_buffer, min_sa_array, min_sa_geo = polygon_misc.buffer_from_routes(transit_routes, transit_param)
    
    #Compute parantransit demand:
    total_met_demand = od_demand_model_pd(min_sa_array, pzones_data, paratransit_od_matrix)
        
    #Compute paratransit cost
    min_paratransit_cost = paratransit_trip_cost * total_met_demand
    
    return total_met_demand, min_paratransit_cost, network_buffer, min_sa_geo


def eval_paratransit_service(sa_array, transit_param):
    """
    Computes the paratransit objective function and service cost, as well
    as the equity objective for the integrated service. 
    
    PARAMETERS:
        + sa_array         -- 2D array representing service area
        + transit_param    -- dictionary with parameters

    OUTPUT:
        + total_met_demand --
        + paratransit_cost --
    
    """
    
    #Get parameters
    pzones_data = transit_param['pzones_data'] #data frame with pzone coordinates
    paratransit_od_matrix = transit_param['paratransit_od_matrix']
    paratransit_trip_cost = transit_param['paratransit_trip_cost']
    
    #Compute parantransit demand:
    total_met_demand = od_demand_model_pd(sa_array, pzones_data, paratransit_od_matrix)
        
    #Compute paratransit cost
    paratransit_cost = paratransit_trip_cost * total_met_demand
    
    return total_met_demand, paratransit_cost