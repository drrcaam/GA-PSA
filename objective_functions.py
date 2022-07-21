# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 10:45:04 2022

@author: Owner
"""
import transit_assignment as ta
import paratransit_model as pm
import equity_objective as eq
import geopandas as gpd


def fixed_transit_OF(network_design, route_set_universe, transit_param):
    """
    This function compute the objective function values necessary to perform
    the design evaluation at the fixed-route design stage. 
    
    In the fixed route design stage, only the minimum paratransit service level
    is determined.
    
    #PARAMETERS:
        + network_design     -- 1D array with the ID of the routes that constitute
                                the transit network
        + route_set_universe -- dictionary with the set of all routes
        + transit_param      -- dictionary with parameters
        
    #OUTPUT: 
        Dictionary that contains the:
            -list with objective function values (fixed route, paratransit, atkison),
            -2D array with minimum paratransit service area
            -remaining budget
            -list with number of population groups covered by fixed transit
        
    
    """
    #Get parameters:
    budget = transit_param['budget']
    project_crs = transit_param['project_crs']
    total_paratransit_demand = transit_param['total_paratransit_demand']
    weight_unmet_demand = transit_param['weight_unmet_demand']
    weight_user_costs   = transit_param['weight_user_costs']
    
    #Generate list with the paths that constitue the transit network
    transit_routes = []
    for route_id in network_design:
        transit_routes.append(route_set_universe[route_id]['path'])
    
    try:
        #Compute trips, cost, and the service area buffer for minimum paratransit
        #service. 
        p_met_demand, min_paratransit_cost, network_buffer, min_sa_geo =\
            pm.min_paratransit_service(transit_routes, transit_param)
        p_objective = total_paratransit_demand - p_met_demand
        
        
        #Perform fixed route transit assignment:
        ft_unmet_demand, fixed_transit_cost, user_costs = ta.assignment(transit_routes, transit_param)
    
        ft_objective = weight_unmet_demand*ft_unmet_demand + weight_user_costs*user_costs
        
    
        #Compute remaining budget
        rem_budget = budget - min_paratransit_cost - fixed_transit_cost
        
        #Compute equity objective
        network_buffer_gpd = gpd.GeoDataFrame(geometry=network_buffer,crs=project_crs)
        sa_gpd =  gpd.GeoDataFrame(geometry=min_sa_geo,crs=project_crs)
        atkinson, bus_pop_group = eq.spatial_proximity_geo(network_buffer_gpd,sa_gpd, transit_param)
        
        return {"OF": [ft_objective, p_objective, atkinson ], \
                "ft_coverage": network_buffer, \
                "sa": min_sa_geo,\
                "rem_budget":rem_budget, \
                "bus_pop_group": bus_pop_group,\
                "transit_routes" : transit_routes,\
                "route_ids": network_design,\
                "p_transit_cost":  min_paratransit_cost}
            
    except:
        return -1 #route network is disconnect or has some problem
    
def paratransit_OF(service_area, service_area_geo, ft_design, id_ft_design, transit_param):
    """
    This function compute the objective function values necessary to perform
    the design evaluation when service area improvements are considered. 
    
    #PARAMETERS:
        + service_area     -- 2D array representing service area
        + service_area_geo -- geoseries object representing service area
        + ft_design        -- dictionary with components of underlying fixed route network
        + id_ft_design     -- id of ft_design
        + transit_param    -- dictionary with parameters
        
    #OUTPUT: 
        Dictionary that contains the:
            -list with objective function values (fixed route, paratransit, atkison),
            -geopandas object with minimum paratransit service area
            -remaining budget
            -list with number of population groups covered by fixed transit
        
    
    """
    #Get parameters
    total_paratransit_demand = transit_param['total_paratransit_demand']
    project_crs = transit_param['project_crs']
    
    #Get data from fixed transit design
    rem_budget = ft_design["rem_budget"]
    min_paratransit_cost = ft_design["p_transit_cost"]
    bus_pop_group = ft_design['bus_pop_group']
    
    #Evaluate paratransit design
    p_met_demand, paratransit_cost = \
        pm.eval_paratransit_service(service_area, transit_param)
        
    p_objective = total_paratransit_demand - p_met_demand
    
    #Compute remaining budget
    rem_budget = rem_budget + min_paratransit_cost - paratransit_cost
    
    #Compute equity objective
    service_area_gpd = gpd.GeoDataFrame(geometry=service_area_geo,crs=project_crs)
    atkinson, paratransit_pop_group  = \
        eq.spatial_proximity_geo_update(service_area_gpd, bus_pop_group, transit_param)
    
    return {"OF": [ft_design["OF"][0], p_objective, atkinson ], \
            "sa": service_area_geo, \
                
            "rem_budget":rem_budget, \
            "paratransit_pop_group": paratransit_pop_group,\
            "transit_routes" : ft_design["transit_routes"],\
            "route_ids": ft_design["route_ids"],\
            "p_transit_cost":  paratransit_cost,\
            "id_ft_design":id_ft_design}
        
    
