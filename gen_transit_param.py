# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 10:45:00 2022
Module generates the parameters related to the transit system.
@author: Owner
"""
import pandas as pd
import numpy as np
import gen_network_v2 as gen_network
import geopandas as gpd
import pickle
import pathlib


def gen():
    """
    [Description]

    OUTPUT:
    + transit_param -- 
    """
    
    transit_param = {}
    
    path_ = str(pathlib.Path(__file__).parent.resolve()) + '\\'
    
    ##############################################
    ########Parameters for fixed route problem
    ##############################################
    
    #Network and zonal data:
    transit_param["set_num_routes"] = [18]
    
    filename = path_ + "road_network_v3.xlsx"
    network_data = pd.read_excel(filename)
   
    filename = path_ + "demographic_data.xlsx"
    zone_data   = pd.read_excel(filename , sheet_name = "zone_centroids")
    
   #Create network
    bus_speed = 3080 #ft/minute;  the units in the projection used are in feet
    road_network, road_node_list, road_edge = gen_network.create_graph_shp(network_data, np.inf, bus_speed)
    transit_param["network_graph"] = road_network
    transit_param["road_edge"] = road_edge
    transit_param["max_node_id"] = road_node_list[-1][2]
    transit_param["transfer_time"] = 10 #minutes
    transit_param["transfer_mark"] = 0.001
    transit_param["minimum_route_length"] = 30000/bus_speed #ft
    transit_param["maximum_route_length"] = 110000/bus_speed #ft
    
    
    transit_param["universe_exists"] = True

    
    #Creating dictionary that contains the coordinates of all nodes
    #The coordinate information is contain in the road_node_list, which contains
    #a list for each vertex=node. Each vertex list has three numbers: 
        #latitude, longitude and the node id.
    #In the geometric operations performed in the models, the 
    #x is longitue and the y is latitude (the nomenclature used by QGIS). This 
    #is why the order of the information is adjusted in the vertex_coordinates
    vertex_coordinates = {}
    for node in road_node_list:
        node_id = node[2]
        node_x = node[1]
        node_y = node[0]
        vertex_coordinates[node_id] = (node_x,node_y)
    transit_param['vertex_coordinates'] = vertex_coordinates
    
    #Create map between network nodes and zones.
    zone_2_node = {}
    rode_node_array = np.array(road_node_list)
    for idx, row in zone_data.iterrows():
        loc = np.argmin((row["x"] - rode_node_array[:,1])**2 + (row["y"] - rode_node_array[:,0])**2)
        zone_2_node[int(row["ID"])] = int(rode_node_array[loc, 2])

    transit_param["zone_2_node"] = zone_2_node
         
    #Demand:
    transit_param["num_zones"] = len(zone_2_node)
    filename = path_ + "trip_table.csv"
    transit_param["OD_flows"] = pd.read_csv(filename)
    transit_param["num_ODs"] = len(transit_param["OD_flows"])
    
    #Behavior:
    transit_param["K_paths"]       = 2
    transit_param["B_time"]        = -0.1 
    
    #Operations
    transit_param["load_factor"]   = 1
    transit_param["budget"]        = 45000000 - 12000000
    transit_param["budget_threshold"] = 0
    transit_param["vehicle_cap"]   = 100
    num_days = 365-52 #operates monday to saturday 
    oper_hours = 8
    oper_exp_per_veh_rev_hour = 130
    transit_param["unit_cost_ft"]  = num_days*oper_hours*oper_exp_per_veh_rev_hour
    transit_param["h_min"]  = 5 #minutes
    transit_param["h_max"]  = 60 #minutes

    #Objectives
    transit_param['weight_unmet_demand'] = 10
    transit_param['weight_user_costs'] = 1

    ##############################################
    ########Parameters for paratransit model
    ##############################################
    
    #Geographic
    transit_param
    transit_param['min_buffer_size'] = 3960 #ft, 3/4 of a mile in feet
    transit_param['project_crs'] = "EPSG:3991"
    transit_param['vertex_orientation'] = 1.0
    # transit_param['centroid_coordinates'] = {}

    
    #Demand:
    filename = path_ + "demographic_data.xlsx"
    pzones_data = pd.read_excel(filename, sheet_name = 'pzones_data') #data frame with pzone coordinates
    pzones_data['GEOID'] = pzones_data['GEOID'].astype(str)
    transit_param['pzones_data'] = pzones_data
        
    filename = path_ + "trip_table_paratransit.csv"
    paratransit_od_matrix = pd.read_csv(filename)
    paratransit_od_matrix["trips"] = 2*paratransit_od_matrix["trips"]
    transit_param['paratransit_od_matrix'] = paratransit_od_matrix
    transit_param['total_paratransit_demand'] = np.sum(paratransit_od_matrix["trips"])
    
    filename = path_ + "zonal_system_v1_3991.shp"
    transit_param["zonal_shp"] = gpd.read_file(filename)
    
    #Operations:
    transit_param['paratransit_trip_cost'] = 45
    
    ##############################################
    ########Parameters for equity function
    ##############################################
    
    transit_param['atkinson_epsilon'] = 1 - 0.75 
    transit_param["pop_group"] =  "target_pop"
    transit_param["spatial_group_column"] = "COUNTY"
    transit_param["county_list"] = list(np.unique(transit_param['pzones_data'][transit_param["spatial_group_column"]]))
    transit_param["pop_container"] = np.zeros((len(transit_param["county_list"]),))
    
    
    #Computing the number of people that live in the regions of interest
    county_list = transit_param["county_list"]
    pop_container = np.zeros((len(county_list),))
    agg_zonal_pop = pzones_data.groupby(transit_param["spatial_group_column"])[transit_param["pop_group"]].sum() 
    
    for k in range(0,len(county_list)):
        if county_list[k] in agg_zonal_pop.index:
            pop_container[k] = agg_zonal_pop.loc[county_list[k]]
        else: 
            pop_container[k] = 0
  
    pop_container = np.append(pop_container,pop_container)
    transit_param['total_pop_group'] = pop_container
    
    #Original coverage
    filehandler = open("original_coverage.obj", 'rb') 
    original_coverage = pickle.load(filehandler)
    filehandler.close()
    transit_param['original_coverage'] = original_coverage
    
    return transit_param
