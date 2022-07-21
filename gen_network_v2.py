# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 10:16:39 2020
Network generation
@author: Owner
"""


import networkx as nx
import numpy as np
import math 

def check_unique_node(node_coordinate, node_list, node_id_counter):
    """
    This function determines if a node is unique based on its coordinates.
    If it is unique it is given a new unique ID, otherwise, it is given it
    is given its existing ID.
    
    INPUTS:
        + node_coordinate -- (lat, long) coordinate for node.
        + node_list       -- list with unique node id and their coordinates
        + node_id_counter -- node counter
    
    OUTPUT: 
        + node_id         --  id for the node_coordinates
        + node_list       -- list with unique node id and their coordinates
        + node_id_counter -- node counter
    """
    num_nodes = len(node_list)
    new = True #flag that stats if the node is new. initially set to True
    #If this is the first node considered
    if num_nodes == 0:
        #The id of the fist node is the number in node_id_counter (which is 0)
        node_id = node_id_counter
        node_list.append([node_coordinate[0], node_coordinate[1], node_id])
        node_id_counter += 1 #update counter
    else:
        #for each node in node list
        for k in range(0, num_nodes):
            #check if the coordinate of intesest is equal to an existing one
            if (node_coordinate[0] == node_list[k][0]) and (node_coordinate[1] == node_list[k][1]):
                node_id = node_list[k][2] #the node id of the node is the id found in the list
                new = False
                break
    
        #If the node is new, assign the new node id, and update the node_list and the node counter.
        if new:
            node_id = node_id_counter
            node_list.append([node_coordinate[0], node_coordinate[1], node_id])
            node_id_counter += 1 #update counter
            
    return node_id, node_list, node_id_counter
    

#Create network
def create_graph_shp(network, max_segment_len, bus_speed):
    """
    This function creates an undirected graph for the networkx using
    link and node information extracted from a line shapefile.
    
    #INPUTS:
        + network         -- array with the link ids, the node ids and the node coordinates
        + max_segment_len -- the maximum length of a segment which in turn controls
                             the number of vertices)
        + bus_speed       -- average speed of buses

    OUTPUT:
        + G               -- networkx graph
        + node_lists      -- list with node information
        + edge_dict       -- dictionary with edge information
    """
    
    #Create empty graph
    G = nx.Graph()
    #Define the variable to track node list and unique node ids
    node_list = []
    edge_dict = {}
    node_id_counter = 0
    #Get the ids for each link
    all_link_id = np.unique(network["id"])
    
    #For each link
    for link_id in all_link_id:
        #Get data for the link
        mask = network["id"] == link_id 
        link_data = network[mask]
        sorted_data = link_data.sort_values(by='vertex_index')
        num_nodes = link_data.shape[0]
        #For each node pair, compute distance
        vertex_1 = (sorted_data.iloc[0]['lat'], sorted_data.iloc[0]['long']) #first finding the ID of the first node
        node_id_1, node_list, node_id_counter = check_unique_node(vertex_1, node_list, node_id_counter)
        
        for k in range(0,num_nodes - 1):
            #Set ID for the arrival vertex
            vertex_2  = (sorted_data.iloc[k+1]['lat'], sorted_data.iloc[k+1]['long'])
           
            #find distance between vertices
            edge_length =  ((vertex_1[0]-vertex_2[0])**2 + (vertex_1[1]-vertex_2[1])**2)**0.5#geodesic( vertex_1, vertex_2).m 
            edge_length = edge_length / bus_speed
            
            #Add edge to graph
            node_id_2, node_list, node_id_counter = check_unique_node(vertex_2, node_list, node_id_counter)
            G.add_edge(node_id_1, node_id_2, weight = edge_length)
            edge_dict[(node_id_1, node_id_2)] = edge_length
            edge_dict[(node_id_2, node_id_1)] = edge_length
                
            #Update the vertex 1
            vertex_1 = vertex_2
            node_id_1 = node_id_2
  
    return G, node_list, edge_dict


#Create network
def create_graph_shp_reduced(network):
    """
    This function creates an undirected graph for the networkx using
    link and node information extracted from a line shapefile.
    
    #INPUTS:
        + network         -- array with the link ids, the node ids and the node coordinates
        + max_segment_len -- the maximum length of a segment which in turn controls
                             the number of vertices)

    OUTPUT: 
        + G               -- networkx graph
        + node_list      -- list with node information
        + edge_dict       -- dictionary with edge information
    """
    
    #Create empty graph
    G = nx.Graph()
    #Define the variable to track node list and unique node ids
    node_list = []
    edge_dict = {}
    node_id_counter = 0
    #Get the ids for each link
    all_link_id = np.unique(network["id"])
    #For each link
    for link_id in all_link_id:
        
        total_edge_length = 0
        
        #Get data for the link
        mask = network["id"] == link_id 
        link_data = network[mask]
        sorted_data = link_data.sort_values(by='vertex_index')
        num_nodes = link_data.shape[0]
        
        #For each node pair, compute distance
        vertex_1 = (sorted_data.iloc[0]['lat'], sorted_data.iloc[0]['long']) #first finding the ID of the first node
        node_id_1, node_list, node_id_counter = check_unique_node(vertex_1, node_list, node_id_counter)
        
        for k in range(0,num_nodes - 1):
            #Set ID for the arrival vertex
            vertex_2  = (sorted_data.iloc[k+1]['lat'], sorted_data.iloc[k+1]['long'])
           
            #find distance between vertices
            edge_length =   ((vertex_1[0]-vertex_2[0])**2 + (vertex_1[1]-vertex_2[1])**2)**0.5#geodesic( vertex_1, vertex_2).m 
            total_edge_length += edge_length

            
            #If in last pair of nodes
            if k == num_nodes - 2:
                #Add edge to graph
                node_id_2, node_list, node_id_counter = check_unique_node(vertex_2, node_list, node_id_counter)
            
                G.add_edge(node_id_1, node_id_2, weight = total_edge_length)
                edge_dict[(node_id_1, node_id_2)] = total_edge_length
                edge_dict[(node_id_2, node_id_1)] = total_edge_length

            #Update the vertex 1
            vertex_1 = vertex_2

  
    return G, node_list, edge_dict



def get_closest_vertex(coordinate, node_list):
    """
    Gets the node closes to a coordiante (at some point this function computed the geodesic distance... but it was too time consuming perhaps)

    PARAMETERS:
        + coordinate -- list/tuple with coodinate
        + node_list  -- vertex list (each element has three points of info; x,y,id)

    OUTPUT:
        + position   -- id of vertex
    """
    min_dist = math.inf
    position = -1
    
    # For each of the other points
    for vertex in node_list:
        temp   =  ((vertex[0]-coordinate[0])**2 + (vertex[1]-coordinate[1])**2)**0.5 #get distance to end coordinate

        if temp < min_dist: #if temp is less than previous closest cummulative distance
            min_dist = temp #update the containers 
            position = vertex[2]
    
    return position
    
    
def network_distance(start, end, node_list, network):
    """This function finds the length of the shortest path between the coordinates start and end.
    
    INPUTS:
        + start       -- starting coordinate
        + end         -- ending coordinate point
        + node_list   -- list with node information (each element has three points of info; x,y,id)
        + network     -- networkx graph of the scooter network
    
    OUTPUT:
        + path_length -- shortest path distance between start and end.
    
    """
    #Find the network nodes closest to the coordinates.
    vertex_start = get_closest_vertex(start, node_list)
    vertex_end   = get_closest_vertex(end, node_list)
    #Get shortest distance:
    path_length = nx.shortest_path_length(network,vertex_start,vertex_end,weight = 'weight')
    return path_length
    

def get_closest_euclid(coordinate, node_list):
    """
    Gets the node closes to a coordiante

    PARAMETERS:
        + coordinate -- list/tuple with coodinate
        + node_list  -- vertex list (each element has three points of info; x,y,id)

    OUTPUT:
        + position   -- id of vertex
    """
    #Initialize
    min_dist = math.inf
    position = -1
    
    # For each of the other points
    for vertex in node_list:
        temp   = ((coordinate[0] -vertex[0])**2 + (coordinate[1] - vertex[1])**2 )**0.5

        if temp < min_dist: #if temp is less than previous closest cummulative distance
            min_dist = temp #update the containers 
            position = vertex[2]
    
    return position