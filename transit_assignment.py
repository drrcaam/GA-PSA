# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 16:07:28 2022

@author: Daniel
"""


#Modules
import numpy as np
import networkx as nx

#Functions
def transit_network_builder(all_routes, road_edge, transfer_mark, max_node_id, K_paths):
    """
    This function creates a networkx graph given a set of routes. Edges
    representing transfers between routes as added.
    
    PARAMETERS:
        + all_routes    -- list with routes. a route is a list of integers. the integers
                           are road node ids that constitute the route.
        + road_edge     -- dictionary that contains all the edges of the road network.
                           keys are tuples with the nodes that define the edge and the
                           values are the edge weight.
        + transfer_mark -- flag given to transfer links. Small number. More cost added when computing OD paths.
        + max_node_id   -- maximum id for node in road network 
        + K-paths       -- number of paths considered by people when deciding travel routes.
    
    OUTPUT:
        + G               -- networkx graph representing the transit network. 
        + link_map        -- 
        + edge_dict       -- dictionary that serves as the link_map in other functions.
                             link_map maps the ids of the edges in G to their original
                             ids which appear on the lists that define each route.
        + all_cost        -- list that stores the total weight of each transit route 
        + link_dictionary -- dictionary that links the transit edges with their original
                            node ids with the position of the link in the link_array
        + link_array      -- arrays of size (K_paths, num of transit links). Does not saves transfer links 
        + link_flows      -- array of size (num of transit links,). Will contain flows on transit links. 
        
    """
    
    G = nx.Graph() # network

    link_map = {} # dictionary of use to map to network links
    all_cost = np.zeros((len(all_routes),)) #array with route costs
    link_dictionary = {} #to map from original link to position in link arrays

    #Objects used to keep track of rode nodes that repeat themselves in the route
    visited_nodes = []
    repeat_map = {}
    
    #Counters
    link_counter = 0 #link counter

    #For each route.
    for route_number in range(len(all_routes)):
        route = all_routes[route_number].copy()
        total_cost = 0

        #Dictionary to keep track fre
        transit_2_road = {}
        for node in route:
            transit_2_road[node] = node
        
        new_nodes = [] #list of nodes unique to this route
        
        for k in range(len(route)-1): 
            v1 = route[k] #first node
            v2 = route[k+1] #second node
            
            #Links with original node ids
            link1 = (route_number, transit_2_road[v1], transit_2_road[v2])
            link2 = (route_number, transit_2_road[v2],  transit_2_road[v1])

            link_dictionary[link1] = link_counter #add segment to dictionary
            link_counter += 1
            link_dictionary[link2] = link_counter #add segment to dictionary
            link_counter += 1

            cost = road_edge[(transit_2_road[v1],transit_2_road[v2])]
            total_cost += cost

            #Check if first node has been visited before          
            if v1 in visited_nodes:
                #Create new unique ID for node:
                temp_id = max_node_id + 1
                
                #Add tranfer links and update repeat map
                G.add_edge(v1, temp_id, weight = transfer_mark) #max_node_id + 1 is the new node id for the node in this route
                
                if v1 in repeat_map:
                    # for node in repeat_map[v1]:
                    #     G.add_edge(node, temp_id, weight = transfer_mark) 
                    
                    #Update
                    repeat_map[v1].append(temp_id)
                else:
                    repeat_map[v1] = [temp_id] #creating list for node
                    
                #Update node id in current route
                route[k] = temp_id
                
                #Update mapping to road network and id
                transit_2_road[temp_id] = v1
                v1 = temp_id
                max_node_id += 1
                
                
            else:
                new_nodes.append(v1)
            
            #Check if first second has been visited before
            if v2 in visited_nodes:
                #Create new unique ID for node:
                temp_id = max_node_id + 1
                
                #Add tranfer links and update repeat map
                G.add_edge(v2, temp_id, weight = transfer_mark) #max_node_id + 1 is the new node id for the node in this route
                
                if v2 in repeat_map:
                    # for node in repeat_map[v2]:
                    #     G.add_edge(node, temp_id, weight = transfer_mark) 
                    
                    #Update
                    repeat_map[v2].append(temp_id)
                else:
                    repeat_map[v2] = [temp_id] #creating list for node
                    
                #Update node id in current route
                route[k+1] = temp_id
                
                #Update id
                transit_2_road[temp_id] = v2
                v2 = temp_id
                max_node_id += 1
                
            else:
                new_nodes.append(v2)            
            
            #Add edges 
            G.add_edge(v1, v2, weight = cost)
            link_map[(v1, v2)] = link1
            link_map[(v2, v1)] = link2
        
        #Ending route analysis
        all_cost[route_number] = total_cost #had parenthesis around it
        
        visited_nodes = visited_nodes + new_nodes

    #Create arrays that will store flows on the transit network
    link_array = np.zeros((K_paths, link_counter))
    link_flows = np.zeros((link_counter,))
  
    return G, link_map, all_cost, link_dictionary, link_array, link_flows

def gen_route_maps(transit_routes, K_paths):
  """
  Outputs link_dictionary and link_array objects that are used to aggregate 
  flows in each route. 

    PARAMETERS:
      + transit_routes  --
      + K_paths         --

    OUTPUT:
      + link_dictionary --
      + link_array      --
      + link-flows      --
  """

  link_dictionary = {}
  counter = 0
  #For each transit route:
  for k in range(len(transit_routes)):
    route = transit_routes[k] #get route
    #For each segment in the route. 
    for j in range(len(route)-1):
      link_dictionary[(k, route[j], route[j+1])] = counter #add segment to dictionary
      counter += 1
      link_dictionary[(k, route[j+1], route[j])] = counter #add segment to dictionary
      counter += 1
      
  link_array = np.zeros((K_paths, counter-1))
  link_flows = np.zeros((counter-1,))

  return link_dictionary, link_array, link_flows
    
def route_assignment(OD, transit_network, OD_demand, K_paths, B_time, link_dictionary, link_array, link_map, transfer_time, transfer_mark):
  """
  Route assignment model. Demand split between K shortest path. 
  Probability of selecting a path is determined using logit model. 

  PARAMETERS: 
    + OD              -- tuple that contains ID of origin and destination zone.
    + transit_network -- networkx object that represents the transit network.
    + OD_demand       -- number of trips from origin zone to destination zone.
    + K_paths         -- number of paths considered in assignment.
    + B_time          -- logit model parameter. Utilty of every path is V = B_time*path_cost. 
    + link_dictionary -- dictionary that maps network link to its location in link array.
    + link_array      -- array with dimensions (K_paths, #links in road network).
    + link_map        --
    + transfer_time   -- number that represents the transfer times between routes
    + transfer_mark   -- number that serves as a flag of a transfer link. 
                         it is edge weight that is deleted from the path cost. 

  OUTPUT:
    + 1-dimensional array with the OD_demand distributed across road network links 
    + user_cost       --

  """
  
  origin = OD[0]
  destination = OD[1]
    
  #Find the K shortest paths
  path_generator = nx.shortest_simple_paths(transit_network, origin, destination, weight = "weight")
  path_cost = []
  num_paths = 0
  #for each path
  for i in range(K_paths):
      
      try: 
        path = next(path_generator)
        path_cost.append(nx.path_weight(transit_network, path, weight = "weight"))
        num_paths += 1

        #Unravel path.
        num_links = len(path)-1
        for j in range(0, num_links):
          #networkx path has "augmented" node ids for the transit network. 
          #here we are finding the original node ids for the original transit segment 

          if (path[j], path[j+1]) in link_map:
            original_link = link_map[(path[j], path[j+1])]  

            #With the original link, we use link_dictionary to get position of the 
            #transit segment on the link_array
            link_loc = link_dictionary[original_link] 

            #Add 1 to mean that the segment will receive flow
            link_array[i, link_loc] = 1

          elif (j != 0) and (j != num_links-1):
            path_cost[i] = path_cost[i] + transfer_time 

          else:
            path_cost[i] = path_cost[i] - transfer_mark

      except:
        break

  #Compute route probabilities
  
  path_cost = np.array(path_cost)
  exp_utility = np.exp(B_time*path_cost) 
  demand_distribution = OD_demand*exp_utility/sum(exp_utility)

  #Compute user costs and distribute flows:
  user_cost = sum(demand_distribution*path_cost)

  return np.sum(demand_distribution[:,None]*link_array[:num_paths,:], axis = 0), user_cost

def headway_setter(transit_routes, link_flows, link_dictionary, load_factor, vehicle_cap, h_min, h_max):
    """
    For each route, the minimum headway is determined so that the load factor is
    satisfied.

    PARAMETERS:
      + transit_routes  --
      + link_flows      --
      + link_dictionary -- 
      + load_factor     --
      + vehicle_cap     --
      + h_min           --
      + h_max           --

    OUTPUT: 
      + all_headways    --
    """
    
    all_headways = np.zeros((len(transit_routes),))
    # segment_flow_dictionary = {}
    
    #For each route:
    for k in range(len(transit_routes)):

        route = transit_routes[k]
        
        #For each segment 
        segment_flows = [] 
        for i in range(0, len(route) - 1):
          link = (k, route[i], route[i+1])
          flow_loc =  link_dictionary[link]
          segment_flows.append(link_flows[flow_loc]) #determine flows
          
        #Compute minimum headway
        if max(segment_flows) > 0:
            all_headways[k] = vehicle_cap*load_factor/max(segment_flows)
        else:
            all_headways[k] = h_max
        
        #Checking headway bounds
        mask = all_headways < h_min
        all_headways[mask] = h_min
        mask = all_headways > h_max
        all_headways[mask] = h_max
    
    return all_headways

def assignment(transit_routes, transit_param):

    """
    [Description]

    PARAMETERS:
      + transit_routes     --
      + transit_param      --

    OUTPUT:
      + unmet_demand       --
      + fixed_transit_cost -- 
      + user_ costs        --
    """
    
    #Get parameters:
    OD_flows       = transit_param["OD_flows"]
    transfer_time = transit_param["transfer_time"]
    transfer_mark = transit_param["transfer_mark"]
    road_edge     = transit_param["road_edge"]
    K_paths       = transit_param["K_paths"]
    B_time        = transit_param["B_time"]
    load_factor   = transit_param["load_factor"]
    vehicle_cap   = transit_param["vehicle_cap"]
    unit_cost_ft  = transit_param["unit_cost_ft"]
    max_node_id   = transit_param["max_node_id"]
    h_min         = transit_param["h_min"]
    h_max         = transit_param["h_max"]
    zone_2_node   = transit_param['zone_2_node']
    
    #Varible to store met demand, unmet demand and user costs:
    # met_demand = 0
    unmet_demand = 0
    user_costs = 0
    
    #Build transit network:
    transit_network, link_map, all_cost, link_dictionary, link_array, link_flows = \
    transit_network_builder(transit_routes, road_edge, transfer_mark, max_node_id, K_paths) 
    
    #Create link dictionary and array, and link_flow arrays
    #link_dictionary, link_array, link_flows = gen_route_maps(transit_routes, K_paths)
    
    #Perform assignment.
    for idx, row in OD_flows.iterrows(): 
        origin_node = zone_2_node[row['origin']]
        destination_node = zone_2_node[row['destination']]
        if ( origin_node in transit_network) and (destination_node in transit_network):
          route_flows, od_user_cost = route_assignment((origin_node, destination_node), transit_network, row['trips'], K_paths, B_time, link_dictionary, link_array.copy(), link_map, transfer_time - transfer_mark, transfer_mark)
          link_flows += route_flows
          user_costs += od_user_cost
          # met_demand += row['trips']
          
        else:
          unmet_demand += row['trips']
    
    #Compute route frequencies and travel disutility
    all_headways = headway_setter(transit_routes, link_flows, link_dictionary, load_factor, vehicle_cap, h_min, h_max)
    
    #Compute fixed transit costs. 
    travel_time = 2*all_cost #multiplication by 2 to capture both travel directions
    fixed_transit_cost = unit_cost_ft*np.sum(travel_time/all_headways)
    
    return unmet_demand, fixed_transit_cost, user_costs#, all_headways, link_flows, unmet_demand



def main():
    #PARAMETERS
    #Network:
    road_edge= {}
    road_edge[(1,3)] = 4
    road_edge[(1,10)] = 2.1
    road_edge[(2,10)] = 2.2
    road_edge[(3,10)] = 2.3
    road_edge[(10,12)] = 20.1
    road_edge[(10,21)] = 20.2
    road_edge[(10,11)] = 20.3
    road_edge[(12,20)] = 40
    road_edge[(20,21)] = 60
    
    transit_param = {}
    transit_param["road_edge"] = road_edge
    transit_param["max_node_id"] = 21
    transit_param["transfer_time"] = 5
    transit_param["transfer_mark"] = 0.001
    
    #Routes
    transit_routes = [[1,10,21], [2,10,11], [3,10,12,20,21], [1,3]]
    
    #Demand:
    import pandas as pd
    
    filename = r'C:\Users\Daniel\Google Drive\Research\Transit_Redesign_NDP\Data\hypothetical_toy_network\toy_network_demand.xlsx'
    transit_param["OD_flows"] = pd.read_excel(filename)
    transit_param["num_ODs"] = len(transit_param["OD_flows"])
    
    #Behavior
    transit_param["K_paths"]       = 2
    transit_param["B_time"]        = -0.1
    
    #Operational paramters:
    transit_param["load_factor"]   = 1
    transit_param["budget"]        = 10
    transit_param["vehicle_cap"]   = 100
    transit_param["bus_speed"]     = 25 
    transit_param["unit_cost_ft"]  = 1
    transit_param["h_min"]  = 0
    transit_param["h_max"]  = 1
    
    #Run model 
    met_demand, fixed_transit_cost, user_costs= \
    assignment(transit_routes, transit_param)
    
    print("Fixed transit cost: ", fixed_transit_cost)
    # print("Headways: ", all_headways)
    # print("Link flows: ", link_flows)
    # print("Unmet demand: ", unmet_demand)
    print("User costs: ", user_costs)

if __name__ == "__main__":
    main()