import numpy as np
import matplotlib.pyplot as plt
import polygon_misc
from shapely.geometry import Polygon

def sa_crossover(parent_1, parent_2, ga_param):
    """Crossover Function
    The next functions implement a crossover operator that swaps segments of service area borders between different service area designs. The basic idea is that a "slices" of vertices can be swapped.
    
    Steps:
    
    Given are two sets of vertices for two service area designs (the parents) and a universal common central location (e.g., the location of a major hospital that must be in any final service area design).
    Compute the angle of all service area vertices for both desings.
    Select a starting angle and an ending angle for the slice.
    For each set of vertices, select the vertex  vs  with angle closest to the starting angle and select the vertex  ve  with angle closest to the ending angle.
    Select and swap the set of vertices bound by the  (vs,ve)  vertices.

    PARAMETERS:
      + parent_1 --
      + parent_2 --
      + ga_param --

    OUTPUT: 
      + child_1 --
      + child_2 --
    """
    #Get parameters
    min_angle_slice = ga_param["min_angle_slice"] 
    angle_slice_span = ga_param["angle_slice_span"] 
    
    # Getting universal common central point
    # Converting polygon vertices to Shapely Polygon to perfrom spatial analysis
    parent1_polygon = Polygon(parent_1)
    parent2_polygon = Polygon(parent_2)
    # Checking if parents intersect. If False, there is no universal common central location and parents are returned
    if parent1_polygon.intersects(parent2_polygon):
        resulting_intersect = parent1_polygon.intersection(parent2_polygon)
        centroid = resulting_intersect.centroid
        central_point = np.array([centroid.x, centroid.y])
    else:
        # Parents dont intersect meaning there is no universal common central location.
        centroid = parent1_polygon.centroid
        central_point = np.array([centroid.x, centroid.y])
    
    #Compute the angle of all service area vertices for both desings.
    coor_difference = parent_1 - central_point
    angle_parent_1 = np.arctan2(coor_difference[:,1], coor_difference[:,0])
    coor_difference = parent_2 - central_point
    angle_parent_2 = np.arctan2(coor_difference[:,1], coor_difference[:,0])
    #print(angle_parent_2)
    
    #Select a starting angle for the slice.
    start_angle = np.random.random()*2*np.pi
    #print(start_angle)
    #Select a ending angle for the slice.
    slice_size = min_angle_slice + angle_slice_span*np.random.random()
    end_angle = start_angle + slice_size
    if end_angle > np.pi:
      end_angle = end_angle - 2*np.pi 
    
    #For each set of vertices, select the vertex  v_s  with angle closest 
    #to the starting angle 
    #and select the vertex  v_e  with angle closest to the ending angle.
    start_vertex_parent_1 = np.argmin(np.abs(start_angle - angle_parent_1)) #v_s
    end_vertex_parent_1   = np.argmin(np.abs(end_angle - angle_parent_1))#v_e
    start_vertex_parent_2 = np.argmin(np.abs(start_angle - angle_parent_2))#v_s
    end_vertex_parent_2   = np.argmin(np.abs(end_angle - angle_parent_2))#v_e
    
    #The vertex sets are ordered in counter_clockwise orientation.
    #If the end vertex has a lower value than the start vertex, then the origin 
    #vertex is in between both vertices, and then the order of the swap must be changed
    if start_vertex_parent_1 > end_vertex_parent_1: #origin in between
      vertex_swap_1 = np.vstack((parent_1[start_vertex_parent_1:,:], parent_1[:end_vertex_parent_1+1,:])) #from start vertex to end, and from the origin to the end vertex
      flag_swap_1 = True
      #print("complex 1")
    else:
      vertex_swap_1 = parent_1[start_vertex_parent_1:end_vertex_parent_1+1,:] 
      flag_swap_1 = False
      #print("normal 1")
    
    if start_vertex_parent_2 > end_vertex_parent_2:
      vertex_swap_2 = np.vstack((parent_2[start_vertex_parent_2:,:], parent_2[:end_vertex_parent_2+1,:])) 
      flag_swap_2 = True
      #print("complex 2")
    else:
      vertex_swap_2 = parent_2[start_vertex_parent_2:end_vertex_parent_2+1,:] 
      flag_swap_2 = False
      #print("normal 2")
    
    #Swap information
    if flag_swap_1 == True:
      child_1 = np.vstack( (vertex_swap_2, parent_1[end_vertex_parent_1:start_vertex_parent_1,:])  )
    else: 
      child_1 = np.vstack( (parent_1[:start_vertex_parent_1,:], vertex_swap_2, parent_1[end_vertex_parent_1+1:,:])  )
    
    if flag_swap_2 == True:
      child_2 = np.vstack( (vertex_swap_1, parent_2[end_vertex_parent_2:start_vertex_parent_2,:])  )
    else: 
      child_2 = np.vstack( (parent_2[:start_vertex_parent_2,:], vertex_swap_1, parent_2[end_vertex_parent_2+1:,:])  )
    
    return child_1, child_2