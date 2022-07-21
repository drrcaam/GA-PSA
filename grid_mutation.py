import numpy as np
import matplotlib.path as mpltPath
import pandas as pd
from operator import add, sub, mul

def rotator(x, y, angle):
    """
    Adjust coordinates given rotation in reference system

    PARAMETERS:
      + x         -- 2D array with x coordinates
      + y         -- 2D array with y coordinates
      + angle     -- angle of rotation

    OUTPUT:
      + rotated_x -- 2D array with x coordinates
      + rotated_y -- 2D array with y coordinates
    """
    rotated_x = x*np.cos(angle) - y*np.sin(angle)
    rotated_y = x*np.sin(angle) + y*np.cos(angle)
    return rotated_x, rotated_y

def pick_starting_vertex(sa_array, reference_point, angle_scope):
  """
  Relative to a reference point, the angle of each polygon vertex is computes.
  A angle range is defined. 
 

  PARAMETERS:
    + sa_array        -- 2D array with the polygon vertices
    + reference_point -- point from where angles are measured.
    + angle_scope     -- number that defines the "cone" of vertex selection
  
  OUTPUT:
    randomly selected vertex
  """
  #Generate indices
  idx = np.arange(sa_array.shape[0])

  #Generate angles
  coor_difference = sa_array - reference_point
  angle = np.arctan2(coor_difference[:,1], coor_difference[:,0])


  while True: 
    #Find angle range from where vertex will be computed. 
    starting_angle = -np.pi + 2*np.pi*np.random.random()
    ending_angle = starting_angle + angle_scope
    
    if ending_angle <= np.pi:
      #make selection of angles.
      mask =(angle >= starting_angle) & (angle <= ending_angle)
      #plt.scatter(sa_array[mask,0], sa_array[mask,1])
      return np.random.choice(idx[mask])

    else:
      ending_angle = ending_angle - 2*np.pi #this angle will be negative
      mask = ((angle <= ending_angle) & (angle >= -np.pi) ) | (angle >= starting_angle)
      #plt.scatter(sa_array[mask,0], sa_array[mask,1])
      return np.random.choice(idx[mask])

def distance_first_vertex(sa_array):
  """
  Function selects the first vertex in the segment to mutate by randomly generating
  a distance (anchor number) within the length of the polygon boundary and selecting the
  vertex closest to that distance.  
  
  PARAMETERS:
    + sa_array        -- 2D array with the polygon vertices

  
  OUTPUT:
    index of vertex
  """
  #Generate indices
  num_vertex = sa_array.shape[0]

  #boundary_distance = np.cumsum((np.sum((sa_array[:num_vertex-1,:] - sa_array[1:,:])**2,axis = 1)))
  boundary_distance = np.cumsum((np.sum(np.absolute((sa_array[:num_vertex-1,:] - sa_array[1:,:])),axis = 1)))#manhattan distance
  spot_in_boundary = np.random.random()*boundary_distance[-1]
  diff = np.absolute(boundary_distance - spot_in_boundary)

  return np.argmin(diff)

def grid_random_walk(mesh_X, mesh_Y,  grid_num_pts, expansion_mode, map_boundary):
  """
  Function that performs the random walk across grid points.

  PARAMETERS:
    + mesh_X         -- x coordinates of grid
    + mesh_Y         -- y coordinates of grid
    + grid_num_pts   -- number of points in marginals of grids
    + expansion_mode -- True if in expansion, False if contraction
    + map_boundary   -- indicates which grid points are inside and outside polygon

  OUTPUT: 
    + new_boundary   -- list with new set of boundary points

  """
  #Initializing new boundary
  new_boundary = []

  #Switching map if contraction.
  if expansion_mode == True:
    map_boundary = map_boundary == False

  for j in range(grid_num_pts):
    #Select points in polygon side of interest.

    mask = map_boundary[:,j] == True
    selected_Y = mesh_Y[mask,j]
    selected_X = mesh_X[mask,j]

    if len(selected_Y) > 0: 
      x_1 = selected_X[0] 
      x_2 = selected_X[-1]
      y_1 = selected_Y[0]
      y_2 = selected_Y[-1]

      #Random point in column, which is a new vertex in the boundary 
      w = np.random.random()
      new_x = x_1 + w*(x_2 - x_1) 
      new_y = y_1 + w*(y_2 - y_1) 

      new_boundary.append([selected_X[0], new_y])

  return np.array(new_boundary)

def grid_random_walk_coordinator(points, sa_array, expansion_mode, grid_parameters):

  """
  [Description]

  PARAMETERS:
    + points          --
    + sa_array        --
    + expansion_mode  --
    + grid_parameters --

  OUTPUT:
    + new_boundary    --
  """

  #Get parameters: 
  box_size_buffer = grid_parameters["box_size_buffer"]
  grid_num_pts = grid_parameters["grid_num_pts"]
  grid_point_2D_space = grid_parameters["grid_point_2D_space"].copy()

  #Generate grid:
  mesh_X, mesh_Y = grid_generator(points, box_size_buffer, grid_num_pts)

  #Determine with if grid points are inside polygon. 
  grid_point_2D_space[:, 0] = mesh_X.reshape((grid_point_2D_space.shape[0],))
  grid_point_2D_space[:, 1] = mesh_Y.reshape((grid_point_2D_space.shape[0],))

  #plt.scatter(grid_point_2D_space[:, 0], grid_point_2D_space[:, 1])

  path = mpltPath.Path(sa_array)
  inside_polygon = path.contains_points(grid_point_2D_space)
  map_boundary = inside_polygon.reshape((grid_num_pts, grid_num_pts)) 

  #Perform random walk.
  new_boundary = grid_random_walk(mesh_X, mesh_Y,  grid_num_pts, expansion_mode, map_boundary)

  return np.array(new_boundary)

def grid_generator(points, box_size_buffer, grid_num_pts):

  """
  Function that generates the grid points

  PARAMETERS:
    + points          -- array with list of pints
    + box_size_buffer -- sets size of grid points box
    + grid_num_pts    -- number of grid points

  OUTPUT:
    + mesh_x          -- 2D array with x coordinates
    + mesh_y          -- 2D array with y coordinates
  """

  diff =  points[-1,:] - points[0,:] 
  angle = np.arctan2( diff[1], diff[0])
  
  if angle >= 0:
    angle = np.pi - angle
  elif angle < 0:
    angle = -1*angle

  #Rotate points  
  adj_points = points - points[0,:]
  rotated_x, rotated_y = rotator(adj_points[:,0], adj_points[:,1], angle)

  #Find minimum and maximums 
  y_min = np.min(rotated_y, axis = 0) - box_size_buffer 
  y_max = np.max(rotated_y, axis = 0) + box_size_buffer

  #Create mesh
  grid_x = np.linspace(rotated_x[0], rotated_x[-1], grid_num_pts)
  grid_y = np.linspace(y_min, y_max, grid_num_pts)
  X, Y = np.meshgrid(grid_x, grid_y)

  #Rotate mesh and re-adjust the x's and the y's
  mesh_x, mesh_y = rotator(X, Y , -1*angle)
  mesh_x = mesh_x + points[0,0]
  mesh_y = mesh_y + points[0,1]

  return mesh_x, mesh_y
  #plt.plot(points[:,0], points[:,1], 'r')
  #plt.scatter(mesh_x, mesh_y)
  #plt.scatter(rotated_x, rotated_y)

def column_inside_polygon(column_X, column_Y, path, expansion_mode, grid_point_2D_space):
  """
  Function determines if points in a column are inside or outside the polygon.
  It returns a boolean array that indicates which points can be visited (True) in the 
  mutation depending on whether we are expanding or contracting the polygon. 

  PARAMETERS:
    + column_X            -- 1D array with x coordinates
    + column_Y            -- 1D array with y coordinates
    + path                -- object that contains the polygon
    + expansion_mode      -- Boolean; determines if operations is expansion or contraction 
    + grid_point_2D_space -- 2D array that contains coordinates

  OUTPUT:
    + map_boundary        --

  """

  grid_point_2D_space[:,0] = column_X
  grid_point_2D_space[:,1] = column_Y
  map_boundary = path.contains_points(grid_point_2D_space)
  
  #Switching map if contraction.
  if expansion_mode == True:
    map_boundary = map_boundary == False
  return map_boundary

def closest_coordinate(all_x, x_i):
  """
  Determines the closest coordinate to a given x

  PARAMETERS:
    + all_x -- 1D array with all the x coordinates
    + x_i   -- coordinate of interest

  OUTPUT:
    position of the coordinate closest to x_i
  """

  return np.argmin(np.absolute(all_x - x_i))

def random_point_selection(column_X, column_Y, map_boundary):
  """
  Function selects the new vertex location

  PARAMETERS:
    + column_X     -- 1D array with x coordinate
    + column_Y     -- 1D array with y coordinate
    + map_boundary -- 2D array with boolean values that say if point is in or out of polygon

  OUTPUT:
    + new_x        -- x coordinate for new vertex
    + new_y        -- y coordinate for new vertex
    + anchor_idx   -- index closest to the selected y coordinate
  """
  
  selected_X = column_X[map_boundary]
  selected_Y = column_Y[map_boundary]

  if len(selected_Y) > 0: 
    x_1 = selected_X[0]
    x_2 = selected_X[-1]
    y_1 = selected_Y[0]
    y_2 = selected_Y[-1]

    #Random point in column, which is a new vertex in the boundary 
    w = np.random.random()
    new_x = x_1 + w*(x_2 - x_1) 
    new_y = y_1 + w*(y_2 - y_1) 

    #Identify the index closest to the selected y coordinate
    anchor_idx = closest_coordinate(selected_Y, new_y)

    return new_x, new_y, anchor_idx

  #This function can return None. Is this right?

def bounded_point_selection(column_X, column_Y, map_boundary, grid_num_pts, anchor_idx, window_length):
  
  """
  [Description]

  PARAMETERS:
    + column_X     -- 1D array with x coordinate
    + column_Y     -- 1D array with y coordinate
    + map_boundary -- 2D array with boolean values that say if point is in or out of polygon
    + grid_num_pts  -- number of points in grid
    + anchor_idx    -- index closest to the selected y coordinate
    + window_length -- number

  OUTPUT:
    + new_x         -- x coordinate for new vertex
    + new_y         -- y coordinate for new vertex
    + anchor_idx    -- index closest to the selected y coordinate
  """

  #Select grid points in space of interest
  selected_X = column_X[map_boundary]
  selected_Y = column_Y[map_boundary]

  if len(selected_Y) > 0: 
    #Determine exploration window. 
    min_idx = max(anchor_idx - window_length, 0)
    max_idx = min(anchor_idx + window_length,  grid_num_pts - 1)

    #if the space defined by the exploration window around the anchor is feasible.
    if map_boundary[min_idx] and map_boundary[max_idx]:
      #Select coordinates
      x_1 = column_X[min_idx]
      x_2 = column_X[max_idx]
      y_1 = column_Y[min_idx]
      y_2 = column_Y[max_idx]
      
    #else, if the space defined by window created by the bottom and ancho idexes is feasible:
    elif map_boundary[min_idx] and map_boundary[anchor_idx]:
      x_1 = column_X[min_idx]
      x_2 = column_X[anchor_idx]
      y_1 = column_Y[min_idx]
      y_2 = column_Y[anchor_idx]

    #else, if the space defined by window created by the bottom and ancho idexes is feasible:
    elif map_boundary[anchor_idx] and map_boundary[max_idx]:
      x_1 = column_X[anchor_idx]
      x_2 = column_X[max_idx]
      y_1 = column_Y[anchor_idx]
      y_2 = column_Y[max_idx]
      
    #Using space that is feasible 
    else: 
      x_1 = selected_X[0]
      x_2 = selected_X[-1]
      y_1 = selected_Y[0]
      y_2 = selected_Y[-1]

    #Random point in column, which is a new vertex in the boundary 
    w = np.random.random()
    new_x = x_1 + w*(x_2 - x_1) 
    new_y = y_1 + w*(y_2 - y_1) 

    #Identify the index closest to the selected y coordinate
    anchor_idx = closest_coordinate(selected_Y, new_y)
    
    return new_x, new_y, anchor_idx

  #This function can return None. Is this right?
  
def grid_c_walk(points, sa_array, expansion_mode, grid_parameters):
  """
  Function that performs the controlled random walk across grid points.

  PARAMETERS:
    + points          -- 2D array with selected polygon vertices to mutate
    + sa_array        -- 2D array with polygon vertices
    + expansion_mode  -- True if in expansion, False if contraction
    + grid_parameters -- dictionary with parameters used in mutation process

  OUTPUT: 
    + new_boundary    -- 2D array with new set of boundary points

  """
  #Get parameters: 
  box_size_buffer = grid_parameters["box_size_buffer"]
  grid_num_pts = grid_parameters["grid_num_pts"]
  grid_point_2D_space = grid_parameters["grid_point_2D_space"].copy()
  p_col_in = grid_parameters["p_col_in"]
  window_length = grid_parameters["window_length"]

  #Generate grid:
  mesh_X, mesh_Y = grid_generator(points, box_size_buffer, grid_num_pts)
  #plt.scatter(mesh_X.reshape((mesh_X.size,)), mesh_Y.reshape((mesh_Y.size,)))
  
  #Create polygon path object
  path = mpltPath.Path(sa_array) 
  
  #PERFORM RANDOM WALK
  
  #Initializing new boundary
  new_boundary = []

  #Find first column
  found = False
  for j in range(grid_num_pts):

    map_boundary = column_inside_polygon(mesh_X[:,j], mesh_Y[:,j], path, expansion_mode, grid_point_2D_space)
    new_point = random_point_selection(mesh_X[:,j], mesh_Y[:,j], map_boundary)
    if new_point != None:
      #plt.scatter(mesh_X[:,j], mesh_Y[:,j], alpha = 0.2)
      new_boundary.append([new_point[0], new_point[1]])
      anchor_idx = new_point[2]
      found = True
      break

  #If feasible space found
  if found == True:
    
    #Find points in the other columns
    for k in range(j+1, grid_num_pts - 1): # whether a column is considered or not is probabilistic, but the last column is alwaya considered.
      if np.random.random() < p_col_in:
        
        map_boundary = column_inside_polygon(mesh_X[:,k], mesh_Y[:,k], path, expansion_mode, grid_point_2D_space)
        new_point = bounded_point_selection(mesh_X[:,k], mesh_Y[:,k], map_boundary, grid_num_pts, anchor_idx, window_length)
        if new_point != None:
          #plt.scatter(mesh_X[:,k], mesh_Y[:,k], alpha = 0.2)
          new_boundary.append([new_point[0], new_point[1]])
          anchor_idx = new_point[2]

    #Find point in last columns:
    
    map_boundary = column_inside_polygon(mesh_X[:,-1], mesh_Y[:,-1], path, expansion_mode, grid_point_2D_space)
    new_point = bounded_point_selection(mesh_X[:,-1], mesh_Y[:,-1], map_boundary, grid_num_pts, anchor_idx, window_length)
    if new_point != None:
      #plt.scatter(mesh_X[:,-1], mesh_Y[:,-1], alpha = 0.2)
      new_boundary.append( [new_point[0], new_point[1]])

  return np.array(new_boundary)

def grid_mutation_coordinator_v2(sa_array, mut_parameters, grid_parameters):
  """
  The function that coordinates the local mutation process. Its main tasks are to:
    a) determine how many mutations are performed.
    b) select which nodes are selected for mutation.

  PARAMETERS:
    + sa_array        -- two dimensional array representing a service area polygon
    + mut_parameters  -- dictionary with mutation parameters
    + grid_parameters -- parameters used to coordinate grid mutations

  OUTPUT:
    + sa_polygon  -- mutated sa_polygon
  """

  #Get parameters
  min_num_mut = mut_parameters["min_num_mut"]
  max_num_mut = mut_parameters["max_num_mut"]
  num_vertices = sa_array.shape[0]
  min_selection_size = int(num_vertices*mut_parameters["min_selection_size"]) #COULD CHANGE!!!!
  max_selection_size = int(num_vertices*mut_parameters["max_selection_size"]) #COULD CHANGE!!!!
  p_expansion = mut_parameters["p_expansion"]
  
  # Probability for expansion to occur
  if np.random.random() < p_expansion:
      expansion_mode = True
  else:
      expansion_mode = False

  #Other parameters
  num_points = sa_array.shape[0]
  x_min, y_min = np.min(sa_array, axis = 0)
  x_max, y_max = np.max(sa_array, axis = 0)
  
  #Determine number of mutations
  num_mut = np.random.randint(min_num_mut, max_num_mut+1)
 
  #Perform mutations
  for k in range(0,num_mut):

    # If array has two or less coordinate points, shape is not a polygon. (Line or point)
    # This is done to avoid program from erroring out, however, this does not solve the
    # problem of having too few vertices as solutions for Genetic Algorithm.
    if sa_array.shape[0] <= 2:
        #print("Grid Mutation: Eliminated too many vertices while mutating.")
        break;
        
    #Select region of polygon that will be mutated
    start_point = distance_first_vertex(sa_array)
    num_points2include = np.random.randint(min_selection_size, max_selection_size) 
    
    #If the point selection is in the middle of the polygon coordinates
    if start_point + num_points2include <  num_points:
        end_point = start_point + num_points2include

        #Select the polygon segment to be mutated:
        points = sa_array[start_point:end_point+1,:].copy()

        #Perform local mutation 
        new_boundary = grid_c_walk(points, sa_array, expansion_mode, grid_parameters)
 
        #Incorporate mutation to polygon. 
        if new_boundary.size != 0:
          sa_array = np.vstack(( sa_array[:start_point,:], new_boundary, sa_array[min(num_points,end_point+1):,:]))
        else:
          sa_array = np.vstack(( sa_array[:start_point,:], sa_array[min(num_points,end_point+1):,:]))
    
    else:
        #Find the end point
        end_point = num_points2include - (num_points - start_point) 

        points = np.vstack((sa_array[start_point:num_points,:], sa_array[:end_point,:]) )  

        #Perform local mutation 
        new_boundary = grid_c_walk(points, sa_array, expansion_mode, grid_parameters)

        #Incorporate mutation to polygon. 
        if new_boundary.size != 0:
          sa_array = np.vstack(( new_boundary, sa_array[end_point:start_point,:]))
        else:
          sa_array = sa_array[end_point:start_point,:]
          
  return sa_array