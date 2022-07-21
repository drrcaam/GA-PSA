def is_clockwise(polygon):
  """ 
  Tells if polygon vertices are in clockwise order or not

  Parameters: 
      + polygon -- Shapely Polygon
      
  OUTPUT:
      Returns True if clockwise, False if counter-clockwise
  """
  # Getting vertices
  coords = list(zip(*polygon.exterior.xy))

  # Implementation of top solution from https://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
  sum = 0
  for i in range(len(coords)):
    if i < len(coords)-1:
      sum = sum + (coords[i+1][0] - coords[i][0])*(coords[i+1][1] + coords[i][1])
  
  return True if sum >= 0 else False