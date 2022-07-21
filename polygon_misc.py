# This file contains multiple functions involving polygons to complement other operations.
# Currently, it's used for testing/protyping.

import geopandas as gpd
from shapely.geometry import Polygon, mapping, MultiPoint, Point, LineString
from shapely.ops import orient
import matplotlib.pyplot as plt
import time
import numpy as np

def polygon_validity(polygon_to_validate):
    """
    Checking Polygon Validity
    A polygon is valid if it is simple, i.e. has no edge "overlap" or intersection.
    In other words, polygon does not intersect with itself.
    For more details: https://en.wikipedia.org/wiki/Simple_polygon

    Note: This function only checks if polygon is simple and uses a polygon as parameter i.e. not a set of vertices.

    PARAMETERS:
        + polygon_to_validate --

    OUTPUT:
        +
    """

    return polygon_to_validate.is_valid.bool()


def extract_vertices(polygons_to_extract_vertices):
    """
    Extracts vertex information from polygon(s)

        PARAMETERS:
            + polygons_to_extract_vertices --

        OUTPUT: List of points as (x, y) tuples
            + point_list --
    """
    # Using "mapping" function to extract vertex information of geometry
    vertices = mapping(polygons_to_extract_vertices)["coordinates"]

    # List to store Point shapely objects
    point_list = []

    for x in vertices:
        for point in x:
            point_list.append(point)

    # Returning list of vertices
    return point_list


def eliminate_gaps(polygon_to_eliminate_gaps):
    """
    [Description]

    PARAMETERS:
        + polygon_to_eliminate_gaps --
    """
    # Extract vertices of exterior polygon only. We don't want the interior vertices
    # as omitting them would eliminate the "gaps"
    eliminated_gaps = Polygon(
        polygon_to_eliminate_gaps.geometry[0].exterior.coords)
    gpd.GeoSeries(eliminated_gaps).plot()
    #plt.show()


def find_centroid(polygon_to_find_centroid):
    """
    [Description]

    PARAMETERS:
        + polygon_to_find_centroid --

    """
    # Geopandas function that returns a GeoSeries of points for each geometric centroid
    centroids = polygon_to_find_centroid.centroid

    # Plotting polygon with centroid
    base = polygon_to_find_centroid.plot()
    centroids.plot(ax=base, color="black")
    plt.show()


def point_in_polygon(polygon, point):
    """
    [Description]

    PARAMETERS:
         + polygon --
         + point   --

    """
    polygon.contains(point)


def test_point_in_polygon():
    """
    [Description]

    OUTPUT:
        +
    """

    # Time execution testing
    start_time = time.time()

    # Creating simple polygon for testing
    s = gpd.GeoSeries(
        [
            Polygon([(0, 0), (1, 1), (0, 1)])
        ])

    # Running function 100 times
    for _ in range(1000):
        # Point in polygon function
        point_in_polygon(s, Point(0.2, 0.6))

    # Printing total time taken by program to run 100 times
    print("Total time: %s" % (time.time() - start_time))


def polygon_from_vertices(vertices):
    """
    Return GeoSeries with polygon from list of vertices.

    PARAMETERS: 
        + vertices --

    OUTPUT:
        + polygon  --
    """
    polygon = gpd.GeoSeries(Polygon(vertices))

    return polygon

def polygon_from_vertices_crs(vertices, project_crs):
    """
    Return GeoSeries with polygon from list of vertices.

    PARAMETERS:
        + vertices    --
        + project_crs --

    OUTPUT: 
        + polygon     --
    """
    polygon = gpd.GeoSeries(Polygon(vertices), crs= project_crs )
              
    return polygon


def area_from_vertices(vertices):
    """
    Function to get area of polygon from vertex information

    PARAMETERS:
        + vertices -- list of points (vertices)

    OUTPUT:
        + Returns the area of the polygon.
    """
    # Getting the polygon from vertex information
    polygon = polygon_from_vertices(vertices)

    # Getting area and returning it
    return polygon.area[0]

def service_area_geo_validity(polygon_vertices, buffer_polygon, project_crs):
    """
    Given a set of vertices, this function checks that the polygon does not intersect
    with itself. Additionally, for the transit network design problem, it checks that
    it does not violate the 3/4 mile buffer.

    Based on Polygon in Polygon (Binary Predicate) Method previously developed.

    PARAMETERS:
        + polygon_vertices -- list of points (vertices)
        + buffer_polygon -- 3/4 mile buffer polygon
        + project_crs --

    OUTPUTS:
        + Returns False and empty list if SA is not valid,
        otherwise, returns True and SA geoseries 
    """
    # If array has two or less coordinate points, shape is not a polygon. (Line or point)
    if polygon_vertices.shape[0] <= 2:
        #plt.plot(polygon_vertices[:,0], polygon_vertices[:,1], 'r')
        #print(polygon_vertices, polygon_vertices.shape[0])
        return False, []
    
    service_area = polygon_from_vertices_crs(polygon_vertices, project_crs)

    # If is not valid (does not intersect with itself) return False
    if not polygon_validity(service_area):
        return False, []

    # Checks if polygon does not violate the 3/4 mile buffer rule
  
    if service_area.contains(buffer_polygon).bool():
        return True, service_area
    else:
        return False, []
    

def service_area_validity(polygon_vertices, buffer_polygon):
    """
    Given a set of vertices, this function checks that the polygon does not interect
    with itself. Additionally, for the transit network design problem, it checks that
    it does not violate the 3/4 mile buffer.

    Based on Polygon in Polygon (Binary Predicate) Method previously developed.

    PARAMETERS:
        + polygon_vertices  -- list of points (vertices)
        + buffer_polygon    -- 3/4 mile buffer polygon

    OUTPUT:
        + Returns boolean; True if valid, False otherwise
    """
    # If array has two or less coordinate points, shape is not a polygon. (Line or point)
    if polygon_vertices.shape[0] <= 2:
        #plt.plot(polygon_vertices[:,0], polygon_vertices[:,1], 'r')
        #print(polygon_vertices, polygon_vertices.shape[0])
        return False
    
    service_area = polygon_from_vertices(polygon_vertices)

    # If is not valid (does not intersect with itself) return False
    if not polygon_validity(service_area):
        return False

    # Checks if polygon does not violate the 3/4 mile buffer rule
  
    min_requirements = service_area.contains(buffer_polygon).bool()

    return min_requirements

def feasible_space(master_polygon, service_area):
    """
    Function that "clips" the service area polygon where it goes beyond feasible space.
    The master polygon defines the acceptable extension of the service area polygon.
    For example, the master polygon cuts away lakes, seas, etc.

    PARAMETERS:
        + master_polygon       -- shapefile/geodataframe of polygon that defines acceptable extension
        + service_area         -- shapefile/geodataframe of service area polygon

    OUTPUT:
        + clipped_service_area -- new clipped service area polygon as a geodataframe
    """
    
    # First checking if service_area is fully contained in master_polygon.
    # If so, there's no need to clip anything. This could avoid unnecessary
    # calculations and operations.
    fully_contained = service_area.contains(master_polygon).bool()

    if fully_contained:
        return service_area
    
    # To clip the service_area where it goes beyond feasible space, we simply
    # find the intersection of both service_area and master_polygon. This would
    # output the service area polygon where it it's allowed to be extended.
    clipped_service_area = gpd.overlay(service_area, master_polygon, how="intersection")

    return clipped_service_area

def create_min_buffer(routes_path):
    """
    WIP

    PARAMETERS: 
        + routes_path   --

    OUTPUT:
        + routes_buffer -- 

    """
    routes = gpd.read_file(routes_path)
    
    # Storing original crs
    routes_crs = routes.crs
    
    # Reprojecting
    routes = routes.to_crs("EPSG:3991")
    
    # Creating buffer polygon of 3/4 mile for constraint (3/4 mile = 3960 ft)
    routes_buffer = routes.buffer(3960)
    
    # Reprojecting buffer back to original CRS
    routes_buffer = routes_buffer.to_crs(routes_crs)
    
    # Combining all buffers to form one single buffer polygon around routes
    routes_buffer = routes_buffer.unary_union
    
    return routes_buffer

def buffer_from_routes(network_routes, transit_param):
    """
    Function generates a buffer from the routes in the transit network.
    
    PARAMETERS:
        + network_routes -- list with the routes in the network.
        + transit_param  -- dictionary with the key parameters
    
    OUTPUT:
        + routes_buffer  -- geopandas object with the route buffer
        + min_sa_array   --
        + min_sa_geo     --
        
    """
    
    #Get parameters
    min_buffer_size = transit_param['min_buffer_size']
    project_crs = transit_param['project_crs']
    vertex_coordinates = transit_param['vertex_coordinates']
    vertex_orientation = transit_param['vertex_orientation']
    
    all_lines = {'geometry': []}
    
    #For each route
    for route in network_routes:
        
        #Collect coordinates
        coordinates = []
        
        for vertex in route:
            coordinates.append(vertex_coordinates[vertex])
        
        #Collect in dictionary as a LineString object.
        all_lines['geometry'].append(LineString(coordinates))
    
    #Create geopandas object with the lines
    gdf = gpd.GeoDataFrame(all_lines, crs= project_crs)
        
    #Generate buffer
    routes_buffer = gdf.buffer(min_buffer_size)

    # Combining all buffers to form one single buffer polygon around routes
    # and orienting vertices to counter_clockwise.
    routes_buffer = orient(routes_buffer.unary_union, sign= vertex_orientation)

    # Extracting only exterior polygon (eliminating inner polygons/holes)
    buffer_vertices = routes_buffer.exterior.coords
    
    # Buffer vertices as a 2D numpy array
    min_sa_array = np.array(buffer_vertices)
    min_sa_geo =  polygon_from_vertices_crs(min_sa_array, project_crs)
    return gpd.GeoSeries(routes_buffer, crs= project_crs), min_sa_array, min_sa_geo    
    
   

# def main():
#     test_polygons = gpd.read_file("../Shapefiles/polygons_valid_invalid.shp")

#     simple_polygon = gpd.read_file("../Shapefiles/PSA_test.shp")

#     bad_polygon = gpd.read_file("../Shapefiles/PSA_bad.shp")

#     routes = gpd.read_file(
#         "../Shapefiles/geo_export_7f388299-10d0-4741-b68e-af3a3f608775.shp")
#     routes = routes.to_crs("EPSG:3991")
#     # Creating buffer polygon of 3/4 mile for constraint (3/4 mile = 3960 ft)
#     routes_buffer = routes.buffer(3960)
#     routes_buffer = gpd.GeoDataFrame(geometry=routes_buffer)
#     routes_buffer = routes_buffer.dissolve()
    
    # Function to check the validity of polygons
    # polygon_validity(test_polygons)

    # Function to extract vertex information
    # extract_vertices(routes_buffer)

    # Function to eliminate gaps (interior polygons) from polygon
    # eliminate_gaps(routes_buffer)

    # Function to find centroid of polygon(s)
    # find_centroid(simple_polygon)

    # Testing execution time of point in polygon function
    # test_point_in_polygon()

    #vertices_test = extract_vertices(routes_buffer)
    #print(type(polygon_from_vertices(vertices_test)))

    # Function to get area from polygon given its vertex information
    # simple_vertices = extract_vertices(simple_polygon)
    # simple_vertices = np.array(simple_vertices)
    # plt.scatter(simple_vertices[:, 0], simple_vertices[:, 1])
    # plt.show()

    # Function to check validity of polygon. (Does not intersect with itself and does not violaate 3/4 mile buffer)
    # invalid_polygon = gpd.read_file("../Shapefiles/Invalid_PSA.shp")
    # simple_vertices = extract_vertices(invalid_polygon)
    # print(service_area_validity(simple_vertices, routes_buffer))

    # Function to clip service area polygon beyond feasible space using master polygon
    # master_polygon = gpd.read_file("../Shapefiles/Master_Polygon.shp")
    # clipped_service_area = feasible_space(master_polygon, simple_polygon)
    # clipped_service_area.plot()
    # plt.show()

# if __name__ == "__main__":
#     main()
