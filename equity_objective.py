# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 17:20:25 2022

@author: Owner
"""
import polygon_misc
import geopandas as gpd
import pandas as pd
import numpy as np

def compute_pop_group(polygon, transit_param):
    """
    Function computes the number of people in each population group that is 
    within the service area,
    
    PARAMETERS:
        + polygon       -- 2D array with the coordinates of the service area
        + transit_param -- dictionary with parameters
        
    OUTPUT:
        + pop_container -- 1D array with the populations groups inside service area
    
    """
    
    #Get parameters 
    zonal_population = transit_param["pzones_data"]
    zonal_shp = transit_param["zonal_data_shp"]
    project_crs = transit_param["project_crs"]
    spatial_group_column = transit_param["spatial_group_column"]
    county_list = transit_param["county_list"]
    pop_container = transit_param["pop_container"].copy()
    pop_group = transit_param["pop_group"]
    
    # Generating the polygon GeoSeries from vertex information [Second most time consuming step]
    polygon_geo = polygon_misc.polygon_from_vertices(polygon)
    polygon_geo.crs = project_crs

    # Converting polygon to a GeoDataFrame to perform intersection
    polygon_gdf = gpd.GeoDataFrame(geometry=polygon_geo,crs=project_crs)
    
    # Intersecting polygon(service area) and zones [THE MOST TIME CONSUMING COMPONENT]
    resulting_intersect = gpd.overlay(zonal_shp, polygon_gdf, how = "intersection",keep_geom_type=False)

    # Caluclating area of each zone of the interesecting polygon
    resulting_intersect["intersect_area"] = resulting_intersect.area

    # For each zone k, calculating the proportion of k's area that fall within SA and adding it as a field/column
    resulting_intersect["Proportion"] = resulting_intersect["intersect_area"] / \
        resulting_intersect["area"]
    
    # Merging zonal population and zonal proportion data
    zonal_pop_prop = pd.merge(zonal_population, resulting_intersect, how='inner', on=['GEOID'])
    
    #Computing size of population of interest by county 
    zonal_pop_prop["pop_size"] = zonal_pop_prop[pop_group] * zonal_pop_prop["Proportion"]
    agg_zonal_pop = zonal_pop_prop.groupby(spatial_group_column)["pop_size"].sum()

    for k in range(0,len(county_list)):
        if county_list[k] in agg_zonal_pop.index:
            pop_container[k] = agg_zonal_pop.loc[county_list[k]]
        else: 
            pop_container[k] = 0
  
    return pop_container

    
def spatial_proximity(sa_array, bus_pop_group, transit_param):
    """
    Compute the Atkinson inequality index based on the change in the sizes of the 
    populations that have spatial access to the services. 
    
    PARAMETERS:
        + sa_array      -- 2D array that contains the current size of the service area.
                           If empty, it means that the analysis is only based on the min_array.
                           If not empty, it is used to compute the size of the populations
                           with access to the paratransit service.
        + bus_pop_group --
        + transit_param -- dictionary with relevant parameters

    OUTPUT: 
        + atkinson      -- number, the atkinson index
    """
    
    #Parameters
    init_pop_group = transit_param['init_pop_group'] #Initial number of users proximal to services
    power = 1 - transit_param['atkinson_epsilon']    
    
    #Compute coverage of new designs
    paratransit_pop_group = compute_pop_group(sa_array, transit_param)
        
    new_pop_group = np.append(paratransit_pop_group,bus_pop_group)
    
    #Compute frac difference:"
    diff = (new_pop_group - init_pop_group)/init_pop_group
    
    #Compute Atkinson index
    atkinson = 1 - (1/np.mean(diff)) * (np.sum(diff**power)/len(diff))**(1/power)
    
    return atkinson


def compute_pop_group_geo(polygon, transit_param):
    """
    Function computes the number of people in each population group that is 
    within the service area,
    
    PARAMETERS:
        + polygon       -- geodataframe object representing service area
        + transit_param -- dictionary with parameters
        
    OUTPUT:
        + pop_container -- 1D array with the populations groups inside service area
    
    """
    
    #Get parameters 
    zonal_population = transit_param["pzones_data"]
    zonal_shp = transit_param["zonal_shp"]
    spatial_group_column = transit_param["spatial_group_column"]
    county_list = transit_param["county_list"]
    pop_container = transit_param["pop_container"].copy()
    pop_group = transit_param["pop_group"]

    # Intersecting polygon(service area) and zones [THE MOST TIME CONSUMING COMPONENT]
    resulting_intersect = gpd.overlay(zonal_shp, polygon, how = "intersection",keep_geom_type=False)

    # Caluclating area of each zone of the interesecting polygon
    resulting_intersect["intersect_area"] = resulting_intersect.area

    # For each zone k, calculating the proportion of k's area that fall within SA and adding it as a field/column
    resulting_intersect["Proportion"] = resulting_intersect["intersect_area"] / \
        resulting_intersect["area"]
    
    # Merging zonal population and zonal proportion data
    zonal_pop_prop = pd.merge(zonal_population, resulting_intersect, how='inner', on=['GEOID'])
    
    #Computing size of population of interest by county 
    zonal_pop_prop["pop_size"] = zonal_pop_prop[pop_group] * zonal_pop_prop["Proportion"]

    agg_zonal_pop = zonal_pop_prop.groupby(spatial_group_column)["pop_size"].sum()

    for k in range(0,len(county_list)):
        if county_list[k] in agg_zonal_pop.index:
            pop_container[k] = agg_zonal_pop.loc[county_list[k]]
        else: 
            pop_container[k] = 0
  
    return pop_container

    
def spatial_proximity_geo(bus_polygon, paratransit_polygon, transit_param):
    """
    Compute the Atkinson inequality index based on the change in the sizes of the 
    populations that have spatial access to the services. 
    
    PARAMETERS:
        + bus_polygon         -- object representing the minimum service area for bus service
        + paratransit_polygon -- object represeting the paratransit service area
        + transit_param       -- dictionary with relevant parameters
        

    OUTPUT:
        + atkinson            -- number; atkinson index
        + bus_pop_group       -- 1D array; population group quantities with bus coverage
    """
    
    #Parameters
    total_pop_group = transit_param['total_pop_group'] 
    power = transit_param['atkinson_epsilon']    
    
    #Compute coverage of new designs
    paratransit_pop_group = compute_pop_group_geo(paratransit_polygon, transit_param)
    bus_pop_group = compute_pop_group_geo(bus_polygon, transit_param)

    joined_pop_group = np.append(paratransit_pop_group,bus_pop_group)
    
    #Compute frac difference:
    coverage_proportion = joined_pop_group/total_pop_group
    
    #Compute Atkinson index
    atkinson = 1 - (1/np.mean(coverage_proportion)) * (np.sum(coverage_proportion**power)/len(coverage_proportion))**(1/power)
    
    return atkinson, bus_pop_group 

def spatial_proximity_geo_update(polygon, bus_pop_group, transit_param):
    """
    Compute the Atkinson inequality index based on the change in the sizes of the 
    populations that have spatial access to the services. 
    
    PARAMETERS:
        + polygon               -- geo_pandas object representing service area
        + bus_pop_group         -- 1D array; population group quantities with bus coverag
        + transit_param         -- dictionary with relevant parameters

    OUTPUT: 
        + atkinson              -- number; atkinson index
        + paratransit_pop_group -- 1D array; population group quantities with paratransit coverag
    """
    
    #Parameters
    total_pop_group = transit_param['total_pop_group'] 
    power = 1 - transit_param['atkinson_epsilon']    
    
    #Compute coverage of new designs
    paratransit_pop_group = compute_pop_group_geo(polygon, transit_param)

    joined_pop_group = np.append(paratransit_pop_group,bus_pop_group)
    
    #Compute frac difference:
    coverage_proportion = joined_pop_group/total_pop_group
    
    #Compute Atkinson index
    atkinson = 1 - (1/np.mean(coverage_proportion)) * (np.sum(coverage_proportion**power)/len(coverage_proportion))**(1/power)
    
    return atkinson, paratransit_pop_group 


def atkison_updater(bus_pop_group, paratransit_pop_group, transit_param):
    """
    Compute the Atkinson index based on previous group information:
        
    PARAMETERS:
        + bus_pop_group         -- 1D array with info on demographic groups covered by ft design
        + paratransit_pop_group -- 1D array with info on demographic groups covered by paratransit
        + transit_param         -- dictionary with parameters:
            
    OUTPUT:
        + atkinson               -- number; atkinson index
    """
    init_pop_group = transit_param['init_pop_group'] #Initial number of users proximal to services
    power = 1 - transit_param['atkinson_epsilon']    
    
    #Merge
    new_pop_group = np.append(paratransit_pop_group,bus_pop_group)
    
    #Compute frac difference:"
    diff = (new_pop_group - init_pop_group)/init_pop_group
    
    #Compute Atkinson index
    atkinson = 1 - (1/np.mean(diff)) * (np.sum(diff**power)/len(diff))**(1/power)
    
    return atkinson