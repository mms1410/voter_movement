#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 16:37:40 2022

@author: svenmaurice
"""
"""
ToDo:
    - discretize polygon area into points
    - population in each constituency -> sample points
    - radius point whole polygone
"""

from os.path import dirname  as dirname
from os import getcwd
from os.path import join as join
#from sys import path as spath
##
import pandas as pd
import numpy as np
import datetime
from dateutil.relativedelta import relativedelta
#
import geopandas as gpd
import shapely
import geopy
##  
import twint

#spath.append(getcwd()) ## in spyder make sure to be really in cwd
#import helper_functions
###############################################################################
PATH_SCRIPT = getcwd()
PATH_BASE = dirname(dirname(PATH_SCRIPT))   
PATH_DATA = join(PATH_BASE, "data")
PATH_FUNCTIONS = join(PATH_BASE, "src", "python", "functions.py")
KEYWORDS = ["#GE2017", "#GeneralElection", "#Tory", "#SNP", "#LibDem", "#Con", "#Lab"]
KEYWORDS = " OR ".join(KEYWORDS)

def hausdorff_distance_km(polygon):
    """
    Calculate hausdorff distance of polygons centroid to its border

    Parameters
    ----------
    polygon : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    exterior = polygon.exterior
    centroid = polygon.centroid
    distances_to_centroid = list()
    for p in list(exterior.coords):
        point = shapely.geometry.Point(p)
        distance = point.distance(centroid)
        distances_to_centroid.append(distance)
        
    idx = np.argmax(distances_to_centroid)
    hausdorff_point = shapely.geometry.Point(list(exterior.coords)[idx])
    distance = geopy.distance.GeodesicDistance(hausdorff_point.coords, centroid.coords).kilometers
    return(distance)
###############################################################################
GE_2017 = datetime.date(2017, 6, 8)
GE_2010 = datetime.date(2010, 5, 6)
GE_2015 = datetime.date(2015, 5, 7)

series_2017 = pd.date_range(end = GE_2017, periods = 365)
series_2010 = pd.date_range(end = GE_2010, periods = 365)
series_2015 = pd.date_range(end = GE_2015, periods = 365)

constituencies = gpd.read_file(join(PATH_DATA,"Westminster_Parliamentary_Constituencies_(December_2019)_Boundaries_UK_BUC.geojson"))

for day_end in series_2017:
    day_start = day_end + relativedelta(days = -1)
    day_start = day_start.strftime("%Y-%m-%d")
    day_end = day_end.strftime("%Y-%m-%d")
    ##
    for idx, constituency in constituencies.iterrows():
        name = constituency["pcon19nm"]
        geom = constituency["geometry"]
        if isinstance(geom, shapely.geometry.polygon.Polygon):
            centroid = geom.centroid
            max_radius = hausdorff_distance_km(geom)
            geo_loc = "" + str(centroid.coords[0][0]) + ", " + str(centroid.coords[0][1]) + ", " + str(max_radius) + "km"
            print(name + "\n" + geo_loc + "\n")
    
# query
c = twint.Config()
c.since = day_start
c.until = day_end
c.Pandas = True
c.Geo = geo_loc
c.Limit = 500000
c.Search = KEYWORDS    
twint.run.Search(c)
