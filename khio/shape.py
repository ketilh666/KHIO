# -*- coding: utf-8 -*-
"""
Functions to write files for export to Petrel or Matlab

    write_xyz: Write co-ordinates and data to xyz file
    writ_irap_grid: Write co-ordinates and data to IRAP ascii grid file

Functions to read files for import from Petrel or Matlab

For many, khio is short for Kunsthøgskolen i Oslo

Created on Fri Dec  4 22:54:00 2020
@author: kehok@equinor.com
"""

import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt

import geopandas as gpd
import shapely.geometry as geo

#------------------------------------------------------
# Write Pandas df to shapefile
#------------------------------------------------------

def df_to_shp(shapefile, df, key_x='lon', key_y='lat', **kwargs):
    """Write data from pd.DataFrame to shapefile
    
    Parameters
    ----------
    shapefile: str. Name of shapefile
    df: pd.DataFrame
    key_x, key_y: str. Keys for coordinates, 
        optional, defaults are key_x='lon', key_y='lat'
    
    kwargs
    ------
    junk: list. List of keys to be dropped before output, default is junk=[]
    crs: str. Reference system, default is crs='EPSG:4326'

    Returns
    -------
    ierr. int
    
    Programmed: KetilH, 25. October 2023
    """

    ierr = 0
    
    # Get the kwargs
    junk = kwargs.get('junk', [])
    crs  = kwargs.get('crs', 'EPSG:4326')
    verbose = kwargs.get('verbose', 0)
    
    # check the file nae
    if shapefile.split('.')[-1] != 'shp':
        shapefile = shapefile + '.shp'
    
    if verbose>0:
        print('khio.df_to_shp: Write df to shapefile')
        print(f' o file = {shapefile}')
        print(f' o crs  = {crs}')
        print(f' o junk = {junk}')
        print(f' o key_x = {key_x}')
        print(f' o key_y = {key_y}')

    # Create shapely geometry object
    geometry = [geo.Point(xy) for xy in zip(df[key_x], df[key_y])]
    df = df.drop([key_x, key_y] + junk, axis=1)
    
    # Write to shapefile
    gdf = gpd.GeoDataFrame(df, crs='EPSG:4326', geometry=geometry)
    gdf.to_file(shapefile, driver='ESRI Shapefile')

    return ierr

