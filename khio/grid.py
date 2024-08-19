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

#---------------------------------------
#  Write IRAP ascii grid
#---------------------------------------

def write_irap_grid(file, x, y, grd, **kwargs):
    """Write data to IRAP ascii grid file 
    
    Parameters
    ----------
    file: str, file name 
    x: array of float, 
        x or longitude, x.shape is (nx)
    y: array of float, 
        y or latitude,  y.shape is (ny)
    grd: 2D array of float 
        data, grd.shape is (ny,nx) 
    
    **kwargs
    --------
    verbose: int, print shit?
        
    Returns
    -------
    ierr, int
    
    Programmed: Ketil Hokstad, Januuary 2015 (Matlab)
                Ketil Hokstad, June 2015 (Matlab)
                KetilH, 4. December 2020
    """

    # Constants
    NHL, LL = 19, 6 # header length, line lengths
    
    # Get kwargs:
    verbose = kwargs.get('verbose', 0)

    # Check if x, y are np.meshgrid arrays:
    if np.prod(x.shape) == np.prod(y.shape) == np.prod(grd.shape):
        x = x[0,:].flatten()
        y = y[:,0].flatten()
    
    # Header:
    head = np.zeros([NHL], dtype=float);
    head[ 0]  = -996         # ?
    head[ 1]  = grd.shape[0] # ny
    head[ 2]  = x[1]-x[0]    # dx
    head[ 3]  = y[1]-y[0]    # dy
    head[ 4]  = x[0]         
    head[ 5]  = x[-1]
    head[ 6]  = y[0]
    head[ 7]  = y[-1]
    head[ 8]  = grd.shape[1] # nx
    head[ 9]  = 0            # ?
    head[10] = x[0]
    head[11] = y[0]

    if verbose>0: print('khio.write_irap_grid: {}'.format(file))
    if verbose>1: print(' o ny = {:d}, nx = {:d}'.format(np.int(head[1]), np.int(head[8])))

    # Flatten and replace nans
    tmp = grd.flatten()
    tmp[np.isnan(tmp)] = 9999900      

    # Map the data to rows of 6 floats per row/line:
    n = np.prod(grd.shape)
    m, k = n // LL, n % LL    
    wrk0 = tmp[0:LL*m].reshape(m,LL)
    wrk1 = tmp[LL*m: ].reshape(1,k)

    if verbose>1: print(' o Got all the data? {}'.format(LL*m + k == n))

    # header formats
    hmt0 = ['%d', '%d', '%-14.6f', '%-14.6f']
    hmt1 = ['%-14.6f', '%-14.6f', '%-14.6f', '%-14.6f']
    hmt2 = ['%d', '%-9.6f', '%-14.6f', '%-14.6f']
    hmt3 = ['%d', '%d', '%1d', '%d', '%d', '%1d', '%d']

    lh0 = len(hmt0)
    lh1 = len(hmt1) + lh0
    lh2 = len(hmt2) + lh1
    lh3 = len(hmt3) + lh2

    # Data formats:
    fmt0 = ['%-14.6f' for jj in range(LL)] # 6 per line
    fmt1 = ['%-14.6f' for jj in range(k)]  # last line

    # Filename or file object?
    if isinstance(file, str):  file = open(file, 'w')

    try:
        # Write a header?
        np.savetxt(file, head[  0:lh0].reshape(1,len(hmt0)), fmt=hmt0)
        np.savetxt(file, head[lh0:lh1].reshape(1,len(hmt1)), fmt=hmt1)
        np.savetxt(file, head[lh1:lh2].reshape(1,len(hmt2)), fmt=hmt2)
        np.savetxt(file, head[lh2:lh3].reshape(1,len(hmt3)), fmt=hmt3)
        # Write data: 6 per line and rest
        np.savetxt(file, wrk0, fmt=fmt0)
        np.savetxt(file, wrk1, fmt=fmt1)
        ierr = 0
    except:
        print('khio.write_irap_frid: Error, no data written to file')
        ierr = 1

    file.close()
    
    return ierr

#------------------------------------
#   Read IRAP ASCII grid
#------------------------------------

def read_irap_grid(fname, **kwargs):
    """Read IRAP ASCII grid exported from Petrel 
    
    Parameters
    ----------
    fname: str. Filename
    
    kwargs
    ------
    verbose: in. Print some shit?
    ret_all: bool. Return x, y, grd, dd if ret_all=True
    
    Returns
    -------
    x, y, grd 
    
    
    Programmed: KEtilH, 20. November 2022
    
    """

    nhl = kwargs.get('nhl', 19)
    verbose = kwargs.get('verbose', 0)
    ret_all = kwargs.get('ret_all', False)

    file = open(fname)
        
    # Read header lines
    hd_list = [file.readline() for jj in range(4)]
    # Get rid of the '\n' at end of line
    hd_list = [hd.strip() for hd in hd_list]

    # Unpack lists
    shd = []
    for hd in hd_list:
        sw_list = hd.split()
        for sw in sw_list:
            shd.append(sw)
            
    dd = {
            'ny': int(shd[1]),
            'dx': float(shd[2]),
            'dy': float(shd[3]),
            'x1': float(shd[4]),
            'x2': float(shd[5]),
            'y1': float(shd[6]),
            'y2': float(shd[7]),
            'nx': int(shd[8]),
            'x0': float(shd[10]),
            'y0': float(shd[11])
        }

    if verbose>0:
        print('khio.read_irap_grid:')
        print(f'   o fname = {fname}')
        print('   o nx, ny = {}, {}'.format(dd['nx'], dd['ny']))
        print('   o dx, dy = {}, {}'.format(dd['dx'], dd['dy']))
        print('   o x0, y0 = {}, {}'.format(dd['x0'], dd['y0']))
                       
    # x and y coordinates
    x = np.linspace(dd['x1'], dd['x2'], dd['nx']) # Easting
    y = np.linspace(dd['y1'], dd['y2'], dd['ny']) # Northing  
            
    # Read the grid data
    ws_list = file.readlines()
    # Get rid of the '\n' at end of line
    ws_list = [ws.strip() for ws in ws_list]
    
    wd_list = []
    for ws in ws_list:
        wd = np.array(ws.split()).astype(float)
        wd_list.append(wd)
     
    # Make the grid:
    grd = np.concatenate(wd_list).reshape(dd['ny'], dd['nx'])

    # Put the nans in
    ind = grd >= 9999900.0
    grd[ind] = np.nan

    if ret_all:
        return x, y, grd, dd
    else:
        return x, y, grd

#------------------------------------
#  Write ESRI grid
#------------------------------------

def write_esri_grid(file, x, y, grd, **kwargs):
    """Write data to IRAP ascii grid file 
    
    Parameters
    ----------
    file: str, file name 
    x: array of float, 
        x or longitude, x.shape is (nx)
    y: array of float, 
        y or latitude,  y.shape is (ny)
    grd: 2D array of float 
        data, grd.shape is (ny,nx) 
    
    **kwargs
    --------
    verbose: int, print shit?
        
    Returns
    -------
    ierr, int
    
    Programmed: Ketil Hokstad, 21. August 2022
                KetilH, 4. December 2020
    """
    
    NODATA_value = kwargs.get('NODATA_value', -9999)

    # grd = gwrk[1540:1550, 0:20]
    # x = grd_pred.x[0:20]
    # y = grd_pred.x[1540:1550]

    nrows, ncols = grd.shape
    xllcorner, yllcorner = x[0], y[0]
    cellsize = x[1]-x[0]
        
    # Replace nans
    ind = np.isnan(grd)
    grd[ind] = NODATA_value
    
    # Filename or file object?
    if isinstance(file, str):  file = open(file, 'w')
       
    try:
        # Write the header
        file.write('nrows {:15d}\n'.format(nrows))    
        file.write('ncols {:15d}\n'.format(ncols))    
        file.write('xllcorner {:14.2f}\n'.format(xllcorner))    
        file.write('yllcorner {:14.2f}\n'.format(yllcorner))    
        file.write('cellsize {:15.2f}\n'.format(cellsize))    
        file.write('NODATA_value {:8d}\n'.format(NODATA_value))    
    
        # Write data
        np.savetxt(file, grd[::-1,:], fmt = '%-.2f')
        #np.savetxt(file, grd, fmt = '%-.3e')

        ierr = 0

    except:
        print('khio.write_esri_grid: Error, no data written to file')
        ierr = 1

    file.close()
    
    # Repair
    grd[ind] = np.nan
    
    return ierr

