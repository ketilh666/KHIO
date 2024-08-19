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

#------------------------------------------------------
# Read/Write xyz ascii file
#------------------------------------------------------

def write_xyz(file, x, y, rd, **kwargs):
    """Write data to xyz ascii file 
    
    Parameters
    ----------
    file: str, file name 
    x: array of float, 
        x or longitude, x.shape is (nx) or (ny*nx) or (ny,nx)
    y: array of float, 
        y or latitude,  y.shape is (ny) or (ny*nx) or (ny,nx)
    rd: list of array of float 
        data to output, rd[:].shape is  (ny*nx) or (ny,nx) 
    skip_nan: bool (default True)
        Do not write nan data if skip_nan is True 
    
    **kwargs
    --------
    head: list of strings, optional
        One header for each output column, len(head) is 2 + len(rd)
    
    verbose: int, print shit?
        
    # Example: Write 2D maps
    write_xyz('file.txt',x,y,[data1, data2])
        
    # Example: Write 3D cube
    write_xyz('file.txt',x,y,[z, data])
    
    Returns
    -------
    ierr, int
    
    Programmed: KetilH, March   2015 (Matlab)
                KetilH, January 2019 (Matlab)
                KetilH, 4. December 2020
    """
    
    # Get kwargs
    verbose = kwargs.get('verbose', 0)
    head = kwargs.get('head', [])
    skip_nan = kwargs.get('skip_nan', True)
    
    # Make a header string
    s1, s2 = ' '*7, ' '*10
    if len(head)>0: 
        head = s1 + s2.join(head)
    
    # Data must be in a list
    if not isinstance(rd,list): rd = [rd]
    nd = len(rd)
    
    if verbose>0: print('khio.write_xyz: {}'.format(file))
    if verbose>1: 
        print(' o len(rd) = {}'.format(nd))
        print(' o rd[0].size = {}'.format(rd[0].size))
    
    # Array sizes:
    mx, my = np.prod(x.shape), np.prod(y.shape)
    mz = [np.prod(rd[jj].shape) for jj in range(nd)]
    
    # Check array sizes
    for jj in range(1,nd):
        if mz[jj] != mz[0]: 
            print('khio.write_xyz: Error mz[0]={}, mz[j]={}'.format(mz[0], mz[jj]))
            return 1
            
    # Make a list of equal-shape flat numpy arrays
    wrk = []
    if mz[0] == my == mx:
        if verbose>1: print (' o array coords')
        wrk.append(x.ravel())
        wrk.append(y.ravel())
        for jj in range(nd): wrk.append(rd[jj].ravel())
    elif mz[0] == my*mx:
        if verbose>1: print (' o vector coords')
        gx, gy = np.meshgrid(x, y)
        wrk.append(gx.ravel())
        wrk.append(gy.ravel())
        for jj in range(nd): wrk.append(rd[jj].ravel())
    else:
        print('khio.write_xyz: Error mx={}, my={}'.format(mx, my))
        return 1
    
    # Remove nans
    if skip_nan:
        if verbose>1: print(' o skip nans')
        if verbose>1: print(' o len(wrk)={}'.format(len(wrk)))
        ind = np.isfinite(wrk[2]) # First data column
        for jj in range(len(wrk)):
            wrk[jj] = wrk[jj][ind]
    
    #print('array sizes:')
    #for ii in range(len(wrk)): print(len(wrk[ii]))
    
    if np.max(np.abs(x)) > 1e3:
        # UTMs
        fmt = ['%15.2f', '%15.2f'] + ['%12.3f' for ii in range(nd)]
    else:   
        # lon/lat
        fmt = ['%15.8f', '%15.8f'] + ['%12.3f' for ii in range(nd)]

    # Filename or file object?
    if isinstance(file, str):  file = open(file, 'w')

    try:
        # Write a header?
        if len(head) > 0:
            np.savetxt(file, np.array([head]), fmt='%s')
        # Write data
        np.savetxt(file, np.vstack(wrk).T, fmt=fmt)
        ierr = 0
    except:
        print('khio.write_xyz: Error, no data written to file')
        ierr = 1

    file.close()
    
    return ierr

def read_xyz(file, **kwargs):
    """Read data from an xyz ascii file
    
    Parameters
    ----------
    file: str or file object
    
    kwargs
    ------
    undef: undefined value, default is undef=-999999
    nhl: int, number of header lines, default is 0
    
    Returns
    -------
    data: array of floats. Data part of the file.
    head: list of str. Headers lines    
    
    Programmed: KetilH, 29 November 2022
    """
    
    undef = kwargs.get('undef', '*')
    nhl = kwargs.get('nhl', 0)

    # Filename or file object?
    if isinstance(file, str):  file = open(file)
    
    # Read the data
    swrk = file.readlines()
    swrk = [ss.strip() for ss in swrk] # Remove newlines 
    file.close()
    
    # Pop the header lines, if any
    head = []
    for jh in range(nhl):
        head.append(swrk.pop(0))
        
    # Make np array from the list and replace 'undef' with nan
    wrk = np.array([ss.split() for ss in swrk]) # np.arrray of objects
    ind = wrk[:,2] == undef
    wrk[ind,2] = np.nan
    data = wrk.astype(float)

    return data, head    

#------------------------------------
#  Rotate (x,y) coordinates
#------------------------------------

def rotate_xy(x, y, x0, y0, phi):
    """Shift and rotate (x,y) coordinate system
    
    Parameters
    ----------
    file: str, file name 
    x: array of float, 
        x or longitude, x.shape is (nx)
    y: array of float, 
        y or latitude,  y.shape is (ny)
    x0: float 
        Shifted origin
    y0: float 
        Shifted origin
    phi: float
        Rotation angle in degrees
        phi>0: counterclockwise rotation of coordinate system
        phi<0: clockwise rotation of coordinate system
            
    Returns
    -------
    xrot, yrot: array of float
        Rotated (x,y) coordinates
    
    Programmed: Ketil Hokstad, 1. February 2021
    """

    # Rotation matrix
    d2r = np.pi/180.0
    axx = np.cos(d2r*phi)
    axy = np.sin(d2r*phi)
    
    # Rotation
    xrot =  axx*(x-x0) + axy*(y-y0)
    yrot = -axy*(x-x0) + axx*(y-y0)
    
    return xrot, yrot

#------------------------------------
#  test the functions
#------------------------------------

if __name__ == '__main__':

    file_1, file_2 = 'junk_xyz_1.txt', 'junk_xyz_2.txt'
    file_3, file_4 = 'junk_xyz_3.txt', 'junk_xyz_4.txt'
    file_5 = 'junk_xyz_5.txt'
    
    file_a, file_b = 'junk_grd_a.txt', 'junk_grd_b.txt'
    
    x = np.array([11, 12, 13, 14, 15], dtype=float)
    y = np.array([1,3,5], dtype=float)
    
    gx, gy = np.meshgrid(x,y)
    
    d = np.array([[101, 102, 103, 104, 105],
                  [201, 202, 203, 204, 205],
                  [301, 302, 303, 304, 305]], dtype=float)

    head  = ['x','y','data']
    head2 = ['x','y','data1','data2']

    ierr_1 = write_xyz(file_1, x, y, d, verbose=1)
    ierr_2 = write_xyz(file_2, x, y, [d], head=head, verbose=1)
    ierr_3 = write_xyz(file_3, gx, gy, [d], head=head, verbose=1)
    ierr_4 = write_xyz(file_4, gx.ravel(), gy.ravel(), d.ravel(), verbose=1)
    ierr_5 = write_xyz(file_5, x, y, [d, 2*d], head=head2, verbose=2)

    ierr_a = write_irap_grid(file_a, x, y, d, verbose=1) 
    ierr_b = write_irap_grid(file_b, gx, gy, d, verbose=2) 
    
    print('ierr_j={}'.format([ierr_1, ierr_2, ierr_3, ierr_4, ierr_5]))
    print('ierr_x={}'.format([ierr_a, ierr_b]))
    
    # Test the rotation
    d2r, r2d = np.pi/180.0, 180.0/np.pi # Radians to/from degrees
    theta = np.array([0, 30, 45, 60, 90], dtype=float)
    rr = 10.0
    x = rr*np.cos(d2r*theta)
    y = rr*np.sin(d2r*theta)
    xrot, yrot = rotate_xy(x, y, 0, 0, 30)
