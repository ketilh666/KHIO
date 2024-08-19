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
import re
import matplotlib.pyplot as plt
import segyio

#------------------------------------------------------
#   SEGY stuff
#------------------------------------------------------

# Get trace headers from segy-file
def parse_trace_headers(segyfile, **kwargs):
    """Parse the segy file trace headers into a pandas dataframe.
    Column names are defined from segyio internal tracefield
    One row per trace
    
    Parameters
    ----------
    segyfile: filename (str) or file object (segyio.open(...))
    
    kwargs
    ------
    n_trace: integer. Number of traces to parse, default is f.tracecount
    
    Returns
    -------
    df: Pandas DataFrame
    
    """
    
    # Make it work for both segy filename (str) and file object
    if isinstance(segyfile, str):
        fid = segyio.open(segyfile, 'r', ignore_geometry=True)
    else:
        fid = segyfile

    n_trace = kwargs.get('n_trace', fid.tracecount)

    # Get all header keys
    headers = segyio.tracefield.keys
    # Initialize dataframe with trace id as index and headers as columns
    df = pd.DataFrame(index=range(1, n_trace + 1), columns=headers.keys())
    # Fill dataframe with all header values
    for kw, bl in headers.items():
        #print(kw, bl)
        df[kw] = fid.attributes(bl)[0:n_trace]

    if isinstance(segyfile, str):
        fid.close()

    return df

# Get ascii header from segy-file
def parse_text_header(segyfile):
    """Format segy text header and output to readable, clean dict
    
    Parameters
    ----------
    segyfile: filename (str) or file object (segyio.open(...))
    
    Returns
    -------
    hd: dict
    """
    
    # Make it work for both segy filename (str) and file object
    if isinstance(segyfile, str):
        fid = segyio.open(segyfile, 'r', ignore_geometry=True)
    else:
        fid = segyfile
    
    raw_header = segyio.tools.wrap(fid.text[0])
    # Cut on C*int pattern
    cut_header = re.split(r'C ', raw_header)[1::]
    # Remove end of line return
    try:
        text_header = [x.replace('\n', ' ') for x in cut_header]
        text_header[-1] = text_header[-1][:-2]
    except: pass
    # Format in dict
    clean_header = {}
    try:
        i = 1
        for item in text_header:
            key = "C" + str(i).rjust(2, '0')
            i += 1
            clean_header[key] = item
    except: pass

    if isinstance(segyfile, str):
        fid.close()

    return clean_header

#---------------------------------------------
#
#---------------------------------------------
    
def segy_empty_header(**kwargs):
    """Create an empty SEGY header with the appropriate SU header words
    
    Parameters
    ----------
    
    Returns
    -------
    hd: dictionary
    
    Programmed:
        KetilH, 3.November 2021
    """

    keys = ['dt', 'fldr', 'cdp', 'cdpt', 'offset', 'cdpx', 'cdpy', 
            'sx','sy', 'gx', 'gy', 'xline', 'iline',
            'gelev','selev', 'sdepth', 'scalel','scalco', 'tstat']
    
    return {key: None for key in keys}

def segy_header(dt, gx, gy, **kwargs):
    """Create a SEGY header with the appropriate header words
    
    Parameters
    ----------
    dt: float. Sample rate
    gx: array of floats, shape=ntrace. x-coordinates
    gy: array of floats, shape=ntrace. y-coordinates
    
    kwargs
    ------
    iline. array of ints. NB must be given
    xline. array of ints. NB must be given
    
    Returns
    -------
    hd: dictionary
    
    Programmed:
        KetilH, 3.November 2021
    """

    # Initialize the dict
    #hd = segy_empty_header()
    hd = {}
        
    # Samplling (in micro secs)
    hd['dt'] = np.int(1000*dt)

    # Scaling factor for co-ordinates
    # hd['scalco'] = kwargs.get('scalco', 0)
    # rscl = 10**(hd['scalco']) # Scaling for coordinates
    rscl = 1.0
    
    # Coordinate headers:
    hd['gx'] = np.array(rscl*gx, dtype=int)
    hd['gy'] = np.array(rscl*gy, dtype=int)

    # More coordinates
    hd['sx'] = np.array(rscl*kwargs.get('sx', gx), dtype=int)
    hd['sy'] = np.array(rscl*kwargs.get('sy', gy), dtype=int)
    hd['cdpx'] = np.array(rscl*kwargs.get('cdpx', gx), dtype=int)
    hd['cdpy'] = np.array(rscl*kwargs.get('cdpy', gy), dtype=int)

    # ilines and xlines    
    hd['iline'] = np.array(kwargs.get('iline'), dtype=int)
    hd['xline'] = np.array(kwargs.get('xline'), dtype=int)
    hd['cdp']   = np.array(kwargs.get('cdp' , hd['iline']), dtype=int)
    hd['fldr']  = np.array(kwargs.get('fldr', hd['xline']), dtype=int)

    hd['offset']  = np.array(kwargs.get('offset', np.zeros_like(gx)), dtype=int)
    hd['cdpt']    = np.array(kwargs.get('cdpt', np.zeros_like(gx)), dtype=int)
    
    hd['tstat'] = np.array(kwargs.get('tstat', 0), dtype=int)
    
    return hd

#---------------------------------------------
#   Read/write SEGY file
#---------------------------------------------

bytes_3d_default = {
    'cdpx': segyio.su.cdpx,
    'cdpy': segyio.su.cdpy,
    'iline': segyio.su.iline,
    'xline': segyio.su.xline 
    }


def segy_read(segyfile, **kwargs):
    """ Read data and headers from a segy file (as unstructured data)
    
    Parameters
    ----------
    segyfile: filename (str) or file object (segyio.open(...))
    
    kwargs
    ------
    frst: int. First trace to read (default is 0)
    last: int. Last trace to read (fefault is f.tracecount)
    step: int. Trace increment (default is 1)
    kplot: bool. PLot seismic section?
    verbose: int. Print shit if verbose>0
    key_x: str. SU header word to label horizontal axis if kplot=True
    amin: float. Min clip amplitude in plot if kplot=True
    amax: float. Max clip amplitude in plot if kplot=True
    perc: float. Clip percentile, used if amin, amax not given and kplot=True
    k_3d: bool. Return 3D headers (cdpx, cdpy, inline, xline) (defaut is k_3d=False)
    bytes_3d: dict. Byte locations of 3D headers
    
    Default locations for 3D headers are
    
        bytes_3d_default = {
            'cdpx': 181,
            'cdpy': 185,
            'iline': 189,
            'xline': 193 
            }
    
    Other headers can be added to the dict bytes_3d to get other specific headers.
    
    Returns
    -------
    dd: Dict with data and headers

    Programmed: KetilH, September 2015
                KetilH, January   2023
    """
    
    # Make it work for both segy filename (str) and file object
    if isinstance(segyfile, str):
        fid = segyio.open(segyfile, 'r', ignore_geometry=True)
    else:
        fid = segyfile

    # Traces to read
    frst = kwargs.get('frst', 0)
    last = kwargs.get('last', fid.tracecount-1)
    step = kwargs.get('step', 1)
    
    # Get 3D geometry?
    k_3d = kwargs.get('k_3d', False)
    bytes_3d = kwargs.get('bytes_3d', bytes_3d_default)
    
    # QC pLot the data?
    kplot = kwargs.get('kplot', False)
    verbose = kwargs.get('verbose', 0)
    key_x = kwargs.get('key_x', 'tracl')
    amin  = kwargs.get('amin', None)
    amax  = kwargs.get('amax', None)
    perc  = kwargs.get('perc', 99.9)

    # Initialize output dict
    dd = {}

    # Header of the first trace
    hd0 = fid.header[0]
    sclt = 1e0
    dd['dt'] = 1e-3*hd0[segyio.su.dt] # Convert us to ms
    dd['ns'] = hd0[segyio.su.ns]

    dd['t'] = sclt*fid.samples
    dd['t0'] = sclt*fid.samples[0]

    # Scaling factor for coordinates
    scalco = hd0[segyio.su.scalco]
    if   scalco < 0:  sclh = -1.0/scalco
    elif scalco > 0:  sclh =  1.0*scalco
    else:             sclh =  1.0
        
        
    scalel = hd0[segyio.su.scalel]
    if   scalel < 0: sclv = -1.0/scalel
    elif scalel > 0: sclv =  1.0*scalel
    else:            sclv =  1.0
        
    ns, dt, t0 = dd['ns'] , dd['dt'] , dd['t0'] 
    if verbose >0:
        print(f'segy_read: ns = {ns}, dt={dt}, t0={t0}')
        print(f'segy_read: scalco, scalel = {scalco}, {scalel}')
        
    dd['gx'] = sclh*np.array([hd[segyio.su.gx] for hd in fid.header[frst:last+1:step]])
    dd['gy'] = sclh*np.array([hd[segyio.su.gy] for hd in fid.header[frst:last+1:step]])
    dd['sx'] = sclh*np.array([hd[segyio.su.sx] for hd in fid.header[frst:last+1:step]])
    dd['sy'] = sclh*np.array([hd[segyio.su.sy] for hd in fid.header[frst:last+1:step]])
  
    dd['gelev'] = sclv*np.array([hd[segyio.su.gelev] for hd in fid.header[frst:last+1:step]])
    dd['selev'] = sclv*np.array([hd[segyio.su.selev] for hd in fid.header[frst:last+1:step]])
    dd['sdepth'] = sclv*np.array([hd[segyio.su.sdepth] for hd in fid.header[frst:last+1:step]])
    dd['swdep'] = sclv*np.array([hd[segyio.su.swdep] for hd in fid.header[frst:last+1:step]])
    dd['gwdep'] = sclv*np.array([hd[segyio.su.gwdep] for hd in fid.header[frst:last+1:step]])

    dd['cdp']    = np.array([hd[segyio.su.cdp]    for hd in fid.header[frst:last+1:step]])
    dd['offset'] = np.array([hd[segyio.su.offset] for hd in fid.header[frst:last+1:step]])

    dd['tracl'] = np.array([hd[segyio.su.tracl] for hd in fid.header[frst:last+1:step]])
    dd['tracr'] = np.array([hd[segyio.su.tracr] for hd in fid.header[frst:last+1:step]])
    dd['fldr']  = np.array([hd[segyio.su.fldr]  for hd in fid.header[frst:last+1:step]])
    dd['tracf'] = np.array([hd[segyio.su.tracf] for hd in fid.header[frst:last+1:step]])
  
    dd['trid'] = np.array([hd[segyio.su.trid] for hd in fid.header[frst:last+1:step]])

    # Trace data
    nd = int((last-frst)/step + 1)
    data = np.zeros((nd, dd['ns']))
    for jj in range(nd):
        #print(jj, frst+jj*step)
        data[jj,:] = fid.trace[frst+jj*step]
        
    dd['data'] = data
    #dd['data'] = data.T
    
    # Optional: Get 3D headers
    if k_3d:
        if verbose>0: print('bytes_3d:')
        for key in bytes_3d.keys():
            if verbose>0: print(f' o key, byte = {key}, {bytes_3d[key]}')
            dd[key] = np.array([hd[bytes_3d[key]] for hd in fid.header[frst:last+1:step]])
        
        dd['cdpx'] = sclh*dd['cdpx']
        dd['cdpy'] = sclh*dd['cdpy']
                    
    # Make it work for both segy filename (str) and file object
    if isinstance(segyfile, str):
        fid.close()

    if kplot:
                
        swrk = str(fid).split('\n')[0]
        if '\\' in swrk:
            title = swrk.split('\\')[-1][:-1]
        else:
            title = swrk.split('\\')[-1][:-1]
                
        perc = np.minimum(perc, 100)
        clip = np.percentile(dd['data'], perc)
        if not amin: amin = -clip
        if not amax: amax =  clip
            
        fig = plt.figure(figsize=(16, 6))
        xtnt = [dd[key_x][0], dd[key_x][-1], (dd['ns']-1)*dd['dt'], 0]
        plt.imshow(dd['data'].T, origin='upper', extent=xtnt, cmap='gray_r')
        plt.clim(amin, amax)
        plt.colorbar()
        plt.axis('auto')
        plt.ylabel('Two-way time [ms]')
        plt.xlabel(f'{key_x} [???]')    
        plt.title(title)
        
        if 'gwdep' in dd.keys(): plt.plot(dd[key_x], dd['gwdep'], 'r-')
                
        dd['fig'] = fig
                
    return dd


def segy_write(segyfile, seis, *args, **kwargs):
    """ Write a cube or line to SEGY file.
    
    Parameters
    ----------
    segyfile: str. SEG-Y filename
    seis: np.array with seismic traces, shape=(ntrace, ns), 
          or dict with data array and headers
    
    args
    ----
    hd: dict with headers if seis is data array only 
    
    kwargs
    ------
    scalco: Coordinate scalar (default is scalco=1)
    scalel: Depth scalar (default is scalel=1)
    verbose: int. Print shit if verbose>0
    
    Returns
    -------
    Nothing
    
    Data array must have shape = (ntrace, nsamp)
    Header words in hd dict must have shape=(ntrace)
    
    Header byte loations in output file are as defined by SU.
    
    Examples
    --------
    1. Data np.array and headers in dict:
        segy_write('filename.sgy', data, hd)

    2. Data and headers in a common dict    
        segy_write('filename.sgy', seis)
        
    3. Applying coordinate scalar:
        segy_write('filename.sgy', seis, scalco=100)
    
    Programmed:
        KetilH, 28. September 2021
        KetilH, 3. November 2021
        ketilH, 2. January 2023
    """

    scalco = int(kwargs.get('scalco', 1))
    scalel = int(kwargs.get('scalel', 1))
    verbose = kwargs.get('verbose', 0)
    
    if len(args)>0:
        data = seis
        hd = args[0]
    else:
        data = seis['data']
        hd = seis
    
    # Number of samples per trace
    ntrace = data.shape[0]
    nsamp  = data.shape[1]
    
    if verbose >0:
        print('segy_write: ntrace, nsamp = {}, {}'.format(ntrace, nsamp))
        print(' o segyfile = {}'.format(segyfile))
        print(' o data min, max = {}, {}'.format(np.min(data), np.max(data)))
        print(' o scalco, scalel = {}, {}'.format(scalco, scalel))

    # Scaling factor for coordinates
    if   scalco < 0:  sclh = -1.0/scalco
    elif scalco > 0:  sclh =  1.0*scalco
    else:             sclh =  1.0
        
    if   scalel < 0: sclv = -1.0/scalel
    elif scalel > 0: sclv =  1.0*scalel
    else:            sclv =  1.0

    # SU keywords output
    su_i_list = ['fldr', 'iline', 'xline', 'cdp', 'offset']
    su_h_list = ['cdpx', 'cdpy', 'sx', 'sy', 'gx', 'gy']
    su_v_list = ['gelev', 'selev', 'sdepth', 'swdep', 'gwdep']
    
    # Check for presence of headers
    su_list = su_i_list + su_h_list + su_v_list
    for key in su_list:
        if not key in hd.keys(): hd[key] = np.zeros(ntrace)

    hdi = {}
    for key in su_i_list: hdi[key] = hd[key].astype(int)
    for key in su_h_list: hdi[key] = (sclh*hd[key]).astype(int)
    for key in su_v_list: 
        print(key)
        hdi[key] = (sclv*hd[key]).astype(int)
        
    # Make a spec struct
    spec = segyio.spec()
    spec.ilines  = np.unique(hd['iline'])
    spec.xlines  = np.unique(hd['xline'])
    spec.samples = np.linspace(0,hd['dt']*(nsamp-1), nsamp, dtype=int)
    spec.sorting = 1
    spec.format  = 1
    
    # Data need to be contiguous float32 C-array
    data = np.ascontiguousarray(data, dtype=np.float32)

    with segyio.create(segyfile, spec) as fut:
        for kk in range(ntrace):
            #print(f'kk={kk}')
            # Set trace header
            fut.header[kk] = {
                    # Original stuff
                    segyio.su.tracl:  kk,
                    segyio.su.fldr:   hdi['fldr'][kk],
                    segyio.su.iline:  hdi['iline'][kk],
                    segyio.su.xline:  hdi['xline'][kk],
                    segyio.su.cdp:    hdi['cdp'][kk],
                    segyio.su.offset: hdi['offset'][kk],
                    segyio.su.cdpx:   hdi['cdpx'][kk],
                    segyio.su.cdpy:   hdi['cdpy'][kk],
                    segyio.su.sx:     hdi['sx'][kk],
                    segyio.su.sy:     hdi['sy'][kk],
                    segyio.su.gx:     hdi['gx'][kk],
                    segyio.su.gy:     hdi['gy'][kk],
                    
                    segyio.su.gelev:  hdi['gelev'][kk],
                    segyio.su.selev:  hdi['selev'][kk],
                    segyio.su.sdepth: hdi['sdepth'][kk],
                    segyio.su.swdep:  hdi['swdep'][kk],
                    segyio.su.gwdep:  hdi['gwdep'][kk],
                    
                    segyio.su.scalco: -scalco,
                    segyio.su.scalel: -scalel,
                    
                    # 3D stuff:
                    segyio.tracefield.keys['CDP_X']: hdi['cdpx'][kk],
                    segyio.tracefield.keys['CDP_Y']: hdi['cdpy'][kk],
                    segyio.tracefield.keys['INLINE_3D']: hdi['iline'][kk],
                    segyio.tracefield.keys['CROSSLINE_3D']: hdi['xline'][kk]
                    }
            # Trace data
            fut.trace[kk] = data[kk,:]
        
        #fut.bin.update(tsort=segyio.TraceSortingFormat.CROSSLINE_SORTING)
        
    return 0
