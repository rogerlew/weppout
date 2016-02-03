from __future__ import print_function

# Copyright (c) 2014, Roger Lew (rogerlew.gmail.com)
# All rights reserved.
#
# The project described was supported by NSF award number IIA-1301792 
# from the NSF Idaho EPSCoR Program and by the National Science Foundation.

"""
This packs GPH (geopotential height) outputs from the Water Erosion and Prediction Project (WEPP) 
model into a HDF5 container.
"""

__version__ = "0.0.1"

# standard library modules
import array
import glob
import os
import sys
import time
import multiprocessing

import h5py
import numpy as np
from wat import read_wat 
# 3rd party dependencies
DAYS = 'Days In Simulation'

def parse_float(x):
    try:
        return float(x)
    except:
        return float('nan')

def read_grp(fname):
    """
    reads a gph text data file


    Parameters
    ----------
      fname : string
        path to the gph file

    Returns
    -------
      (meta, data) : tuple(dict, dict)
        It is returning this way to support multiprocessing
    """
    global DAYS
    uint_types = [DAYS,
                  'Current crop type',                          
                  'Current residue on ground type',             
                  'Previous residue on ground type',            
                  'Old residue on ground type',                 
                  'Current dead root type',                     
                  'Previous dead root type',                    
                  'Old dead root type']

    meta = {}
    data = None
    header = []

    meta['fname'] = fname

    basename = os.path.basename(fname)
    meta['id'] = ''.join([L for L in basename if L in '0123456789'])
    
    with open(fname, 'rb') as fid:
        for i, line in enumerate(fid.readlines()):
            line_as_list = line.strip().split()

            if len(line_as_list) == 0:
                continue

            elif line_as_list[0][0] == '#':
                continue

            elif line_as_list[0] == 'int':
                meta[line_as_list[1]] = int(line_as_list[2])
                    
            elif line_as_list[0] == 'float':
                meta[line_as_list[1]] = float(line_as_list[2])

            elif line_as_list[0] == 'char':
                continue

            elif line_as_list[0][0] == '{':
                cname = line.strip()[1:-1].replace(r'kg/m',      r'kg*m**-1')   \
                                          .replace(r'kg/m**2',   r'kg*m**-2')   \
                                          .replace(r'kg/m**3',   r'kg*m**-3')   \
                                          .replace(r'kg/m**4',   r'kg*m**-4')   \
                                          .replace(r'mm/hr',     r'mm*hr**-1')  \
                                          .replace(r'mm/h',      r'mm*hr**-1')  \
                                          .replace(r'm/day',     r'm*day**-1')  \
                                          .replace(r'g/cc',      r'g*cc**-1')   \
                                          .replace(r'kg-s/m**4', r'kg-s*m**-4') \
                                          .replace(r's/m',       r's*m**-1')    \
                                          .replace(r'Irrigation_volume_supplied/unit_area',
                                                   r'Irrigation_volume_supplied*unit_area**-1')
                header.append(cname)

            else:
                if len(header) == len(line_as_list):
                        
                    # if we are here and data == None we need to initialize the data dictionary
                    if data == None:
                        data = {}
                        for cname in header:
                            typecode = ('f', 'h')[any([cname==s for s in uint_types])]
                            data[cname] = array.array(typecode)

                    for (cname, string) in zip(header, line_as_list):
                        if any([cname==s for s in uint_types]):
                            value = int(string)
                        else:
                            value = parse_float(string)

                        #if cname == DAYS:

                        #    if value in set(data[DAYS]):
                        #        break

                        data[cname].append(value)

                else:
                    raise Exception('Failed to parse line %i, unexpected number of columns.'%(i+1))

    num_ofes = meta['NumOFEs']

    # pack the table data into numpy arrays
    for (cname, v) in data.items():
        dtype = (np.float32, np.int16)[any([cname==s for s in uint_types])]
        x = np.array(v, dtype=dtype)
        L = len(x) - 2 # this is to not include min and max lines at end of file
        days_in_sim = L / num_ofes
        data[cname] = np.reshape(x[:L], (days_in_sim, num_ofes))

        ##if cname == DAYS:
        ##    data[cname] = data[cname][:, 0]

    return (meta, data)

if __name__ == '__main__':
    path = r'E:\WSU_Puyallup\Research_Stuffs_05202010\Post-doc_Idaho\3_Lake_Fernan\WEPP_Online_GIS\Upstream\2_Simulations_845_hills\Undisturbed\Simulations\Single_OFE - Old'
    ignore = ['pw0_gph.txt']
#CAS
    ignore2 = ['pw0_wat.txt'] #Just gave a name for channel water balance to ignore --> check what is the usual name?
    
    # parallel worker pool
    numcpus = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numcpus)

    fid = h5py.File(r'C:\Python27x64\Lib\site-packages\WEPP_outputs\WEPPout.hdf5', 'w')
    
    #
    # Read GPH
    #
    fnames = []
    for fname in glob.glob(os.path.join(path, '*_gph.txt')):
        if not any([fn in fname for fn in ignore]):
            fnames.append(fname)

    for fname in glob.glob(os.path.join(path, 'grph_*.txt')):
        if not any([fn in fname for fn in ignore]):
            fnames.append(fname)
#            print fname
#            read_grp(fname)

    # this launches the batch processing of the grp files
    packed_tuples = pool.imap(read_grp, fnames)

    # this reads the files without multiprocessing
    ##packed_tuples = []
    ##for fname in fnames:
    ##    print(fname)
    ##    packed_tuples.append(read_grp(fname))

    # this provides feedback as the sets of files complete. Using imap
    # guarentees that the files are in the same order as fnames but
    # delays receiving feedback

    gph = fid.create_group("gph")

    for i, (meta, data) in enumerate(packed_tuples):
        print('    {:<43}{:10}'.format(fnames[i], meta['id']))
        sub = gph.create_group(meta['id'])

        for k,v in meta.items():
            sub.attrs.create(k, v)

        for cname, v in data.items():
            sub.create_dataset(cname, compression="gzip", compression_opts=9, data=v)
            
    del packed_tuples
#CAS
#
# Read water files
#
    fnames = []
    for fname in glob.glob(os.path.join(path, '*_wat.txt')):
        if not any([fn in fname for fn in ignore2]):
            fnames.append(fname)

    packed_tuples = pool.imap(read_wat, fnames)

    ##packed_tuples = []
    ##for fname in fnames:
    ##    print(fname)
    ##    packed_tuples.append(read_wat(fname))

    wat = fid.create_group("wat")

    for i, (meta, data) in enumerate(packed_tuples):
        print('     {:<43}{:10}'.format(fnames[i],meta['id']))
        sub = wat.create_group(meta['id'])

        for cname, v in data.items():
            if cname != "Date":
                sub.create_dataset(cname, compression="gzip", compression_opts=9, data=v)

    del packed_tuples



#CAS end
    #
    # Read Channel
    #


    fid.close()
