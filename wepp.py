from __future__ import print_function

# Copyright (c) 2014, Roger Lew (rogerlew.gmail.com)
# All rights reserved.
#
# The project described was supported by NSF award number IIA-1301792 
# from the NSF Idaho EPSCoR Program and by the National Science Foundation.

"""
This script produces Keyhole Markup Language (KML) visualizations of
the Water Erosion and Prediction Project (WEPP) model outputs
"""
__version__ = "0.0.1"

# standard library modules
import csv
import datetime
import os
import math
import multiprocessing
import sys
import subprocess
import time
import shutil
import zipfile

from copy import copy
from collections import namedtuple
from glob import glob
#from xml.dom import minidom

# 3rd party dependencies
import shapefile
import h5py
import numpy as np
from numpy.testing import assert_array_equal
import matplotlib.pyplot as plt
from matplotlib.colors import rgb_to_hsv, hsv_to_rgb
import pyproj
from shapely import ops
from shapely.geometry import Polygon, MultiPolygon, LinearRing

from matplotlib.colors import Normalize, ListedColormap, \
                              LinearSegmentedColormap

from osgeo import gdal, ogr, osr
from osgeo.gdalconst import *

from weppout import Wat

georefTiff = r'C:\Python27x64\Lib\site-packages\WEPP_outputs\subcatchment\Fernan_16bit_wBathymetry3.tif'
hillshade = r'inputs\hillshade.tif'

#
# WEPP specific data container class and functions used to pull data from kml
#

class Watershed:
    def __init__(self, topaz_wepp_map, shp, prj, dbf, hdf5_fn=None):
        
    
        # Before WEPP another routine is used to identify the hillslopes from the 
        # DEM. This program is Topaz. Topaz and WEPP both assign integer IDs to the 
        # hillslopes but they are different. The WEPP IDs are consecutive from 1.
        # the Topaz IDs are non-consecutive. The subcatchment's in the shape file
        # produced by the WEPP analysis codes the Topaz IDs in the attribute table
        # under the "DN" field name. Here we are reading a space delimited table to
        # get the mapping betwen Topaz IDs and Wepp IDs
    
        with open(topaz_wepp_map) as fid:
            rdr = csv.DictReader(fid, delimiter=' ')
            t2w = dict([(int(L['TOPAZ']), int(L['WEPP'])) for L in rdr])
            w2t = dict([(v,k) for (k,v) in t2w.items()])

            # this assures one-to-one mapping
            assert len(list(t2w.keys())) == len(list(w2t.keys()))

        self.hdf5 = None
        if hdf5_fn is not None and os.path.exists(hdf5_fn):
            self.hdf5 = h5py.File(hdf5_fn)

        # subcatchments.shp polygon contains projected coordinates. We want to
        # obtain geographic coordinates. This can be accomplished with pyproj
        prjText = open( prj, 'r').read()
        srs = osr.SpatialReference()
        srs.ImportFromWkt( prjText )
        inProj = pyproj.Proj( srs.ExportToProj4() )
        outProj = pyproj.Proj(init='epsg:4326')
        transform = lambda x, y : pyproj.transform(inProj, outProj, x, y)

        # use shapefile (pyshp) to unpack shapes from shp 
        subcatchments = []
        channels = []

        sf = shapefile.Reader(shp, dbf)
        for shapeRecs in sf.shapeRecords():
            shp, rec = shapeRecs.shape, shapeRecs.record
            x, y = zip(*shp.points)
            lng, lat = transform(x, y)
            topaz = rec[0]
            wepp = t2w.get(topaz, None)
            if wepp:
                subcatchments.append(Subcatchment(topaz, wepp, list(zip(lng,lat)), self))
            else:
                channels.append(Channel(list(zip(lng,lat))))

        self.n = len(subcatchments)
        self.t2w = t2w
        self.w2t = w2t
        self.inProj = inProj
        self.outProj = outProj
        self.transform = transform
        self.subcatchments = subcatchments
        self.channels = channels

        self._findExtents()

#        self.buildMask()
        
    def __iter__(self):
        for sub in self.subcatchments:
            yield sub

    def __getitem__(self, index):
        return self.subcatchments[index]

    def getWeppID(topazID):
        return self.t2w[topazID]
    
    def getTopazID(weppID):
        return self.w2t[weppID]

    def _findExtents(self):
        mp = MultiPolygon([Polygon(sub.data) for sub in self.subcatchments])
        
        # order left, top, right, bottom for gdal 
        self.extents = [mp.bounds[0], mp.bounds[3], mp.bounds[2], mp.bounds[1]]

        mp = None
               
    def buildMask(self):
        global georefTiff
        indx_sub = findProjectionIndices(self.subcatchments)
        indx_chn = findProjectionIndices(self.channels)

        ds = gdal.Open(georefTiff)

        x = np.zeros((ds.RasterYSize, ds.RasterXSize), dtype=np.uint8)
        x[indx_sub] = 1
        x[indx_chn] = 1

        self.mask = x < 1
        
        indx_sub = indx_chn = x = ds = None

    def timeseries(self, outputtype, pv, weppids=None):
        fig = plt.figure()

        dates = []
        values = []
        if outputtype == 'wat':
            for sub in self:
                if weppids is None:
                    values.append(sub.wat.data[pv])
                elif sub.wepp in weppids:
                    values.append(sub.wat.data[pv])

            d = self.subcatchments[0].wat.data
            dates = [datetime.datetime(y,m,d) for y,m,d in zip (d['Y'], d['M'], d['D'])]

        dates = np.array(dates)
        values = np.array(values)
        stdev = np.std(values, axis=0).flatten()
        values =np.mean(values, axis=0).flatten()
        ci = stdev * 1.96
        ci_u = values + ci
        ci_l = values - ci

        plt.plot(dates, values)
        plt.fill_between(dates, ci_l, ci_u,  lw=0, alpha=0.4, facecolor='b')

        return fig
        
         
    
class Channel:
    """
    Represents a Channel
    """
    def __init__(self, data):
        self.data = data

        # will define these later
        # a list of indices in projected, georeferenced
        # raster that coorespond to the channel
        self.indices = None 

    def __str__(self):
        return 'dn: %i\nlen(data): %s'%(self.dn, self.n)

class Subcatchment:
    """
    Represents a Subcatchment
    """
    def __init__(self, topaz, wepp, data, watershed):
        self.topaz = topaz
        self.wepp = int(wepp)
        self.data = data
        self.n = len(data)
        
        # will define these later
        # a list of indices in projected, georeferenced
        # raster that coorespond to the hillslope
        self.indices = self.findProjectionIndices() 

        self.wat = None

        if watershed.hdf5 is not None:
            if 'wat' in watershed.hdf5.keys():
                watershed.hdf5['wat']
                self.wat = Wat()
                self.wat.read_hdf5(self.wepp, watershed.hdf5['wat'][str(self.wepp)])

    def __str__(self):
        return 'dn: %i\nlen(data): %s'%(self.dn, self.n)

    def centroid(self):
        centroid = Polygon(self.data).centroid
        return (centroid.x, centroid.y)
    
    def findProjectionIndices(self):
        global georefTiff

        # This function is based loosely off of Frank's tests for
        # gdal.RasterizeLayer.
        # https://svn.osgeo.org/gdal/trunk/autotest/alg/rasterize.py

        # open the reference 
        # we use this to find the size, projection,
        # spatial reference, and geotransform to
        # project the subcatchment to
        ds = gdal.Open(georefTiff)
    
        pszProjection = ds.GetProjectionRef()
        if pszProjection is not None:
            srs = osr.SpatialReference()
            if srs.ImportFromWkt(pszProjection ) == gdal.CE_None:
                pszPrettyWkt = srs.ExportToPrettyWkt(False)
            
        geoTransform = ds.GetGeoTransform()

        # initialize a new raster in memory
        driver = gdal.GetDriverByName('MEM')
        target_ds = driver.Create('', 
                                  ds.RasterXSize, 
                                  ds.RasterYSize, 
                                  1, gdal.GDT_Byte )
        target_ds.SetGeoTransform(geoTransform)
        target_ds.SetProjection(pszProjection)
    
        # close the reference
        ds = None

        # Create a memory layer to rasterize from.
        rast_ogr_ds = ogr.GetDriverByName('Memory') \
                            .CreateDataSource( 'wrk' )
        rast_mem_lyr = rast_ogr_ds.CreateLayer( 'poly', srs=srs )

        # Add a polygon.
        coords = ','.join(['%f %f'%(lng,lat) for lng,lat in self.data])
        wkt_geom = 'POLYGON((%s))'%coords
        feat = ogr.Feature( rast_mem_lyr.GetLayerDefn() )
        feat.SetGeometryDirectly( ogr.Geometry(wkt = wkt_geom) )
        rast_mem_lyr.CreateFeature( feat )

        # Run the rasterization algorithm
        err = gdal.RasterizeLayer( target_ds, [1], rast_mem_lyr,
                                    burn_values = [255] )
        rast_ogr_ds = None
        rast_mem_lyr = None

        band = target_ds.GetRasterBand(1)
        data = band.ReadAsArray()
    
        # find nonzero indices and return
        return np.nonzero(data)


## this is outside of Watershed to make it easier to parallelize
#def findProjectionIndices(subcatchments):
#    global georefTiff
#    if isinstance(subcatchments, Subcatchment):
#        subcatchments = [subcatchments]

#    # This function is based loosely off of Frank's tests for
#    # gdal.RasterizeLayer.
#    # https://svn.osgeo.org/gdal/trunk/autotest/alg/rasterize.py


#    # open the reference 
#    # we use this to find the size, projection,
#    # spatial reference, and geotransform to
#    # project the subcatchment to
#    ds = gdal.Open(georefTiff)
    
#    pszProjection = ds.GetProjectionRef()
#    if pszProjection is not None:
#        srs = osr.SpatialReference()
#        if srs.ImportFromWkt(pszProjection ) == gdal.CE_None:
#            pszPrettyWkt = srs.ExportToPrettyWkt(False)
            
#    geoTransform = ds.GetGeoTransform()

#    # initialize a new raster in memory
#    driver = gdal.GetDriverByName('MEM')
#    target_ds = driver.Create('', 
#                              ds.RasterXSize, 
#                              ds.RasterYSize, 
#                              1, gdal.GDT_Byte )
#    target_ds.SetGeoTransform(geoTransform)
#    target_ds.SetProjection(pszProjection)
    
#    # close the reference
#    ds = None

#    for sub in subcatchments:
#        # Create a memory layer to rasterize from.
#        rast_ogr_ds = ogr.GetDriverByName('Memory') \
#                         .CreateDataSource( 'wrk' )
#        rast_mem_lyr = rast_ogr_ds.CreateLayer( 'poly', srs=srs )

#        # Add a polygon.
#        coords = ','.join(['%f %f'%(lng,lat) for lng,lat in sub.data])
#        wkt_geom = 'POLYGON((%s))'%coords
#        feat = ogr.Feature( rast_mem_lyr.GetLayerDefn() )
#        feat.SetGeometryDirectly( ogr.Geometry(wkt = wkt_geom) )
#        rast_mem_lyr.CreateFeature( feat )

#        # Run the rasterization algorithm
#        err = gdal.RasterizeLayer( target_ds, [1], rast_mem_lyr,
#                                   burn_values = [255] )
#        rast_ogr_ds = None
#        rast_mem_lyr = None

#    band = target_ds.GetRasterBand(1)
#    data = band.ReadAsArray()
    
#    # find nonzero indices and return
#    return np.nonzero(data)

#def makeChoropleth(pv, ws, x, epoch, builddir, shading=True, kmloverlay=True):
#    global hillshade

#    print(hillshade)
        
#    print('    Building %s' % str(epoch))
#    i0 = epoch.i0
#    iend = epoch.iend
#    begin = epoch.begin
#    end = epoch.end

#    dst_fname = r'%s/%s.tif'%(builddir, begin)

#    # open the reference 
#    # we use this to find the size, projection,
#    # spatial reference, and geotransform to
#    # project the subcatchment to
#    hill_ds = gdal.Open(hillshade)
#    assert hill_ds is not None

#    hillband = hill_ds.GetRasterBand(1)
#    hillbandnodatavalue = hillband.GetNoDataValue()
#    bytes = hillband.DataType
    
#    assert bytes <= 2

#    hillbandmax = 2.0 ** (bytes*8.0) - 1.0
    
#    pszProjection = hill_ds.GetProjectionRef()
#    if pszProjection is not None:
#        srs = osr.SpatialReference()
#        if srs.ImportFromWkt(pszProjection ) == gdal.CE_None:
#            pszPrettyWkt = srs.ExportToPrettyWkt(False)
            
#    geoTransform = hill_ds.GetGeoTransform()

#    # initialize a new raster in memory
#    driver = gdal.GetDriverByName('GTiff')
#    dst_ds = driver.Create(dst_fname, 
#                             hill_ds.RasterXSize, 
#                             hill_ds.RasterYSize, 
#                             4, gdal.GDT_Byte )
#    dst_ds.SetGeoTransform(geoTransform)
#    dst_ds.SetProjection(pszProjection)

#    # find rgb values for new raster based on subcatchment projection
#    # indices and matplotlib colormap in pv
#    data = np.zeros((hill_ds.RasterYSize, hill_ds.RasterXSize, 3), dtype=np.float)
#    for sub in ws.subcatchments:
#        xInd, yInd = sub.indices
#        val, color = pv.interpret(x[sub.wepp-1, i0:iend])
#        r,g,b,a = color # color channels are scaled [0-1]
        
#        data[xInd, yInd, 0] = r
#        data[xInd, yInd, 1] = g
#        data[xInd, yInd, 2] = b

#    if shading:
#        hill_data = np.array(hillband.ReadAsArray(), np.float) 
#        hill_data /= hillbandmax # scale between 0-1

#        # bake hillside shading into raster
#        #
#        # We know these rasters will be pretty small so we
#        # don't have to convert them line by line like
#        # hsv_merge.py
#        hsv = rgb_to_hsv(data)
         
#        #replace v with hillshade
#        hsv[:,:,2] = hill_data

#        #convert back to 8-bit per channel RGB
#        data = hsv_to_rgb(hsv)
        
#    # Write color data to Gtiff
#    data = np.array(data*255.0, dtype=np.uint8)

#    for i in xrange(3):
#        band = dst_ds.GetRasterBand(i+1)
#        band.WriteArray(data[:,:,i])

#    # Add alpha channel
#    alpha = np.ones((hill_ds.RasterYSize, hill_ds.RasterXSize), dtype=np.uint8)
#    alpha *= np.uint8(alpha*255.0)
    
#    # Clip to just the watershed
#    alpha[ws.mask] = 0
#    band = dst_ds.GetRasterBand(4)
#    band.WriteArray(alpha)
    
#    # close the raster datasets
#    hill_ds = None
#    dst_ds = None

#    # Build kml overlay
#    if kmloverlay:
#        g2t = GDAL2Tiles( ['-k', '-z', '10-15',
#                           dst_fname, 
#                           r'%s/epoch%03i-%03i'%(builddir, i0, iend) ])
#        g2t.process()

#        time.sleep(1)
##
## KML templates for writing file. 
## The structure is pretty basic so simple is better
## 
   
#opentag = '''\
#<?xml version="1.0" encoding="utf-8"?>
#    <kml xmlns="http://www.opengis.net/kml/2.2">
#      <Document>
#        <name>{name}</name>
#        <description></description>
#        <Style>
#          <ListStyle id="hideChildren">
#            <listItemType>checkHideChildren</listItemType>
#          </ListStyle>
#        </Style>
#'''

#networklink = '''\
#        <NetworkLink>
#          <TimeSpan>
#            <begin>{begin}</begin>
#            <end>{end}</end>
#          </TimeSpan>
#          <name>{epoch}/11/359/1333.png</name>
#          <Region>
#            <LatLonAltBox>
#              <north>47.75409797968003</north>
#              <south>47.63578359086485</south>
#              <east>-116.71874999999999</east>
#              <west>-116.89453125000000</west>
#            </LatLonAltBox>
#            <Lod>
#              <minLodPixels>128</minLodPixels>
#              <maxLodPixels>-1</maxLodPixels>
#            </Lod>
#          </Region>
#          <Link>
#            <href>11/359/1333.kml</href>
#            <viewRefreshMode>onRegion</viewRefreshMode>
#            <viewFormat/>
#          </Link>
#        </NetworkLink>
#        <NetworkLink>
#          <TimeSpan>
#            <begin>{begin}</begin>
#            <end>{end}</end>
#          </TimeSpan>
#          <name>11/360/1333.png</name>
#          <Region>
#            <LatLonAltBox>
#              <north>47.75409797968003</north>
#              <south>47.63578359086485</south>
#              <east>-116.54296875000000</east>
#              <west>-116.71874999999999</west>
#            </LatLonAltBox>
#            <Lod>
#              <minLodPixels>128</minLodPixels>
#              <maxLodPixels>-1</maxLodPixels>
#            </Lod>
#          </Region>
#          <Link>
#            <href>{epoch}/11/360/1333.kml</href>
#            <viewRefreshMode>onRegion</viewRefreshMode>
#            <viewFormat/>
#          </Link>
#        </NetworkLink>    
#'''

#closetag = '''\
#          </Document>
#    </kml>
#'''

#class PV:
#    """
#    represents the parameters for a process variable
#    """
#    def __init__(self, name, lim, outname, colormap,
#                 alpha=1.0, aggregator=None):
#        if not aggregator:
#            aggregator = np.mean

#        self.name = name
#        self._setlim(lim) # lim is implemented at as a property to ensure
#                          # range and the norm_func stay updated
                                   
#        self.outname = outname
#        self.colormap = colormap # blue, green, red
#                           # this is a from KML

#        cm = plt.get_cmap(colormap)

#        self._cm = cm
#        self.alpha = alpha
#        self.aggregator = aggregator
    
#    def _setlim(self, lim):
#        self._lim = lim # expecting [minvalue->0, maxvalue->maxbit]
#                        # in physical units e.g. kg/s

#        self._rng = lim[1] - lim[0] # range of value in its physical units
#        assert self._rng > 0

#        # when called normalized data into the [0.0, 1.0] interval
#        self.norm_func = Normalize(lim[0], lim[1], clip=True)

#    def _getlim(self):
#        return self._lim

#    lim = property(_getlim, _setlim)

#    def aggregate(self, values):
#        return self.aggregator(values)

#    def interpret(self, values):
#        norm_func = self.norm_func
#        cm = self._cm
#        alpha = self.alpha

#        val = self.aggregator(values)
        
#        # color bytes are ordered alpha, blue, green, red (abgr) for KML
#        rgba = list(cm(norm_func(val)))[:-1] + [alpha]
##        rgba = [int(round(f*255.0)) for f in rgba]
##        hex_str = ''.join([u'%02x'%int(round(f*255.0)) for f in abgr])

#        return val, rgba

#def timeseries_ensemble(x, pv, date0, epoch0, numepochs, fname, 
#                        alpha=0.01, dpi=None):
#    if not dpi:
#        dpi = 3840.0 / 16.0
        
#    # figure out the starting and ending indices for the 2nd axis of x
#    i0 = epoch0 * dt
#    iend = (epoch0 + numepochs) * dt

#    # create a figure
#    plt.figure(figsize=(16,9))
#    plt.subplots_adjust(left=0.03, bottom=0.03, right=0.99, top=0.96)

#    # loop over subcatchments and plot dv
#    for j in xrange(x.shape[0]):
#        y = np.arange(i0+1, iend+1, 1.0)
#        z = x[j, i0:iend].flatten()
#        assert y.shape == z.shape
#        plt.plot(y, z, 'b', alpha=alpha)

#    # adjust x-axis limits and set tick labels
#    plt.xlim([i0-dt, iend+dt])

#    # find first of months
#    dates = [(i, date0 + datetime.timedelta(i)) for i in xrange(i0+1, iend+1)]
#    dates = [[i, str(date)] for (i,date) in dates if str(date)[-2:] == '01']
#    [ticks, labels] = zip(*dates)
#    plt.xticks(ticks, labels)
    
#    # set title
#    plt.title(pv.name)

#    # save figure
#    plt.savefig(fname, dpi=dpi)
#    plt.close()
        
Epoch = namedtuple('Epoch', ['i0', 'iend'])

def makeChoroplethHdf5(hillshade, pv, ws, builddir, epoch=None, cmap=None, ymin=0, ymax=10, alpha=1, ffmeg_frame_num=None):
        
    if epoch is None:
        epoch = Epoch(0, -1)
    i0 = epoch.i0
    iend = epoch.iend

    if cmap is None:
        cmap = plt.get_cmap('viridis')
    
    norm_f = Normalize( ymin, ymax, clip=True)

    if ffmeg_frame_num is not None:
        dst_fname =  r'%s/%s_%05i.tif' % ( builddir, pv, ffmeg_frame_num)
    else:
        dst_fname = r'%s/%s_%05i-%05i.tif' % ( builddir, pv, i0, iend)

    # open the reference 
    # we use this to find the size, projection,
    # spatial reference, and geotransform to
    # project the subcatchment to
    hill_ds = gdal.Open(hillshade)
    assert hill_ds is not None

    hillband = hill_ds.GetRasterBand(1)
    hillbandnodatavalue = hillband.GetNoDataValue()
    bytes = hillband.DataType
    
    assert bytes <= 2

    hillbandmax = 2.0 ** (bytes*8.0) - 1.0
    
    pszProjection = hill_ds.GetProjectionRef()
    if pszProjection is not None:
        srs = osr.SpatialReference()
        if srs.ImportFromWkt(pszProjection ) == gdal.CE_None:
            pszPrettyWkt = srs.ExportToPrettyWkt(False)
            
    geoTransform = hill_ds.GetGeoTransform()

    # initialize a new raster in memory
    driver = gdal.GetDriverByName('GTiff')
    dst_ds = driver.Create(dst_fname, 
                             hill_ds.RasterXSize, 
                             hill_ds.RasterYSize, 
                             4, gdal.GDT_Byte )
    dst_ds.SetGeoTransform(geoTransform)
    dst_ds.SetProjection(pszProjection)

    # find rgb values for new raster based on subcatchment projection
    # indices and matplotlib colormap in pv
    data = np.zeros((hill_ds.RasterYSize, hill_ds.RasterXSize, 4), dtype=np.float)
    for sub in ws.subcatchments:
        xInd, yInd = sub.indices
        val = np.mean(sub.wat.data[pv][i0:iend, :])
        normed_val = norm_f(val)
        r,g,b,a = list(cmap(val))[:-1] + [alpha]

        data[xInd, yInd, 0] = r
        data[xInd, yInd, 1] = g
        data[xInd, yInd, 2] = b
        data[xInd, yInd, 3] = a

    # Write color data to Gtiff
    data = np.array(data*255.0, dtype=np.uint8)

    for i in xrange(4):
        band = dst_ds.GetRasterBand(i+1)
        band.WriteArray(data[:,:,i])

    # close the raster datasets
    hill_ds = None
    dst_ds = None

if __name__ == '__main__':

    sub_shp = r'C:\Python27x64\Lib\site-packages\WEPP_outputs\subcatchment\subcatchments.shp'
    sub_prj = r'C:\Python27x64\Lib\site-packages\WEPP_outputs\subcatchment\subcatchments.prj'
    sub_dbf = r'C:\Python27x64\Lib\site-packages\WEPP_outputs\subcatchment\subcatchments.dbf'
    hdf5_fn = r'C:\Python27x64\Lib\site-packages\WEPP_outputs\WEPPout.hdf5'

    topaz_wepp_map = r'C:\Python27x64\Lib\site-packages\WEPP_outputs\TopazWEPP\Fernan_ck_input_files_clearcuts.txt'
    

    # build watershed instance representing subcatchments
    ws = Watershed(topaz_wepp_map, 
                   shp=sub_shp, prj=sub_prj, dbf=sub_dbf,
                   hdf5_fn=hdf5_fn)
                   
    for sub in ws:
        print(sub.topaz, sub.wepp, sub.centroid())

    print(ws[0].wat.data.keys())
    
    fig = ws.timeseries('wat', 'latqcc (mm)')
    fig.savefig(r'C:\Python27x64\Lib\site-packages\WEPP_outputs\latqcc_timeseries.png')

    epochs = [Epoch(i, i+7) for i in xrange(0, 50*52*7, 7)]

    for i, epoch in enumerate(epochs):
        print('making frame %04i' % i)
        makeChoroplethHdf5(georefTiff, 'latqcc (mm)',  ws, 
                           r'C:\Python27x64\Lib\site-packages\WEPP_outputs', 
                           epoch, ffmeg_frame_num=i)

