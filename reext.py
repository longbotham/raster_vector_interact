##[Example scripts]=group
##Input_raster=raster
##Input_vector=vector
##Transform_vector_to_raster_CRS=boolean
##Output_table=output table

import os
from osgeo import gdal, ogr, osr
#from processing.core.TableWriter import TableWriter
#from processing.core.GeoAlgorithmExecutionException import \
#        GeoAlgorithmExecutionException
#from processing.tools.raster import *

Input_raster = "./data/extract_test.tif"
Input_vector = "./data/extract_test.shp"

def mapToPixel( mX, mY, geoTransform ):
  if geoTransform[ 2 ] + geoTransform[ 4 ] == 0:
    pX = ( mX - geoTransform[ 0 ] ) / geoTransform[ 1 ]
    pY = ( mY - geoTransform[ 3 ] ) / geoTransform[ 5 ]
  else:
    pX, pY = applyGeoTransform( mX, mY, invertGeoTransform( geoTransform ) )
  return int( pX + 0.5 ), int( pY + 0.5 )

def pixelToMap( pX, pY, geoTransform ):
  mX, mY = applyGeoTransform( pX, pY, geoTransform )
  return mX, mY

def applyGeoTransform( inX, inY, geoTransform ):
  outX = geoTransform[ 0 ] + inX * geoTransform[ 1 ] + inY * geoTransform[ 2 ]
  outY = geoTransform[ 3 ] + inX * geoTransform[ 4 ] + inY * geoTransform[ 5 ]
  return outX, outY

def invertGeoTransform( geoTransform ):
  det = geoTransform[ 1 ] * geoTransform[ 5 ] - geoTransform[ 2 ] * geoTransform[ 4 ]

  if abs( det ) < 0.000000000000001:
    return

  invDet = 1.0 / det

  # compute adjoint and divide by determinate
  outGeoTransform = [ 0, 0, 0, 0, 0, 0 ]
  outGeoTransform[ 1 ] = geoTransform[ 5 ] * invDet
  outGeoTransform[ 4 ] = -geoTransform[ 4 ] * invDet

  outGeoTransform[ 2 ] = -geoTransform[ 2 ] * invDet
  outGeoTransfrom[ 5 ] = geoTransform[ 1 ] * invDet

  outGeoTransform[ 0 ] = ( geoTransform[ 2 ] * geoTransform[ 3 ] - geoTransform[ 0 ] * geoTransform[ 5 ] ) * invDet
  outGeoTransform[ 3 ] = ( -geoTransform[ 1 ] * geoTransform[ 3 ] + geoTransform[ 0 ] * geoTransform[ 4 ] ) * invDet

  return outGeoTransform


raster = gdal.Open(Input_raster)

rasterBaseName = os.path.splitext(os.path.basename(Input_raster))[0]

bandCount = raster.RasterCount
rasterXSize = raster.RasterXSize
rasterYSize = raster.RasterYSize
geoTransform = raster.GetGeoTransform()
rasterCRS = osr.SpatialReference()
rasterCRS.ImportFromWkt(raster.GetProjectionRef())

vector = ogr.Open(Input_vector, False)
layer = vector.GetLayer(0)
featureCount = layer.GetFeatureCount()
if featureCount == 0:
    raise GeoAlgorithmExecutionException(
            'There are no features in input vector.')

vectorCRS = layer.GetSpatialRef()

columns = []
featureDefn = layer.GetLayerDefn()
for i in xrange(featureDefn.GetFieldCount()):
    fieldDefn = featureDefn.GetFieldDefn(i)
    columns.append([fieldDefn.GetNameRef()])

layer.ResetReading()
feature = layer.GetNextFeature()
while feature is not None:
    for i in xrange(featureDefn.GetFieldCount()):
        fieldDefn = featureDefn.GetFieldDefn(i)
        if fieldDefn.GetType() == ogr.OFTInteger:
            columns[i].append(feature.GetFieldAsInteger(i))
        elif fieldDefn.GetType() == ogr.OFTReal:
            columns[i].append(feature.GetFieldAsDouble(i))
        else:
            columns[i].append(feature.GetFieldAsString(i))
    feature = layer.GetNextFeature()

current = 0
total = bandCount + featureCount * bandCount

#if Transform_vector_to_raster_CRS:
#    coordTransform = osr.CoordinateTransformation(vectorCRS, rasterCRS)
#    if coordTransform is None:
#        raise GeoAlgorithmExecutionException(
#                'Error while creating coordinate transformation.')

columnName = rasterBaseName[:8]
for i in xrange(bandCount):
    current += 1
#    progress.setPercentage(int(current * total))

    rasterBand = raster.GetRasterBand(i + 1)
    try:
        data = rasterBand.ReadAsArray()
    except:
        raise GeoAlgorithmExecutionException(
                'Error reading raster data. File might be too big.')
    layer.ResetReading()
    feature = layer.GetNextFeature()
    col = []
    col.append(columnName + '_' + str(i + 1))
    while feature is not None:
        current += 1
#        progress.setPercentage(int(current * total))

        geometry = feature.GetGeometryRef()
        x = geometry.GetX()
        y = geometry.GetY()
#        if Transform_vector_to_raster_CRS:
#            pnt = coordTransform.TransformPoint(x, y, 0)
#            x = pnt[0]
#            y = pnt[1]
        (rX, rY) = mapToPixel(x, y, geoTransform)
        if rX >= rasterXSize or rY >= rasterYSize:
            feature = layer.GetNextFeature()
            continue
        value = data[rY, rX]
        col.append(value)

        feature = layer.GetNextFeature()

    rasterBand = None
    columns.append(col)

import pdb; pdb.set_trace()
raster = None
vector.Destroy()

writer = TableWriter(Output_table, 'utf-8', [])
row = []
for i in xrange(len(columns[0])):
    for col in columns:
        row.append(col[i])
    print col
    writer.addRecord(row)
    row[:] = []
