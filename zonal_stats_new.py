"""
Zonal Statistics
Vector-Raster Analysis

Copyright 2013 Matthew Perry

Changes by AsgerPetersen:
    - Use band nodata value if user does not override it
    - Allow calculating statistics for all geometry types
    - Allow buffer to be applied before calculation
    - Allow user to select which measures to calculate
    - Write output to ogr dataset
    - argparse
    - Allow calculation of percentiles

usage: zonal_stats.py [-h] [--buffer distance] [--nodata value] [-f FORMAT]
                      [-nln name] [-dsco [DSCO [DSCO ...]]]
                      [-lco [LCO [LCO ...]]] [--preload]
                      [--stats [{min,max,avg,std,sum,cnt,med} [{min,max,avg,std,sum,cnt,med} ...]]]
                      [--hist bins] [--perc [percentiles [percentiles ...]]]
                      [--prefix PREFIX]
                      inputvector raster outputvector [layer]

positional arguments:
  inputvector           Vector source
  raster                Raster file
  outputvector          Output vector datasource. See OGR documentation
  layer                 Layer name from input data source. See OGR
                        documentation

optional arguments:
  -h, --help            show this help message and exit
  --buffer distance     Buffer geometry by this distance before calculation
  --nodata value        Use this nodata value instead of value from raster
  -f FORMAT, --format FORMAT
                        Output format name. See OGR documentation
  -nln name, --outlayer name
                        Output layer name. See OGR documentation
  -dsco [DSCO [DSCO ...]]
                        Datasource creation options. See OGR documentation
  -lco [LCO [LCO ...]]  Layer creation options. See OGR documentation
  --preload             Preload entire raster into memory instead of a read
                        per vector feature
  --stats [{min,max,avg,std,sum,cnt,med} [{min,max,avg,std,sum,cnt,med} ...]]
                        Measure to calculate for each feature
  --hist bins           Calculate histogram with number of bins
  --perc [percentiles [percentiles ...]]
                        Calculate percentiles
  --prefix PREFIX       Column prefix
"""
from osgeo import gdal, ogr
from osgeo.gdalconst import *
import numpy as np
import sys
import argparse
gdal.PushErrorHandler('CPLQuietErrorHandler')
ogr.UseExceptions()

def bbox_to_pixel_offsets(gt, bbox):
    originX = gt[0]
    originY = gt[3]
    pixel_width = gt[1]
    pixel_height = gt[5]
    x1 = int((bbox[0] - originX) / pixel_width)
    x2 = int((bbox[1] - originX) / pixel_width) + 1

    y1 = int((bbox[3] - originY) / pixel_height)
    y2 = int((bbox[2] - originY) / pixel_height) + 1

    xsize = x2 - x1
    ysize = y2 - y1
    return (x1, y1, xsize, ysize)

ALLSTATS = ['min', 'max', 'avg', 'std', 'sum', 'cnt', 'med']
def zonal_stats(input_vector, raster_path, output_vector,
                    input_layer = None,
                    output_driver = 'ESRI shapefile',
                    output_layer= '',
                    output_dsco = [],
                    output_lco = [],
                    nodata_value=None,
                    column_prefix = '',
                    global_src_extent=False,
                    buffer_distance=0,
                    stats=ALLSTATS,
                    histbins=0,
                    percentiles = []
                    ):
    rds = gdal.Open(raster_path, GA_ReadOnly)
    assert rds, "Could not open raster"
    rb = rds.GetRasterBand(1)
    rgt = rds.GetGeoTransform()

    if nodata_value:
        # Override with user specified nodata
        nodata_value = float(nodata_value)
        rb.SetNoDataValue(nodata_value)
    else:
        # Use nodata from band
        nodata_value = float(rb.GetNoDataValue())
    # Warn if nodata is NaN as this will not work with the mask (as NaN != NaN)
    #assert nodata_value == nodata_value, "Cannot handle NaN nodata value"

    if buffer_distance:
        buffer_distance = float(buffer_distance)

    stats = [s.lower() for s in stats]
    #print "Calculating the following measures: ", stats

    vds = ogr.Open(input_vector, GA_ReadOnly)  # TODO maybe open update if we want to write stats
    assert(vds)
    vlyr =  vds.GetLayerByName(input_layer) if input_layer else vds.GetLayer(0)
    vdefn = vlyr.GetLayerDefn()

    # Calculate (potentially buffered) vector layer extent
    vlyr_extent = vlyr.GetExtent()
    if buffer_distance:
        expand_by = [-buffer_distance, buffer_distance, -buffer_distance, buffer_distance]
        vlyr_extent = [a + b for a, b in zip(vlyr_extent, expand_by)]

    # Create output lyr
    dstdrv = ogr.GetDriverByName( output_driver )
    dstds = dstdrv.CreateDataSource( output_vector, output_dsco )
    dstlyr = dstds.CreateLayer(output_layer, vlyr.GetSpatialRef(), vdefn.GetGeomType(), output_lco )

    for fld_index in xrange( vdefn.GetFieldCount() ):
        src_fd = vdefn.GetFieldDefn( fld_index )
        fd = ogr.FieldDefn( src_fd.GetName(), src_fd.GetType() )
        fd.SetWidth( src_fd.GetWidth() )
        fd.SetPrecision( src_fd.GetPrecision() )
        dstlyr.CreateField( fd )

    for s in stats:
        fd = ogr.FieldDefn( column_prefix + s, ogr.OFTReal )
        dstlyr.CreateField( fd )

    for p in percentiles:
      fd = ogr.FieldDefn( column_prefix + 'p' + str(p), ogr.OFTReal )
      dstlyr.CreateField( fd )

    for binnum in xrange(histbins):
        fd = ogr.FieldDefn( column_prefix + 'b' + str(binnum), ogr.OFTReal )
        dstlyr.CreateField( fd )


    dstfdefn = dstlyr.GetLayerDefn()

    # create an in-memory numpy array of the source raster data
    # covering the whole extent of the vector layer
    if global_src_extent:
        # use global source extent
        # useful only when disk IO or raster scanning inefficiencies are your limiting factor
        # advantage: reads raster data in one pass
        # disadvantage: large vector extents may have big memory requirements
        src_offset = bbox_to_pixel_offsets(rgt, vlyr_extent)
        src_array = rb.ReadAsArray(*src_offset)

        # calculate new geotransform of the layer subset
        new_gt = (
            (rgt[0] + (src_offset[0] * rgt[1])),
            rgt[1],
            0.0,
            (rgt[3] + (src_offset[1] * rgt[5])),
            0.0,
            rgt[5]
        )

    mem_drv = ogr.GetDriverByName('Memory')
    driver = gdal.GetDriverByName('MEM')

    # Loop through vectors

    skippednulgeoms = False
    total = vlyr.GetFeatureCount(force = 0)
    vlyr.ResetReading()
    count = 0
    feat = vlyr.GetNextFeature()
    while feat is not None:
        count = count + 1
        if count % 100 == 0:
            sys.stdout.write("\r{0} of {1}".format(count, total))
            sys.stdout.flush()
        if feat.GetGeometryRef() is None:
            # Null geometry. Write to dst and continue
            if not skippednulgeoms:
                print "\nWarning: Skipping nullgeoms\n"
                skippednulgeoms = True
            feat = vlyr.GetNextFeature()
            continue
        mem_feat = feat.Clone()
        mem_type = mem_feat.GetGeometryRef().GetGeometryType()
        if buffer_distance:
            mem_type = ogr.wkbPolygon
            mem_feat.SetGeometryDirectly( mem_feat.GetGeometryRef().Buffer(buffer_distance) )

        if not global_src_extent:
            # use local source extent
            # fastest option when you have fast disks and well indexed raster (ie tiled Geotiff)
            # advantage: each feature uses the smallest raster chunk
            # disadvantage: lots of reads on the source raster
            src_offset = bbox_to_pixel_offsets(rgt, mem_feat.geometry().GetEnvelope())
            src_array = rb.ReadAsArray(*src_offset)

            # calculate new geotransform of the feature subset
            new_gt = (
                (rgt[0] + (src_offset[0] * rgt[1])),
                rgt[1],
                0.0,
                (rgt[3] + (src_offset[1] * rgt[5])),
                0.0,
                rgt[5]
            )

        # Create a temporary vector layer in memory
        mem_ds = mem_drv.CreateDataSource('out')
        mem_layer = mem_ds.CreateLayer('mem_lyr', None, mem_type)
        mem_layer.CreateFeature(mem_feat)

        # Rasterize it
        rvds = driver.Create('', src_offset[2], src_offset[3], 1, gdal.GDT_Byte)
        rvds.SetGeoTransform(new_gt)
        gdal.RasterizeLayer(rvds, [1], mem_layer, burn_values=[1])
        rv_array = rvds.ReadAsArray()

        # Mask the source data array with our current feature
        # we take the logical_not to flip 0<->1 to get the correct mask effect
        # we also mask out nodata values explictly
        masked = np.ma.MaskedArray(
            src_array,
            mask=np.logical_or(
                src_array == nodata_value,
                np.logical_not(rv_array)
            )
        )

        import ipdb; ipdb.set_trace()

        # Destination feature
        dstfeat = ogr.Feature( dstfdefn )
        dstfeat.SetFrom( feat )
        for s in stats:
            if s == 'min':
                dstfeat.SetField(column_prefix + s, float(masked.min()))
            elif s == 'avg':
                dstfeat.SetField(column_prefix + s, float(masked.mean()))
            elif s == 'max':
                dstfeat.SetField(column_prefix + s, float(masked.max()))
            elif s == 'std':
                dstfeat.SetField(column_prefix + s, float(masked.std()))
            elif s == 'sum':
                dstfeat.SetField(column_prefix + s, float(masked.sum()))
            elif s == 'cnt':
                dstfeat.SetField(column_prefix + s, float(masked.count()))
            elif s == 'med':
                dstfeat.SetField(column_prefix + s, float(np.ma.extras.median(masked)))
            else:
                raise Exception("Unknown stats option: " + s)
        if histbins or percentiles:
          # We need a compressed version
          compressed = masked.compressed()
          if histbins:
              hist = np.histogram( compressed, histbins )
              for binnum in xrange(histbins):
                  dstfeat.SetField(column_prefix + 'b' + str(binnum), float( hist[0][binnum] ))
          if percentiles and compressed.size:
            values = np.percentile( compressed, percentiles )
            for p,v in zip( percentiles, values):
              dstfeat.SetField(column_prefix + 'p' + str(p), float( v ))
        dstlyr.CreateFeature( dstfeat )

        rvds = None
        mem_ds = None
        feat = vlyr.GetNextFeature()

    vds = None
    rds = None

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("inputvector",
                        help="Vector source")
    parser.add_argument("raster",
                        help="Raster file")
    parser.add_argument("outputvector",
                        help="Output vector datasource. See OGR documentation")
    parser.add_argument("layer", nargs='?', default=None,
                        help="Layer name from input data source. See OGR documentation")
    parser.add_argument("--buffer", type=float, default=0, metavar="distance",
                        help="Buffer geometry by this distance before calculation")
    parser.add_argument("--nodata", type=float, default=None, metavar="value",
                        help="Use this nodata value instead of value from raster")
    parser.add_argument("-f", "--format", default="ESRI shapefile",
                        help="Output format name. See OGR documentation")
    parser.add_argument("-nln", "--outlayer", default='', metavar="name",
                        help="Output layer name. See OGR documentation")
    parser.add_argument("-dsco" , nargs='*', action='append', default= [],
                        help="Datasource creation options. See OGR documentation")
    parser.add_argument("-lco", nargs='*', action='append', default = [],
                        help="Layer creation options. See OGR documentation")
    parser.add_argument("--preload", action="store_true",
                        help="Preload entire raster into memory instead of a read per vector feature")
    parser.add_argument("--stats", nargs='*', type=str, default=ALLSTATS, choices=ALLSTATS,
                        help="Measure to calculate for each feature")
    parser.add_argument("--hist", type=int, default=0, metavar="bins",
                        help="Calculate histogram with number of bins")
    parser.add_argument("--perc", nargs='*', type=int, default=[], metavar="percentiles",
                        help="Calculate percentiles")
    parser.add_argument("--prefix", default='',
                        help="Column prefix")


    args = parser.parse_args()

    zonal_stats(args.inputvector, args.raster, args.outputvector,
                    input_layer = args.layer,
                    output_driver = args.format,
                    output_layer = args.outlayer,
                    output_dsco = args.dsco,
                    output_lco = args.lco,
                    nodata_value = args.nodata,
                    column_prefix = args.prefix,
                    global_src_extent = args.preload,
                    buffer_distance = args.buffer,
                    stats = args.stats,
                    histbins = args.hist,
                    percentiles = args.perc
                    )
