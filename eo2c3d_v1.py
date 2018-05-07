# *****************************************************************************
# NAME: eo2c3d.py
# DEV: J. Heath Harwood
# DATE: 23  March 2018
#
#
# DESCRIPTION:  Script to will convert the EO output file from IPASCO+
#               to a comma delimited text file for use in Correlator 3D
#
#
# DEPENDENCIES: OGR/GDAL, Python Standard Library, easygui, utm
#
# SOURCE(S):    - Modified from RCD30_SubsetRead_TBX.py for ArcToolbox
#               - Modified doConvex Hull - https://gis.stackexchange.com/questions/49741/python-ogr-convert-a-point-geometry-into-polygon
#				- http://gdal.org/python/osgeo.ogr-module.html
#				- https://gis.stackexchange.com/questions/82935/ogr-layer-intersection
#				- https://pcjericks.github.io/py-gdalogr-cookbook/vector_layers.html
# 				
#*****************************************************************************

'''
Change log:
2018-March-23 J. Heath Harwood, tool is operational

TODO:
 -

'''

###>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.IMPORT STATEMENTs.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>###

import ogr, os, sys, glob, easygui, csv, string, shapefile, utm, traceback,fnmatch
import numpy as np
import ogr, osr
from osgeo import ogr

###>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.CONSTANTS.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>###

# MGRS Hard Coded Path
#mgrsDir = r'D:\heath\eo2c3d\_mgrs_grids\NAD83_UTM'
mgrsDir = r'D:\2018_PR\processing\_camera\eo2c3d\_mgrs_grids\WGS84_GEO'
esriDriver = ogr.GetDriverByName("ESRI Shapefile")

###>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.FUNCTIONS.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>###


def eos2shp(eoFileList):
    """
        Takes a list of IPAS CO+ Exterior Orientation (EO) files and creates a comma delimited EO file for
        Correlator 3D and Agisoft Photoscan.  This function also create shapefiles for the EO files as well as a
        Tile Definition File (*.tdf) for clipping the RGBN mosiac to MGRS Grids

        :param Applanix style EO *.txt file that is space delimited with header (15 rows) and footer of 9 rows

        :returns: CSV comma delimited *.txt file, *.shp (and associated files), *tdf file

    """
    # Get list of EO Files
    eoFiles = []
    for inEoFiles in eoFileList:
        eoFile = inEoFiles.split('\\')[-1]
        eoFiles.append(eoFile)
    nEos = eoFiles.__len__()
    print '     Number of input EO Files:  ' + str(nEos)

    # Parse the EO Files in the list
    for eo in eoFiles:
        eoPath = os.getcwd()
        global eoName
        eoName = eo
        print eoName
        eoShp = eoPath + '\\' + eoName.strip('txt') + 'shp'
        eoPoly = eoPath + '\\' + eoName.split('.')[0] + '_eoMBGPoly' + '.shp'
        eoBufferShp = eoPath + '\\' + eoName.split('.')[0] + '_eoBuff500mPoly' + '.shp'
        eoPrj = eoPath + '\\' + eoName.strip('txt') + 'prj'
        eoTxt = eoPath + '\\' + eoName.strip('.txt') + '_for_c3d.txt'
        eoAPS = eoPath + '\\' + eoName.strip('.txt') + '_for_photoscan.txt'
        print '     RCD30 EO File:  ' + eo

        # Do a quick scan and grab the UTM Zone; this won't give us the latitude number but allows us to set the
        # Projection for the EO Shapefile
        with open(eo, 'r') as f:
            for line in f:
                if 'Reference System Information:' in line:
                    global getZone
                    getZone = line[34:36]
                    print getZone

        # Open file and read the data into an array (str,int,float,float,float,float,float,float,float,float)
        eoFile = open(eo, 'r')

        # Read columns into the array
        eoData = np.genfromtxt(eoFile, skip_header=15, skip_footer=9, dtype=None,
                                  usecols=(0, 1, 3, 4, 5, 6, 7, 8, 9, 10),
                                  names=['image', 'event_num', 'easting', 'northing', 'ellip', 'omega', 'phi', 'kappa',\
                                         'lat', 'lon'])

        # Store the size of the eoData array
        eoSize = np.size(eoData)

        # Read all individual columns of data into their own array
        getImage = eoData['image']
        getEventNum = eoData['event_num']
        getEast = eoData['easting']
        getNorth = eoData['northing']
        getEllip = eoData['ellip']
        getOmega = eoData['omega']
        getPhi = eoData['phi']
        getKappa = eoData['kappa']
        getLat = eoData['lat']
        getLon = eoData['lon']

        # Write a projection file so we can view the shp in ArcGIS or QGIS
        eoDataSource = esriDriver.CreateDataSource(eoShp)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4269)
        #srs.SetUTM(int(getZone),1)
        print srs

        # Create the EO shapefile layer
        eoShpLayer = eoDataSource.CreateLayer(eoName, srs, ogr.wkbPoint)

        # Create attribute fields for the shapefile
        nameField = ogr.FieldDefn("image", ogr.OFTString)
        nameField.SetWidth(100)
        eoShpLayer.CreateField(nameField)
        print nameField
        eoShpLayer.CreateField(ogr.FieldDefn("event_num", ogr.OFTInteger))
        eoShpLayer.CreateField(ogr.FieldDefn("easting"), ogr.OFTReal)
        eoShpLayer.CreateField(ogr.FieldDefn("northing"), ogr.OFTReal)
        eoShpLayer.CreateField(ogr.FieldDefn("ellip"), ogr.OFTReal)
        eoShpLayer.CreateField(ogr.FieldDefn("omega"), ogr.OFTReal)
        eoShpLayer.CreateField(ogr.FieldDefn("phi"), ogr.OFTReal)
        eoShpLayer.CreateField(ogr.FieldDefn("kappa"), ogr.OFTReal)
        eoShpLayer.CreateField(ogr.FieldDefn("lat"), ogr.OFTReal)
        eoShpLayer.CreateField(ogr.FieldDefn("lon"), ogr.OFTReal)

        # Loop through the eoData Array and set the point fields for the shapefile
        for pts in range(0, eoSize):
            # Create the feature
            feature = ogr.Feature(eoShpLayer.GetLayerDefn())

            # create the WKT for the feature using Python string formatting
            wkt = "POINT(%f %f)" % (float(getLon[pts]), float(getLat[pts]))

            # Create the point from the Well Known Txt
            point = ogr.CreateGeometryFromWkt(wkt)

            # Set the feature geometry using the point
            feature.SetGeometry(point)
            feature.SetField("image", getImage[pts])
            feature.SetField("event_num", getEventNum[pts])
            feature.SetField("easting", getEast[pts])
            feature.SetField("northing", getNorth[pts])
            feature.SetField("ellip", getEllip[pts])
            feature.SetField("omega", getOmega[pts])
            feature.SetField("phi", getPhi[pts])
            feature.SetField("kappa", getKappa[pts])
            feature.SetField("lat", getLat[pts])
            feature.SetField("lon", getLon[pts])

            # Create the feature in the layer (shapefile)
            eoShpLayer.CreateFeature(feature)
            # Dereference the feature
            feature = None

        # Save and close the data source
        eoDataSource = None
        print '     Completed importing ' + eo + '.\n\n'

        # Lets figure out the UTM Zone to determine which grids we need
        boundingBox = getBoundary(eoShp)
        print boundingBox
        lowerLat = boundingBox[2]
        lowerLon = boundingBox[1]
        upperlat = boundingBox[3]
        upperLon = boundingBox[1]

        # Get the UTM Zones
        utmZones = getZones(lowerLat, lowerLon, upperlat, upperLon)
        # returns 'These are your utm zones. ['18', 'T', '18', 'S']' as an example
        print 'These are your utm zones.' + str(utmZones)

        # Find the grid files for clipping
        gridFiles = getGrids(mgrsDir, utmZones)
        print gridFiles

        # Buffer Distance must be in degree units
        bufferDist = 0.005000
        # Create buffer from EO points and MBG Convex Hull Geometry
        createBuffer(eoShp, eoBufferShp, bufferDist, eoPoly)

        # Loop through grids and clip those grids to the eo buffer file
        for grid in gridFiles:
            # Clip buffer to MGRS Grids
            print "These are the grid files " + str(gridFiles)
            print "This grid file in list " + grid
            curGrid = mgrsDir + '\\' + grid
            print "This is the current grid " + curGrid
            eoGridShp = eoPath + '\\' + eoName.split('.')[0] + '_' + grid[0:12] + '.shp'
            zone = grid[9:11]
            print curGrid
            #print eoGridShp
            clip2shp(eoBufferShp, zone, curGrid, eoGridShp)
            shp2tdf(eoGridShp)

        # Export the eoData to eoTXT for c3d
        fieldNamesC3D = ["image", "event_num", "easting", "northing", "ellip", "omega", "phi", "kappa", "lat", "lon"]
        with open(eoTxt, 'wb') as f:
            dw = csv.DictWriter(f, fieldNamesC3D)
            for row in eoData:
                dw.writerow(dict(zip(fieldNamesC3D, row)))

        print "Done writing the formatted EO file for Simactive Correlator 3D"

        # Export the eoData to eoAPS for Agisoft Photoscan
        with open(eoAPS, 'wb') as f:
            for row in eoData:
                image = str(row[0])
                lat = float(row[8])
                lon = float(row[9])
                ellipHt = float(row[4])
                omega = float(row[5])
                phi = float(row[6])
                kappa = float(row[7])
                dataAPS = ('%s,%s,%s,%s,%s,%s,%s \n') % (image,lat,lon,ellipHt,omega,phi,kappa)
                f.write(dataAPS)

        # Close the file
        f.close()
        print "Done writing the formatted EO file for Agisoft Photoscan"

def getZones(lowLat,lowLon,upLat,upLon):

    zones = []
    # Convert the lat and lon to UTM and determine zone
    upperUTMZone = utm.from_latlon(upLat, upLon)
    lowerUTMZone = utm.from_latlon(lowLat, lowLon)

    # Feed the zone of the upper and lower bounding coordinates to a variable
    upZLetter = upperUTMZone[-1]
    lwZLetter = lowerUTMZone[-1]
    upZ = upperUTMZone[-2]
    lwZ = lowerUTMZone[-2]
    # Compare the strings of the longitudinal zones to see if they match
    if upZ == lwZ:
        print "Zones match, no need to look at another grid\n"
        print "The zone is " + str(upZ)
        # Check the latitudinal zone letter to see if they match
        if upZLetter == lwZLetter:
            print "\nZone letter is " + upZLetter
            zones = [str(upZ), upZLetter, str(lwZ), lwZLetter]
            return zones
        else:
            print "\nZone letters are different\n"
            print "Script will use grid " + str(upZ) + upZLetter + " and " + str(lwZ) + lwZLetter + ".\n"
            zones = [str(upZ), upZLetter, str(lwZ), lwZLetter]
            return zones
    else:
        print "Zones do not match\n"
        print "Script will use grid " + str(upZ) + upZLetter + " and " + str(lwZ) + lwZLetter + ".\n"
        zones = [str(upZ), upZLetter, str(lwZ), lwZLetter]
        return zones

def getGrids(mgrsDir, utmZones):
    # Get list of files in MGRS Grids directory
    mgrsList = os.listdir(mgrsDir)
    #print mgrsList

    # Get the zones and concatenate for comparison to see what files we need
    # to use to compare or event files
    get1stZone = utmZones[0:2]
    firstZone = "".join(get1stZone)
    get2ndZone = utmZones[2:]
    secondZone = "".join(get2ndZone)

    # Test whether the lat/lon zone numbers and letters are the same
    # will return the shapefile name for the zone in question
    if firstZone == secondZone:
        pattern = 'MGRS_1km_' + firstZone
        for file in mgrsList:
            if fnmatch.fnmatch(file, pattern+'*.shp'):
                mgrs1 = file
                print "This is the MGRS Grid File the script will use " + mgrs1
                return [mgrs1]
    # If they are different assign each zone to a varible and return the too files
    else:
        pattern1 = 'MGRS_1km_' + firstZone + '*.shp'
        pattern2 = 'MGRS_1km_' + secondZone + '*.shp'
        for file in mgrsList:
            if fnmatch.fnmatch(file, pattern1):
                mgrs1 = file
            if fnmatch.fnmatch(file, pattern2):
                mgrs2 = file
        print "These are the MGRS Grid Files the script will use " + mgrs1 + ' and ' + mgrs2

        return [mgrs1, mgrs2]

def createBuffer(inShp, outputBuffer, bufferDist, eoPoly):

    # Open EO Shapefile
    source = esriDriver.Open(inShp, 1)
    layerEO = source.GetLayer(0)

    doConvexHull(inShp, eoPoly)

    eoPolySource = esriDriver.Open(eoPoly)
    layerEOPoly = eoPolySource.GetLayer(0)

    # Create output EO buffer
    if os.path.exists(outputBuffer):
        esriDriver.DeleteDataSource(outputBuffer)
    outputBufferds = esriDriver.CreateDataSource(outputBuffer)
    bufferLayer = outputBufferds.CreateLayer(outputBuffer, geom_type=ogr.wkbPolygon)
    featureDefn = bufferLayer.GetLayerDefn()

    eoFeature = layerEOPoly.GetFeature(0)
    eoPolyGeom = eoFeature.GetGeometryRef()
    bufferEOPoly = eoPolyGeom.Buffer(bufferDist)

    outFeature = ogr.Feature(featureDefn)
    outFeature.SetGeometry(bufferEOPoly)
    bufferLayer.CreateFeature(outFeature)
    outFeature = None

def getBoundary(inShp):
    #inShp = shp.split('\\')[-1]
    source = esriDriver.Open(inShp, 1)
    layer = source.GetLayer(0)
    extent = layer.GetExtent()
    print extent
    #feature = layer.GetFeature(0)
    #geom = feature.GetGeometryRef()
    #extent = geom.GetExtent()
    return extent

def doConvexHull(infile, outfile):
    inH = ogr.Open(infile, 0)
    if inH is None:
        print "Could not open file {0}. Exit.".format(infile)
    layer = inH.GetLayer()

    # get all polygons
    geom_collection = ogr.Geometry(ogr.wkbGeometryCollection)
    for feature in layer:
        geometry = feature.GetGeometryRef()
        geom_collection.AddGeometry(geometry)

    convexHull = geom_collection.ConvexHull()

    if os.path.exists(outfile):
        esriDriver.DeleteDataSource(outfile)
    ds = esriDriver.CreateDataSource( outfile )
    if ds is None:
        print "Could not create file {0}".format( outfile)

    # fields
    fldDfn = ogr.FieldDefn('id', ogr.OFTInteger)
    fldDfn.SetWidth(4)

    spatRef = osr.SpatialReference()
    spatRef.ImportFromEPSG(4326)

    lyrname = "convexHull_${0}".format( layer.GetName() )
    lyr = ds.CreateLayer( lyrname, spatRef, ogr.wkbPolygon )
    lyr.CreateField(fldDfn)

    thisFeature = ogr.Feature( lyr.GetLayerDefn() )
    thisFeature.SetGeometry( convexHull )
    thisFeature.SetField('id',1)
    lyr.CreateFeature( thisFeature )

    ds.Destroy()

def clip2shp(inBuffFile, zone, inMGRS, outGridFile):
    inBuff = ogr.Open(inBuffFile, 0)
    if inBuff is None:
        print "Could not open file {0}. Exit.".format(inBuffFile)
    layerBuff = inBuff.GetLayer()

    inGrid = ogr.Open(inMGRS, 0)
    if inGrid is None:
        print "Could not open file {0}. Exit.".format(inMGRS)
    layerGrid = inGrid.GetLayer()

    if os.path.exists(outGridFile):
        esriDriver.DeleteDataSource(outGridFile)
    gridDataSource = esriDriver.CreateDataSource(outGridFile)
    if gridDataSource is None:
        print "Could not create file {0}".format(outGridFile)

    # Get the input spatial reference from the Grid file
    inSpatRef = osr.SpatialReference()
    inSpatRef.ImportFromEPSG(4269)

    # Set the output file from WGS84 Geo to WGS84 UTM
    outSpatRef = osr.SpatialReference()
    outSpatRef.ImportFromEPSG(getEPSG(zone)),1

    # Transform the coordinates
    coordTrans = osr.CoordinateTransformation(inSpatRef,outSpatRef)

    # Create a layer and new shapefile
    layerName = inMGRS.split('\\')[-2]
    lyr = gridDataSource.CreateLayer(layerName, outSpatRef, ogr.wkbPolygon)
    # Copy fields from MGRS Grid layer
    layerGridDef = layerGrid.GetLayerDefn()
    for i in range(layerGridDef.GetFieldCount()):
        lyr.CreateField(layerGridDef.GetFieldDefn(i))

    # Test the buffer and grid intersection and add feature/fields to output
    for i in range(layerBuff.GetFeatureCount()):
        feature1 = layerBuff.GetFeature(i)
        geometry1 = feature1.GetGeometryRef()
        #envelope1 = geometry1.GetEnvelope()
        for i in range(layerGrid.GetFeatureCount()):
            feature2 = layerGrid.GetFeature(i)
            geometry2 = feature2.GetGeometryRef()
            # Create fields for MGRS grid
            attEasting = feature2.GetField("EASTING")
            attNorthing = feature2.GetField("NORTHING")
            attkmSQ_ID = feature2.GetField("kmSQ_ID")
            attGZD = feature2.GetField("GZD")
            attShape_Leng = feature2.GetField("Shape_Leng")
            attMGRS = feature2.GetField("MGRS")
            attMGRS_10km = feature2.GetField("MGRS_10km")
            attShape_Le_1 = feature2.GetField("Shape_Le_1")
            attShape_Area = feature2.GetField("Shape_Area")
            if geometry2.Intersects(geometry1):
                #intersection = geometry1.Intersection(geometry2)
                dstfeature = ogr.Feature(lyr.GetLayerDefn())
                geometry2.Transform(coordTrans)
                dstfeature.SetGeometry(geometry2)
                #print "Get geometry information {0}".format(dstfeature.GetGeometryRef())
                # Create fields for MGRS grid
                dstfeature.SetField("EASTING", attEasting)
                dstfeature.SetField("NORTHING", attNorthing)
                dstfeature.SetField("kmSQ_ID", attkmSQ_ID)
                dstfeature.SetField("GZD", attGZD)
                dstfeature.SetField("Shape_Leng", attShape_Leng)
                dstfeature.SetField("MGRS", attMGRS)
                dstfeature.SetField("MGRS_10km", attMGRS_10km)
                dstfeature.SetField("Shape_Le_1", attShape_Le_1)
                dstfeature.SetField("Shape_Area", attShape_Area)
                # if dstfeature.IsFieldSetAndNotNull("MGRS") is True:
                #     print "Field is set and not null"
                #     print "Get field information {0}".format(dstfeature.GetFieldAsString("MGRS"))
                # else:
                #     print "The field is empty"
                lyr.CreateFeature(dstfeature)
                dstfeature.Destroy()

def shp2tdf(newGridShp):

    fn = newGridShp
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(fn, 1)
    layer = dataSource.GetLayer(0)

    # # Get Field Names
    # msg = "Enter the Project Prefix and Data Type"
    # title = "Modify Field Names"
    # fieldNames = ["Project Prefix", "Data Type"]
    # fieldValues = []  # we start with blanks for the values
    # fieldValues = easygui.multenterbox(msg, title, fieldNames)

    # # make sure that none of the fields was left blank
    # while 1:
    #     if fieldValues == None: break
    #     errmsg = ""
    #     for i in range(len(fieldNames)):
    #         if fieldValues[i].strip() == "":
    #             errmsg = errmsg + ('"%s" is a required field.\n\n' % fieldNames[i])
    #     if errmsg == "": break  # no problems found
    #     fieldValues = easygui.multenterbox(errmsg, title, fieldNames, fieldValues)
    fieldValues = ['2018_FEMA_PR', 'RGBN']
    print "Reply was:", fieldValues

    # Get shapefile name
    fnSplit = fn.split('/')[-1:]
    print "The shapefile name is ", fnSplit[0]

    # Create the Tile Definition File (TDF)
    tdFile = open((fnSplit[0] + '.tdf'), 'w')
    tdfHdr = """**************************************************************

SimActive Tile Definition File
Copyright(c) 2009 SimActive Inc.
All Rights reserved.

**************************************************************
#Tile Name	First Corner(x,y)	Second Corner(x,y)"""

    tdFile.writelines(tdfHdr)

    # Loop through layer to get features and feature geometry
    print "Writing TDF lines to file"
    for i in range(layer.GetFeatureCount()):
        feature = layer.GetFeature(i)
        name = fieldValues[0] + '_' + feature.GetField("MGRS") + '_' + fieldValues[1]
        geometry = feature.GetGeometryRef()
        envelope = geometry.GetEnvelope()
        #print envelope
        # utmMinXMaxY = utm.from_latlon(envelope[3], envelope[0])
        # print utmMinXMaxY
        # utmMaxXMinY = utm.from_latlon(envelope[2], envelope[1])
        # print utmMaxXMinY
        # print "minX: %.3f, maxY: %.3f, maxX: %.3f, minY: %.3f" %(round(utmMinXMaxY[0]),round(utmMinXMaxY[1]),round(utmMaxXMinY[0]),round(utmMaxXMinY[1]))
        print "minX: %.3f, maxY: %.3f, maxX: %.3f, minY: %.3f" % ((envelope[0]), (envelope[3]), (envelope[1]), (envelope[2]))
        tdfLines = ('\n%s %.3f %.3f %.3f %.3f') % (name, round(envelope[0]), round(envelope[3]), round(envelope[1]), round(envelope[2]))
        # tdfLines = ('\n%s %.3f %.3f %.3f %.3f\n') % (name, round(utmMinXMaxY[0]), round(utmMinXMaxY[1]), round(utmMaxXMinY[0]), round(utmMaxXMinY[1]))
        print tdfLines+ '\n'
        tdFile.writelines(tdfLines)

    tdFile.close()

def getEPSG(zone):
    if zone == '20':
        epsg = 26919
        return epsg
    elif zone == '19':
        epsg = 26919
        return epsg
    elif zone == '18':
        epsg = 26918
        return epsg
    elif zone == '17':
        epsg = 26917
        return epsg
    elif zone == '16':
        epsg = 26916
        return epsg
    elif zone == '15':
        epsg = 26915
        return epsg
    elif zone == '14':
        epsg = 26914
        return epsg
    elif zone == '13':
        epsg = 26913
        return epsg
    elif zone == '12':
        epsg = 26912
        return epsg
    elif zone == '11':
        epsg = 26911
        return epsg
    elif zone == '10':
        epsg = 26910
        return epsg
    elif zone == '9':
        epsg = 26909
        return epsg
    elif zone == '8':
        epsg = 26908
        return epsg
    elif zone == '7':
        epsg = 26907
        return epsg
    elif zone == '6':
        epsg = 26906
        return epsg
    elif zone == '5':
        epsg = 26905
        return epsg
    elif zone == '4':
        epsg = 26904
        return epsg

###>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.MAIN.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>###

if __name__ == '__main__':
    eoFileList = easygui.fileopenbox(msg="Select IPAS CO+ EO Files", title= "EO Selection", default='*', filetypes=".txt", multiple=True)
    eos2shp(eoFileList)
