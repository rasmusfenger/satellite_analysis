import os
import numpy as np
import gdal
import subprocess
import glob
from shapefunctions import *
import sys
sys.path.append('/Users/rasmus/PycharmProjects/fieldsamples')
from stattest import ttest_manual

def refl2dn(inFile, outFile, quantificationValue=65535):
    ds = gdal.Open(inFile)
    drv = gdal.GetDriverByName('GTiff')
    outTif = drv.Create(outFile, ds.RasterXSize, ds.RasterYSize, ds.RasterCount, gdal.GDT_UInt16)
    outTif.SetGeoTransform(ds.GetGeoTransform())
    outTif.SetProjection(ds.GetProjection())
    for band in range(ds.RasterCount):
        band += 1
        array = ds.GetRasterBand(band).ReadAsArray().astype(np.float32)
        array = np.where(array < 0, 0, np.round(array*quantificationValue,decimals=0))
        outTif.GetRasterBand(band).WriteArray(array)
    outTif = None
    #print 'processed: ', outFile


def dn2refl(inFile, outFile, quantificationValue):
    ds = gdal.Open(inFile)
    drv = gdal.GetDriverByName('GTiff')
    outTif = drv.Create(outFile, ds.RasterXSize, ds.RasterYSize, ds.RasterCount, gdal.GDT_Float32)
    outTif.SetGeoTransform(ds.GetGeoTransform())
    outTif.SetProjection(ds.GetProjection())
    for band in range(ds.RasterCount):
        band += 1
        array = ds.GetRasterBand(band).ReadAsArray().astype(np.float32)
        array = array / quantificationValue
        outTif.GetRasterBand(band).WriteArray(array)
    outTif = None
    #print 'processed: ', outFile

def clip_raster(inFile, outFile, clipFile):
    # Get coordinates of clipFile
    src = gdal.Open(clipFile)
    ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
    lrx = ulx + (src.RasterXSize * xres)
    lry = uly + (src.RasterYSize * yres)

    # Build command
    command = ['gdal_translate',
               '-of', 'GTiff',
               '-projwin', str(ulx), str(uly), str(lrx), str(lry),
               inFile,
               outFile]
    # -ot', 'Uint16'
    cmd = ''
    for c in command:
        cmd += c + ' '
    # run command
    subprocess.check_call(cmd,shell=True)
    #print 'processed: ', outFile

def clip_match_raster(inFile, outFile, clipFile):
    # Get coordinates of clipFile
    src = gdal.Open(clipFile)
    ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
    lrx = ulx + (src.RasterXSize * xres)
    lry = uly + (src.RasterYSize * yres)

    # Build command
    command = ['gdal_translate',
               '-of', 'GTiff',
               '-b', '1',
               '-b', '2',
               '-b', '3',
               '-b', '4',
               '-projwin', str(ulx), str(uly), str(lrx), str(lry),
               '-tr', str(xres), str(yres),
               '-r', 'average',
               inFile,
               outFile]
    # -ot', 'Uint16'
    cmd = ''
    for c in command:
        cmd += c + ' '
    # run command
    subprocess.check_call(cmd,shell=True)
    #print 'processed: ', outFile

def getTimestamp(x):
    head,tail = os.path.split(x)
    timestamp = tail.split('_')[7]
    return(timestamp)

def timebands(inFileList, outFolder):
    inFileList = sorted(inFileList, key=getTimestamp)
    ds1 = gdal.Open(inFileList[0])
    bandCount = ds1.RasterCount
    drv = gdal.GetDriverByName('GTiff')
    # write timestamp to txt file
    with open(os.path.join(outFolder, 'band_timestamp.txt'), 'w') as txt:
        firstRun = True
        for band in range(bandCount):
            band += 1
            outFile = os.path.join(outFolder, 'band_' + str(band) + '.tif')
            outTif = drv.Create(outFile, ds1.RasterXSize, ds1.RasterYSize, len(inFileList), gdal.GDT_Float32)
            outTif.SetGeoTransform(ds1.GetGeoTransform())
            outTif.SetProjection(ds1.GetProjection())
            writeTime = True
            for num in range(len(inFileList)):
                if firstRun and writeTime:
                    #head,tail = os.path.split(inFileList[num])
                    #timestamp = tail.split('_')[7]
                    timestamp = getTimestamp(inFileList[num])
                    txt.write(timestamp + '\n')
                ds = gdal.Open(inFileList[num])
                array = ds.GetRasterBand(band).ReadAsArray().astype(np.float32)
                outTif.GetRasterBand(num+1).WriteArray(array)
            outTif = None
            writeTime = False
            firstRun = False

def extractband(inFile, band):
    outFile = inFile[:-4] + '_' + str(band) + '.tif'
    ds = gdal.Open(inFile)
    array = ds.GetRasterBand(band).ReadAsArray()
    drv = gdal.GetDriverByName('GTiff')
    outTif = drv.Create(outFile, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_UInt16)
    outTif.SetGeoTransform(ds.GetGeoTransform())
    outTif.SetProjection(ds.GetProjection())
    outTif.GetRasterBand(1).WriteArray(array)
    outTif = None

def resample(inFile, outRes):
    outFile = inFile[:-4] + '_' + str(outRes) + 'm.tif'
    # Build command
    command = ['gdal_translate',
               '-of', 'GTiff',
               '-b', '1',
               '-b', '2',
               '-b', '3',
               '-b', '4',
               '-tr', str(outRes), str(outRes),
               '-r', 'average',
               inFile,
               outFile]
    cmd = ''
    for c in command:
        cmd += c + ' '
    # run command
    subprocess.check_call(cmd,shell=True)
    #print 'processed: ', outFile

def raster_extract(imgFile, polyShp, all_touched='all_touched'):
    # prepare dictionary
    meanList = []
    stdList = []
    countList = []
    ds = gdal.Open(imgFile)
    for band in range(ds.RasterCount):
        band += 1
        stats = extract_from_img(polyShp, imgFile, band, all_touched)
        meanList.append(stats[0]['mean'])
        stdList.append(stats[0]['std'])
        countList.append(stats[0]['count'])
    return meanList, stdList, countList

def extract_plot_values(site, sensor, folder, all_touched=True, print_p=False):
    # Extract values within plots from images
    # test cultural vs natural
    viFolder = os.path.join('/Volumes/RASMUS_1/Satellite/remains_sites',site,sensor,folder)
    imgList = glob.glob(viFolder + '/*.tif')

    results_natural = {}
    results_cultural = {}

    for img in imgList:
        # calculate natural background
        inShpNat = '/Volumes/RASMUS_1/Satellite/analysis/shapefiles/'+site+'_plots_natural_'+sensor+'.shp'
        natural = extract_from_img(inShpNat, img, band=1, all_touched=all_touched)
        natMean = []
        for elem in natural:
            natMean.append(elem['mean'])
        # calculate cultural
        inShpCul = '/Volumes/RASMUS_1/Satellite/analysis/shapefiles/'+site+'_plots_cultural_'+sensor+'.shp'
        cultural = extract_from_img(inShpCul, img, band=1, all_touched=all_touched)
        culMean = []
        for elem in cultural:
            culMean.append(elem['mean'])

        imgSplit = img.split('_')
        vi = imgSplit[len(imgSplit)-1][:-4]
        results_natural[vi] = natMean
        results_cultural[vi] = culMean
    pList = []
    for key in results_natural:
        if print_p:
            print key
        t, p = ttest_manual(results_natural[key], results_cultural[key], equal_var=False, print_it=print_p)
        pList.append(p)
    return results_natural,results_cultural,pList

def extract_poly_values(site, sensor, folder, all_touched=True, print_p=False):
    from scipy.stats import ttest_ind_from_stats
    # Extract values within plots from images
    # test cultural vs natural
    viFolder = os.path.join('/Volumes/RASMUS_1/Satellite/remains_sites',site,sensor,folder)
    imgList = glob.glob(viFolder + '/*.tif')

    results_natural = {}
    results_cultural = {}

    for img in imgList:
        # calculate natural background
        inShpNat = '/Volumes/RASMUS_1/Satellite/analysis/shapefiles/'+site+'_natural_'+sensor+'_poly.shp'
        natural = extract_from_img(inShpNat, img, band=1, all_touched=all_touched)
        # calculate cultural
        inShpCul = '/Volumes/RASMUS_1/Satellite/analysis/shapefiles/'+site+'_cultural_'+sensor+'_poly.shp'
        cultural = extract_from_img(inShpCul, img, band=1, all_touched=all_touched)

        imgSplit = img.split('_')
        vi = imgSplit[len(imgSplit)-1][:-4]
        results_natural[vi] = natural[0]
        results_cultural[vi] = cultural[0]
    pList = []

    for key in results_natural:
        if print_p:
            print key
        t, p = ttest_ind_from_stats(results_natural[key]['mean'], results_natural[key]['std'],
                                          results_natural[key]['count'], results_cultural[key]['mean'],
                                          results_cultural[key]['std'], results_cultural[key]['count'], equal_var=False)
        pList.append(p)
    return results_natural,results_cultural,pList