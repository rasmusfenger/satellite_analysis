import gdal
import numpy as np
import os
import glob


site = 'kangeq'
sensor = 'sequoia'
channel = False
overwrite = False

inFolder = os.path.join('/Volumes/RASMUS_1/Satellite/remains_sites', site, sensor)
# pleiades
#inFile = glob.glob(inFolder + '/*_ms_toa_subset.tif')[0]
#inFile = glob.glob(inFolder + '/*_ms_atmcorr_subset.tif')[0]
# sequoia
inFile = glob.glob(inFolder + '/*_seqirrad_empcorr.tif')[0]
# Sentinel-2
#inFile = '/Volumes/RASMUS_1/Satellite/remains_sites/sandnes/S2/toa/S2A_OPER_MSI_L1C_TL_SGS__20160728T202436_A005742_T22WES_subset_refl.tif'
#inFile = '/Volumes/RASMUS_1/Satellite/remains_sites/iffiartafik/S2/toa/S2A_OPER_MSI_L1C_TL_SGS__20160728T202436_A005742_T22WES_subset_refl.tif'
#inFile = '/Volumes/RASMUS_1/Satellite/remains_sites/qoornoq/S2/toa/S2A_OPER_MSI_L1C_TL_SGS__20160728T202436_A005742_T22WDS_subset_refl.tif'
#inFile = '/Volumes/RASMUS_1/Satellite/remains_sites/ersaa/S2/toa/S2A_OPER_MSI_L1C_TL_SGS__20160728T202436_A005742_T22WDS_subset_refl.tif'
#inFile = '/Volumes/RASMUS_1/Satellite/remains_sites/kangeq/S2/toa/S2A_OPER_MSI_L1C_TL_SGS__20160814T202041_A005985_T22WDS_subset_refl.tif'

if channel:
    outFolder = os.path.join(inFolder, 'channel')
else:
    outFolder = os.path.join(inFolder, 'VI')

def oneband(band):
    return band

def evi(blue,red,nir):
    evi = 2.5*((nir-red)/(nir + 6*red - 7.5*blue + 1))
    return evi

def greenndvi(green,nir):
    greenndvi = (nir - green) / (nir + green)
    return greenndvi

def ndvi(red,nir):
    ndvi = (nir - red) / (nir + red)
    return ndvi

def sr(red,nir):
    sr = nir/red
    return sr

def msr(red,nir):
    msr = red/(nir/red+1)**0.5
    return msr

def mtvi2(green,red,nir):
    mtvi2 = 1.5*(1.2*(nir-green)-2.5*(red-green))*(2*nir+1)**2-(6*nir-5*red)-0.5
    return mtvi2

def rdvi(red,nir):
    rdvi = (nir - red) / (nir + red)**0.5
    return rdvi

def irg(green,red):
    irg = red-green
    return irg

def pvi(red,nir):
    pvi = (((0.851*red+0.355*nir)-red)**2 + ((0.355*red + 0.148*nir)-nir)**2)**(0.5)
    return pvi

def rvi(red,nir):
    rvi = red/nir
    return rvi

def arvi(blue,red,nir):
    prb = red - 1*(blue-red)
    arvi = (nir - prb)/(nir + prb)
    return arvi

def msavi(red,nir):
    msavi = (2*nir+1-((2*nir+1)**2 - 8*(nir-red))**0.5)/2
    return msavi

def gemi(red,nir):
    n = (2*(nir**2-red**2)+1.5*nir+0.5*red)/(nir+red+0.5)
    gemi = n*(1-0.25*n)*(red-0.125)/(1-red)
    return gemi

def sarvi(blue,red,nir):
    prb = red - 1*(blue-red)
    sarvi = (1 + 0.5) * (nir - prb)/(nir + prb + 0.5)
    return sarvi

def osavi(red,nir):
    osavi = (nir-red)/(nir+red+0.16)
    return osavi

def dvi(red,nir):
    dvi = nir-red
    return dvi

def srndvi(red, nir):
    srndvi = (nir**2-red)/(nir+red**2)
    return srndvi

# following VI's are not included for now

def srndvi(red,nir):
    srndvi = (nir**2-red)/(nir+red**2)
    return srndvi

def gbndvi(blue,green,nir):
    gbndvi = (nir-(green+blue)/(nir+(green+blue)))
    #GBNDVI=[NIR-(Green+Blue)]/ [NIR+(Green+Blue)]
    return gbndvi

def rbndvi(blue,red,nir):
    rbndvi = (nir - (red + blue) / (nir + (red + blue)))
    #RBNDVI = [NIR - (Red + Blue)] / [NIR + (Red + Blue)]
    return rbndvi

def bluendvi(blue,nir):
    ndviblue = (nir - blue) / (nir + blue)
    return ndviblue

def ndvigreenred(red,green):
    ndvigreenred = (green - red) / (green + red)
    return ndvigreenred

def multiply(blue,green,red,nir):
    ndvi = (nir - red) / (nir + red)
    multiply = np.where(ndvi>0.35,blue*green*red*nir,0)
    return multiply

def array2tif(array, inFile, outFile):
    ds = gdal.Open(inFile)
    drv = gdal.GetDriverByName('GTiff')
    outTif = drv.Create(outFile, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Float32)
    outTif.SetGeoTransform(ds.GetGeoTransform())
    outTif.SetProjection(ds.GetProjection())
    outTif.GetRasterBand(1).WriteArray(array)
    outTif.GetRasterBand(1).SetNoDataValue(-999)
    outTif = None

def viArrayCalc(alg,blue,green,red,nir,rededge):
    if alg == 'blue':
        vi = oneband(blue)
    elif alg == 'green':
        vi = oneband(green)
    elif alg == 'red':
        vi = oneband(red)
    elif alg == 'rededge':
        vi = oneband(rededge)
    elif alg == 'nir':
        vi = oneband(nir)
    elif alg == 'evi':
        vi = evi(blue, red, nir)
    elif alg == 'greenndvi':
        vi = greenndvi(green, nir)
    elif alg == 'ndvi':
        vi = ndvi(red, nir)
    elif alg == 'sr':
        vi = sr(red, nir)
    elif alg == 'msr':
        vi = msr(red,nir)
    elif alg == 'mtvi2':
        vi = mtvi2(green, red, nir)
    elif alg == 'rdvi':
        vi = rdvi(red, nir)
    elif alg == 'irg':
        vi = irg(green, red)
    elif alg == 'pvi':
        vi = pvi(red,nir)
    elif alg == 'rvi':
        vi = rvi(red,nir)
    elif alg == 'arvi':
        vi = arvi(blue,red,nir)
    elif alg == 'dvi':
        vi = dvi(red,nir)
    elif alg == 'srndvi':
        vi = srndvi(red,nir)
    elif alg == 'msavi':
        vi = msavi(red, nir)
    elif alg == 'gemi':
        vi = gemi(red, nir)
    elif alg == 'sarvi':
        vi = sarvi(blue, red, nir)
    elif alg == 'osavi':
        vi = osavi(red, nir)
    return vi

def correct(array):
    array = np.where(array <= 0.00001, 0.00001, array)
    array = np.where(array > 1., 1., array)
    return array

def calcVI(inFile, outFolder, sensor):
    viListPleiades = ['evi',
                      'greenndvi',
                      'ndvi',
                      'sr',
                      'msr',
                      'mtvi2',
                      'rdvi',
                      'irg',
                      'pvi',
                      'rvi',
                      'msavi',
                      'arvi',
                      'gemi',
                      'sarvi',
                      'osavi',
                      'dvi',
                      'srndvi']

    chListPleiades = ['blue',
                      'green',
                      'red',
                      'nir']

    viListSequoia = ['greenndvi',
                     'ndvi',
                     'sr',
                     'msr',
                     'mtvi2',
                     'rdvi',
                     'irg',
                     'pvi',
                     'rvi',
                     'msavi',
                     'gemi',
                     'osavi',
                     'dvi',
                     'srndvi']

    chListSequoia = ['green',
                     'red',
                     'rededge',
                     'nir']

    viListSony = ['ndvigreenred',
                  'irg']

    chListSony = ['blue',
                  'green',
                  'red']

    ds = gdal.Open(inFile)
    if sensor == 'pleiades' or sensor == 'spot':
        blue = ds.GetRasterBand(1).ReadAsArray().astype(np.float32)
        green = ds.GetRasterBand(2).ReadAsArray().astype(np.float32)
        red = ds.GetRasterBand(3).ReadAsArray().astype(np.float32)
        nir = ds.GetRasterBand(4).ReadAsArray().astype(np.float32)
        rededge = None
        if channel:
            viList = chListPleiades
        else:
            viList = viListPleiades
    elif sensor == 'sequoia':
        blue = None
        green = ds.GetRasterBand(1).ReadAsArray().astype(np.float32)#/32768.
        red = ds.GetRasterBand(2).ReadAsArray().astype(np.float32)#/32768.
        rededge = ds.GetRasterBand(3).ReadAsArray().astype(np.float32)#/32768.
        nir = ds.GetRasterBand(4).ReadAsArray().astype(np.float32)#/32768.
        if channel:
            viList = chListSequoia
        else:
            viList = viListSequoia

    elif sensor == 'sony':
        red = ds.GetRasterBand(1).ReadAsArray().astype(np.float32)/255
        green = ds.GetRasterBand(2).ReadAsArray().astype(np.float32)/255
        blue = ds.GetRasterBand(3).ReadAsArray().astype(np.float32)/255
        rededge = None
        nir = None
        if channel:
            viList = chListSony
        else:
            viList = viListSony

    elif sensor == 'S2':
        blue = ds.GetRasterBand(2).ReadAsArray().astype(np.float32)
        green = ds.GetRasterBand(3).ReadAsArray().astype(np.float32)
        red = ds.GetRasterBand(4).ReadAsArray().astype(np.float32)
        nir = ds.GetRasterBand(8).ReadAsArray().astype(np.float32)
        rededge = None
        if channel:
            viList = chListPleiades
        else:
            viList = viListPleiades

    # correct data. Set minimum to 0.0001
    blue = correct(blue)
    green = correct(green)
    red = correct(red)
    nir = correct(nir)
    rededge = correct(rededge)

    root, ext = os.path.splitext(inFile)
    head, tail = os.path.split(root)
    if overwrite:
        for the_file in os.listdir(outFolder):
            file_path = os.path.join(outFolder, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
            except Exception as e:
                print(e)

    for vi in viList:
        outFile = os.path.join(outFolder, tail + '_' + vi + '.tif')
        if os.path.isfile(outFile) == False:
            viArray = viArrayCalc(vi,blue,green,red,nir,rededge)
            array2tif(viArray, inFile, outFile)
            print 'Created: ', outFile
calcVI(inFile, outFolder, sensor)


print 'finished'