import os
from fiona import collection
from fiona.crs import from_epsg
from shapely.geometry import Point, mapping, shape
from rasterstats import zonal_stats

def make_buf(inShp, buf, outFile=None):
    # Creates polygon shapefile as buffers from input point shapefile
    # inShp   = points shapefile
    # buf     = buffer distance in meter
    # outFile = polygon shapefile
    epsg_32622 = 'PROJCS["WGS_1984_UTM_Zone_22N",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-51],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["Meter",1]]'
    dir = os.path.dirname(os.path.abspath(__file__))
    tempFolder = os.path.join(dir, 'temp')
    if not os.path.isdir(tempFolder):
        os.mkdir(tempFolder)

    if not outFile:
        head,tail = os.path.split(inShp)
        outFile = os.path.join(tempFolder, tail)

    if not os.path.isfile(outFile):
        with collection(inShp, 'r') as input:
            schema = input.schema.copy()
            schema['geometry'] = 'Polygon'
            with collection(outFile, "w", crs=from_epsg(32622), driver="ESRI Shapefile", schema=schema) as output:
                for point in input:
                    prop = point['properties']
                    output.write({'properties': prop,
                                    'geometry': mapping(shape(point['geometry']).buffer(buf))})
        prjfile = outFile[:-3] + 'prj'
        with open(prjfile, 'w') as prj:
            prj.write(epsg_32622)
        print 'Shapefile created: ' + outFile
    return outFile

def extract_from_img(polyShp, imgFile, band=1, all_touched=True):
    # extracts zonal statistics from band 1 in raster using polygon shapefile as zones
    stats = zonal_stats(polyShp, imgFile, band=band, all_touched=all_touched, stats=['mean', 'std', 'min', 'max', 'count'])
    return stats

'''
import glob
#viFolder = '/Volumes/RASMUS_1/Satellite/remains_sites/ersaa/pleiades/VI'
#viFolder = '/Volumes/RASMUS_1/Satellite/remains_sites/ersaa/spot/VI'
#viFolder = '/Volumes/RASMUS_1/Satellite/remains_sites/ersaa/sequoia/VI'

viFolder = '/Volumes/RASMUS_1/Satellite/remains_sites/iffiartafik/pleiades/VI'
#viFolder = '/Volumes/RASMUS_1/Satellite/remains_sites/iffiartafik/spot/VI'
#viFolder = '/Volumes/RASMUS_1/Satellite/remains_sites/iffiartafik/sequoia/VI'

sensor = 'pleiades'
imgList = glob.glob(viFolder + '/*.tif')

results_natural = {}
results_cultural = {}
results_contrast = {}
#print 'image\ttype\tmean\tstd\tcontrast'
for img in imgList:
    # calculate natural background
    inShpNat = '/Volumes/RASMUS_1/Satellite/analysis/shapefiles/iffiartafik_plots_natural.shp'
    polyShpBack = make_buf(inShpNat, 10)
    natural = extract_from_img(polyShpNat, img, band=1)

    # calculate cultural
    inShpCul = '/Volumes/RASMUS_1/Satellite/analysis/shapefiles/iffiartafik_plots_cultural.shp'
    polyShp = make_buf(inShp, 10)
    cultural = extract_from_img(polyShp, img, band=1)

    imgSplit = img.split('_')
    vi = imgSplit[len(imgSplit)-1][:-4]
    results_natural[vi] = natural[0]
    results_cultural[vi] = cultural[0]

import numpy as np
import matplotlib.pyplot as plt

N = len(results_natural)
label = []
mean_natural = []
std_natural = []
max_natural = []
min_natural = []
mean_cultural = []
std_cultural = []
max_cultural = []
min_cultural = []

for key in results_natural:
    label.append(key)
    mean_natural.append(results_natural[key]['mean'])
    std_natural.append(results_natural[key]['std'])
    max_natural.append(results_natural[key]['max'])
    min_natural.append(results_natural[key]['min'])
for key in results_cultural:
    mean_cultural.append(results_cultural[key]['mean'])
    std_cultural.append(results_cultural[key]['std'])
    max_cultural.append(results_cultural[key]['max'])
    min_cultural.append(results_cultural[key]['min'])

ind = np.arange(N)  # the x locations for the groups
width = 0.35       # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind, mean_natural, width, color='r', yerr=std_natural)

rects2 = ax.bar(ind + width, mean_cultural, width, color='y', yerr=std_cultural)

# add some text for labels, title and axes ticks
ax.set_ylabel('VI')
ax.set_title('Vegetation indices ' + sensor)
ax.set_xticks(ind + width / 2)
ax.set_xticklabels(label)

ax.legend((rects1[0], rects2[0]), ('natural', 'cultural'))


def autolabel(rects):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%d' % int(height),
                ha='center', va='bottom')
#autolabel(rects1)
#autolabel(rects2)

fig2, ax2 = plt.subplots()

min = np.where(np.array(min_cultural) < np.array(min_natural), np.array(min_cultural), np.array(min_natural))
max = np.where(np.array(max_cultural) > np.array(max_natural), np.array(max_cultural), np.array(max_natural))

mean_cul_shifted = np.array(mean_cultural) - min
mean_nat_shifted = np.array(mean_natural) - min
max_shifted = max - min

contrast = (mean_cul_shifted - mean_nat_shifted) / max_shifted
ax2.bar(ind, contrast, width)
ax2.set_xticks(ind)
ax2.set_xticklabels(label)

print mean_natural
print mean_cultural
#print mean_cul_shifted
#print mean_nat_shifted
#print max_shifted
plt.show()
'''