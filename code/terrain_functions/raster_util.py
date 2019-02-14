# ===============================================================================
# Copyright 2019 Jan Hendrickx and Gabriel Parrish
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ===============================================================================
import gdal
import os
from numpy import array, reshape
# ============= standard library imports ========================

def convert_raster_to_array(input_raster_path, raster=None, band=1):
    """
    Convert .tif raster into a numpy numerical array.

    :rtype: object
    :param input_raster_path: Path to raster.
    :param raster: Raster name with \*.tif
    :param band: Band of raster sought.

    :return: Numpy array.
    """
    # print "input raster path", input_raster_path
    # print "raster", raster
    p = input_raster_path
    if raster is not None:
        p = os.path.join(p, raster)

    # print "filepath", os.path.isfile(p)
    # print p
    if not os.path.isfile(p):
        print('Not a valid file: {}'.format(p))

    raster_open = gdal.Open(p)
    ras = array(raster_open.GetRasterBand(band).ReadAsArray(), dtype=float)
    return ras


def convert_array_to_raster(output_path, arr, geo, output_band=1):
    driver = gdal.GetDriverByName('GTiff')
    out_data_set = driver.Create(output_path, geo['cols'], geo['rows'],
                                 geo['bands'], geo['data_type'])
    out_data_set.SetGeoTransform(geo['geotransform'])
    out_data_set.SetProjection(geo['projection'])

    output_band = out_data_set.GetRasterBand(output_band)
    print('size of array', arr.shape)
    print('cols {} and rows {}'.format(geo['cols'], geo['rows']))
    print('geotransform {}'.format(geo['geotransform']))
    output_band.WriteArray(arr, 0, 0)
    del out_data_set, output_band

    if not os.path.isfile(output_path):
        print("Not a valid file: '{}' - Raster could not be written!".format(output_path))
        return

# todo - modify this to take an individual path
def get_raster_geo_attributes(root):
    """
    Creates a dict of geographic attributes from any of the pre-processed standardized rasters.

    :param root: Path to a folder with pre-processed standardized rasters.
    :return: dict of geographic attributes.
    """
    # statics = [filename for filename in os.listdir(statics_path) if filename.endswith('.tif')]
    # file_name = statics[0]
    file_name = next((fn for fn in sorted(os.listdir(root)) if fn.endswith('.tif')), None)
    dataset = gdal.Open(os.path.join(root, file_name))

    band = dataset.GetRasterBand(1)
    raster_geo_dict = {'cols': dataset.RasterXSize, 'rows': dataset.RasterYSize, 'bands': dataset.RasterCount,
                       'data_type': band.DataType, 'projection': dataset.GetProjection(),
                       'geotransform': dataset.GetGeoTransform(), 'resolution': dataset.GetGeoTransform()[1]}
    return raster_geo_dict


def remake_arr(arr, shp):
    """
    :param arr: numpy array
    :param shp: int or tuple of ints
    :return: 2D array
    """
    arr = reshape(arr, shp)
    return arr


