import os
import numpy as np

# TODO - MAKE SURE THAT YOU USE PYMETRIC-WRITTEN FUNCTIONS TO READ THE RASTERS AS ARRAYS DONT REPEAT FUNCTIONS ALREADY WRITTEN BY DRI


def calc_solar_topo_image_variables(DOY, LatDeg, SlopeDeg, AspectDeg, local_time, Lz, Lm ):
    """"""
    # This function calculates the cosine and sine for variables that remain constant for a given Landsat image.
    #
    # The solar variables are: LatDeg (latitude), DOY (day of year)
    # The topographical variables are SlopeDeg (slope of a pixel) and AspectDeg (aspect of a pixel)
    #


    # CONVERT FROM DEGREES TO RADIANS
    LatRad = LatDeg * np.pi / 180
    DeclRad = 0.409 * np.sin((2 * np.pi / 365 * DOY) - 1.39)
    DeclDeg = DeclRad * 180 / np.pi
    SlopeRad = SlopeDeg * np.pi / 180
    AspectRad = AspectDeg * np.pi / 180
    # CALCULATE SIN AND COSIN FOR IMAGE VARIABLES (1 OR 2D)
    cos_lat = np.cos(LatRad)
    sin_lat = np.sin(LatRad)
    cos_decl = np.cos(DeclRad)
    sin_decl = np.sin(DeclRad)
    cos_slope = np.cos(SlopeRad)
    sin_slope = np.sin(SlopeRad)
    cos_aspect = np.cos(AspectRad)
    sin_aspect = np.sin(AspectRad)

    # SOLAR CONSTANT ???
    bSc = 2 * np.pi * (DOY - 81) / 364
    Sc = 0.1645 * np.sin(2 * bSc) - 0.1255 * np.cos(bSc) - 0.025 * np.sin(bSc)
    solar_time = (local_time + 0.6667 * (Lz - Lm) + Sc)

    # Hour Angle (w) Lower case omega in Allen 2006
    HourAngleRad = np.pi / 12 * ((local_time + 0.6667 * (Lz - Lm) + Sc) - 12)
    cos_hourangle = np.cos(HourAngleRad)
    sin_hourangle = np.sin(HourAngleRad)

    # The variables a, b, and c are needed for the integration of cos_theta over a day or part of a day
    a = sin_decl * cos_lat * sin_slope * cos_aspect - sin_decl * sin_lat * cos_slope
    b = cos_decl * cos_lat * cos_slope + cos_decl * sin_lat * sin_slope * cos_aspect
    c = cos_decl * sin_aspect * sin_slope

    # From Allen 2010 todo - elaborate on comment
    cos_theta_unadj = - a + b * cos_hourangle + c * sin_hourangle
    cos_theta_adj = cos_theta_unadj / cos_slope
    if cos_theta_adj < 0.1:
        cos_theta_adj = 0.1
    if cos_theta_adj > 10.0:
        cos_theta_adj = 10.0

    # cos_theta_horizontal_pixel at satellite overpass with Eq.[4] in Allen (2006)
    # solar incidence angle at satellite overpass on horizontal pixel is solar_incidence_angle_theta
    # solar elevation angle at satellite  overpass on horizontal pixel is solar_elevation_angle
    # solar azimuth angle at satellite overpass is solar_azimuth_angle with Eq. 11.5 in Campbell and Norman (1998)
    #    where solar_azimuth_angle is calculated with respect to due south, increasing in the counter clockwise direction
    #    so 90 degrees is east. Afternoon azimuth angles can be calculated by taking 360 degrees minus the solar_azimuth angle
    # solar_azimuth_angle_Landsat (also used by NOAA) is measured clockwise starting at 0 in the north after Allen, 2006

    cos_theta_horizontal_pixel = sin_decl * sin_lat + cos_decl * cos_lat * cos_hourangle
    solar_incidence_angle_theta = np.arccos(cos_theta_horizontal_pixel) * 180 / np.pi
    solar_elevation_angle = 90.0 - solar_incidence_angle_theta
    cos_solar_azimuth_angle = -(sin_decl - cos_theta_horizontal_pixel * sin_lat) / \
                              (cos_lat * np.sin(solar_incidence_angle_theta * np.pi / 180))
    solar_azimuth_angle = np.arccos(cos_solar_azimuth_angle) * 180 / np.pi

    # correcting the landsat solar azimuth angle based on hour angle rad less than zero (morning) vs greater than zero (evening)
    if HourAngleRad <= 0.000:
        solar_azimuth_angle_landsat = 180 - solar_azimuth_angle
    solar_azimuth_angle_landsat = 180 + solar_azimuth_angle

    return LatRad, DeclRad, DeclDeg, SlopeRad, AspectRad, cos_lat, sin_lat, cos_slope, sin_slope, \
    cos_aspect, sin_aspect, sin_decl, cos_decl, bSc, Sc, HourAngleRad, cos_hourangle, sin_hourangle, a, b, c, cos_theta_unadj, \
    cos_theta_adj, cos_theta_horizontal_pixel, solar_incidence_angle_theta, \
           solar_elevation_angle, solar_azimuth_angle, solar_azimuth_angle_landsat, solar_time


def calculate_integration_limits(DOY, LatDeg, SlopeDeg, AspectDeg, local_time, Lz, Lm, raster_mode=False):
    """"""
    # Derived variables
    #
    # DeclRad = Declination of the earth (radians; positive in summer in northern hemisphere)
    # HourAngleRad = hour angle (radians)
    # HourAngleRad: 0 at solar noon, negative in the morning, positive in afternoon
    # Sc = seasonal correction for solar time (hours)
    #
    # cos_theta_unadj = unadjusted cosine of the solar incidence angle
    # cos_theta_adj = cos_theta_unadj/cos_slope to convert cos_theta_unadj to a horizontal equivalent
    #     so that final calculation for solar energy is expressed as energy per unit of horizontal area
    #     The minimum value for cos_theta_adj is limited to 0.1 to account for the energy
    #     contribution by diffuse solar radiation even if beam radiation is blocked by the terrain.
    #     IMPORTANT: cos_theta_adj is NOT a cosine, it's the ratio of two cosine values.
    #                Therefore, its value can vary from zero for flat pixels to 3.0 for slope
    #                70 degree, to 10.0 for slope 84.266.
    # threshold_cos_theta_adj = 10.0 is used to prevent extreme values.
    # Projects with a focus on steep slopes may want to change this threshold value.

    if raster_mode:
        # TODO - read in each raster parameter as a numpy array and continue
        pass


    LatRad, DeclRad, DeclDeg, SlopeRad, AspectRad, cos_lat, sin_lat, cos_slope, sin_slope, cos_aspect, sin_aspect, \
    sin_decl, cos_decl, bSc, Sc, HourAngleRad, cos_hourangle, sin_hourangle, a, b, c, cos_theta_unadj, cos_theta_adj, \
    cos_theta_horizontal_pixel, solar_incidence_angle_theta, solar_elevation_angle, solar_azimuth_angle, \
    solar_azimuth_angle_landsat, solar_time = \
        calc_solar_topo_image_variables(DOY, LatDeg, SlopeDeg, AspectDeg, local_time, Lz, Lm)

    if not raster_mode:
        print('\n',
              'Output Function calc_solar_topo_image_variables(DOY, LatDeg, SlopeDeg, AspectDeg, local_time, Lz, Lm)')
        print('\n', 'DOY = {:.0f}'.format(DOY), '    ', 'LatDeg = {:.2f}'.format(LatDeg), '   ',
              'LatRad = {:.4f}'.format(LatRad),
              '   ', 'local_time = {:.2f}'.format(local_time), '  ', 'solar_time = {:.2f}'.format(solar_time), '  '
              , 'Lz = {:.0f}'.format(Lz), '  ', 'Lm = {:.0f}'.format(Lm))
        print('\n', 'solar_incidence_angle_theta = {:.1f}'.format(solar_incidence_angle_theta), \
              'solar_elevation_angle = {:.1f}'.format(solar_elevation_angle),
              '    ', 'solar_azimuth_angle_landsat = {:.1f}'.format(solar_azimuth_angle_landsat), '    ',
              'solar_azimuth_angle = {:.1f}'.format(solar_azimuth_angle))
        print('\n',
              'slope azimuth METRIC: 0 degree for slopes oriented due south, ±180 degrees due north, -90 degrees due east,'
              ' +90 degrees due west', '\n',
              'solar azimuth Landsat: 0 degree for sun in north, 90 degrees sun in east, '
              '180 degrees sun in south, 270 degrees sun in west')
        print('\n', 'SlopeDeg = {:.0f}'.format(SlopeDeg), '    ', 'AspectDeg = {:.0f}'.format(AspectDeg), '    ',
              'cos_theta_adj = {:.4f}'.format(cos_theta_adj), '    ', 'cos_theta_unadj ={:.4f}'.format(cos_theta_unadj),
              '    ', 'cos_slope ={:.4f}'.format(cos_slope), '    ',
              'cos_theta_horizontal_pixel ={:.4f}'.format(cos_theta_horizontal_pixel))
        print('\n', 'DeclDeg ={:.4f}'.format(DeclDeg), '    ', 'DeclRad ={:.4f}'.format(DeclRad), '    ',
              'HourAngleRad = {:.4f}'.format(HourAngleRad), '    ', 'Sc = {:.4f}'.format(Sc), '    ',
              'bSc = {:.4f}'.format(bSc), '    ', 'a = {:.4f}'.format(a), '    ', 'b = {:.4f}'.format(b), '    ',
              'c = {:.4f}'.format(c))




if __name__ == "__main__":


    # =============== 1-D MODE INPUTS ==========================

    # The solar input variables are: DOY, local_time, Lz and Lm

    # Input variables for Eq. [13-1] in Allen et al. (2010)
    #
    # DOY = day of year
    # AspectRad = aspect in radians: 0 for due south, + or - pi, - pi/2 for east, pi/2 for west
    # Lz = longitude of the center of the local time zone (degrees west of Greenwich
    #    Lz =75, 90, 105, and 120 degrees for Eastern, Central, Rocky Mountain, and Pacific time zones
    # Lm = longitude of the center of the satellite image (degrees west of Greenwich).
    #
    # local_time = local standard clock time for the satellite overpass (hours: 14:30, local t =14.5)

    # ALWAYS SCALAR VALUES
    DOY = 171
    local_time = 12.0
    Lz = -105
    Lm = -105.0

    # SCALAR VALUES ONLY IN 1D
    LatDeg = 35.1056
    SlopeDeg = 80
    AspectDeg = 180

    # =============== RASTER MODE INPUTS =======================

    # We may not need the DEM.
    dem_path = '/Users/dcadol/Desktop/academic_docs_II/m_mountain_metric/Los_Lunas_AOI/dem_clipped/los_lunas_aoi_dem_clipped.tif'

    # RASTER EQUIVALENT TO SlopeDeg
    slope_path = ''
    # RASTER EQUIVALENT TO AspectDeg
    aspect_path = ''
    # RASTER EQUIVALENT TO LatDeg
    latitude_path = ''

    # todo - make it optional to run a raster mode or a 1-D mode.

    calculate_integration_limits(DOY, LatDeg, SlopeDeg, AspectDeg, local_time, Lz, Lm, raster_mode=False)

    # # The way the function would be called to run the function in raster mode
    # calculate_integration_limits(DOY, latitude_path, slope_path, aspect_path, local_time, Lz, Lm, raster_mode=True)