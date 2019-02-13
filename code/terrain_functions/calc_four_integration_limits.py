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
import os
import drigo
import numpy as np
from matplotlib import pyplot as plt
# ============= standard library imports ========================


# TODO - MAKE SURE THAT YOU USE PYMETRIC-WRITTEN FUNCTIONS TO READ THE RASTERS AS ARRAYS DONT REPEAT FUNCTIONS ALREADY WRITTEN BY DRI

# This function applies when Eq [7] from Allen 2006 is met.
def calc_two_integration_limits(LatRad, DeclRad, a, b, c, raster_mode=False):
    # This function is based on Allen (2006)
    # The function is for the northern hemisphere below the arctic circle, i.e. there are no 24 hour days or nights
    # See Appendix A in Allen (2006)

    # STEP A
    # Calculate the sunset and sunrise angles for horizontal surfaces. Eq. [8]
    # Allen (2006) uses often [minus omega_s] for sunrise angle; this function uses the variable omega_rise_hor,
    #     the minus sign is already part of this variable.
    # The sunset and sunrise angles for horizontal pixels calculated with Eq. [8] differ from the sunset and sunrise
    #     angles calculated using Eqs. [13a and 13b]. Eq. [8] calculates the times when the top edge of the sun reaches
    #     the horizon in the morning or when the top edge of the sun has passed the horizon in the evening. Eqs. [13a
    #     13b calculate the time when the sun center reaches the horizon in the morning or when the sun center has
    #     passed the horizon in the evening.
    omega_set_hor = np.arccos(-np.tan(DeclRad) * np.tan(LatRad))
    omega_rise_hor = -np.arccos(-np.tan(DeclRad) * np.tan(LatRad))
    # Calculate cos_theta with Eq. [14] using slope and aspect of a specific pixel to produce cosines of potential
    #     solar incidence angles when the time angles equal those of the time of horizontal sunset and horizontal
    #     sunrise.
    cos_theta_omega_set_hor = - a + b * np.cos(omega_set_hor) + c * np.sin(omega_set_hor)
    cos_theta_omega_rise_hor = - a + b * np.cos(omega_rise_hor) + c * np.sin(omega_rise_hor)
    #

    if not raster_mode:
        print('\n', '\n', 'Output Function calc_integration_limits(LatRad, DeclRad, a, b, c)', '\n', '\n', 'STEP A')
        print(' omega_rise_hor = {:.5f}'.format(omega_rise_hor), '              ',
              'omega_set_hor = {:.5f}'.format(omega_set_hor))
        print(' cos_theta_omega_rise_hor = {:.5f}'.format(cos_theta_omega_rise_hor), '    ',
              'cos_theta_omega_set_hor = {:.5f}'.format(cos_theta_omega_set_hor))

    # STEP D - Section (i)
    #     In Eqs. [13a, b] the expression (b^2 + C^2 - a^2) under the square root sign must be limited to  > 0 for
    #         numerical stability. Therefore, if (b^2 + C^2 - a^2) is 0 or less, it is set equal to 0.001.
    quadratic_function = b**2 + c**2 - a**2
    print('\n', 'STEP D - Section (i)', '\n', 'Before check on positive value', '      ', 'quadratic_function =  {:.6f}'.format(quadratic_function))
    if quadratic_function <= 0.0:
        quadratic_function = 0.001
    print(' After check on positive value', '       ', 'quadratic_function =  {:.6f}'.format(quadratic_function))
    print('\n', 'STEP B: Determine the beginning integration limit at sun rise',
          '\n','\n', 'STEP B - Section (i)')
    #
    # STEP B: Determine the beginning integration limit at sun rise
    # STEP B - Section (i)
    # Calculate the sine of the sunrise angle on a specific pixel using Eq. [13a]
    sin_omega_rise_pixel = (a * c - b * np.sqrt(quadratic_function)) / (b**2 + c**2)

    if not raster_mode:
        print(' Before check on sin values within ±1    sin_omega_rise_pixel = {:.5f}'.format(sin_omega_rise_pixel))

    if sin_omega_rise_pixel < -1.0:
        sin_omega_rise_pixel = -1.0
    if sin_omega_rise_pixel > 1.0:
        sin_omega_rise_pixel = 1.0

    if not raster_mode:
        print(' After check on sin values within ±1     sin_omega_rise_pixel = {:.5f}'.format(sin_omega_rise_pixel))
        print('\n', 'STEP B - Section (ii)')

    # STEP B - Section (ii)
    omega_rise_pixel = np.arcsin(sin_omega_rise_pixel)

    if not raster_mode:
        print(' omega_rise_pixel = {:.5f}'.format(omega_rise_pixel))

        print('\n', 'STEP B - Section (iii)')

    # STEP B - Section (iii)
    # Calculate cosine of theta using omega_rise_pixel using Eq. [14]
    #cos_theta_omega_rise_pixel_eq14 = - a + b * np.cos(omega_rise_pixel) + c * np.sin(omega_rise_pixel)
    cos_theta_omega_rise_pixel = - a + b * np.cos(omega_rise_pixel) + c * np.sin(omega_rise_pixel)

    if not raster_mode:
        print(' Before check on cos_theta_omega_rise_pixel in STEP B - Section (iv)     '
              'cos_theta_omega_rise_pixel = {:.5f}'.format(cos_theta_omega_rise_pixel))

        # STEP B - Section (iv)
        print('\n','STEP B - Section (iv)')
        print(' IF cos_theta_omega_rise_hor <= cos_theta_omega_rise_pixel AND cos_theta_omega_rise_pixel < 0.001')
        print('                {:.9f}'.format(cos_theta_omega_rise_hor),
              '                      {:.9f}'.format(cos_theta_omega_rise_pixel), '                          0.001')

    if cos_theta_omega_rise_hor <= cos_theta_omega_rise_pixel and cos_theta_omega_rise_pixel < 0.001:
        omega_rise_pixel_24 = omega_rise_pixel
        print('\n','IF-statement is TRUE so that omega_rise_pixel_24 = omega_rise_pixel')
        print(' first omega_rise_pixel_24 = {:.5f}'.format(omega_rise_pixel_24))
    else:
        print('\n','IF-statement is FALSE so that a new candidate for omega_rise_pixel needs to be selected')
        print('\n','STEP B - Section (v)')
        print(' The new candidate is called omega_rise_pixel_x = - np.pi - omega_rise_pixel')

    # STEP B - Section (v)
        omega_rise_pixel_x = - np.pi - omega_rise_pixel
        print(' omega_rise_pixel_x = {:.5f}'.format(omega_rise_pixel_x),
              ' omega_rise_pixel = {:.5f}'.format(omega_rise_pixel))
        print(' sin_omega_rise_pixel_x = {:.5f}'.format(np.sin(omega_rise_pixel_x)),
              ' sin_omega_rise_pixel = {:.5f}'.format(sin_omega_rise_pixel))

    # STEP B - Section (v-a)
        cos_theta_omega_rise_pixel_x = - a + b * np.cos(omega_rise_pixel_x) + c * np.sin(omega_rise_pixel_x)
        print(' cos_theta_omega_rise_pixel_x =  {:.5f}'.format(cos_theta_omega_rise_pixel_x),'  ',
              ' omega_rise_hor =  {:.5f}'.format(omega_rise_hor))
        if cos_theta_omega_rise_pixel_x > 0.001:
            omega_rise_pixel_24 = omega_rise_hor
            print(' IF-statement in B_v_a TRUE: omega_rise_pixel_24 =  {:.5f}'.format(omega_rise_pixel_24))
    # Would it not be better to use instead of omega_rise_hor the value for omega_rise_pixel_horizontal, i.e. slope = 0?
    # Section B_v_b
        elif cos_theta_omega_rise_pixel_x <= 0.001 and omega_rise_pixel_x <= omega_rise_hor:
            omega_rise_pixel_24 = omega_rise_hor
            print(' IF-statement in B_v_b TRUE: omega_rise_pixel_24 =  {:.5f}'.format(omega_rise_pixel_24))
        # Section B_v_c
        elif cos_theta_omega_rise_pixel_x <= 0.001 and omega_rise_pixel_x > omega_rise_hor:
            omega_rise_pixel_24 = - np.pi - omega_rise_pixel
            print(' IF-statement in B_v_c TRUE: omega_rise_pixel_24 =  {:.5f}'.format(omega_rise_pixel_24))
    print('check for omega_rise_pixel_24','   ','omega_rise_pixel_24 =  {:.5f}'.format(omega_rise_pixel_24))
    # STEP B - Section (vi)
    if omega_rise_pixel_24 < omega_rise_hor:
        omega_rise_pixel_24 = omega_rise_hor
        print('\n','IF-statement in B_vi TRUE: omega_rise_pixel_24 =  {:.5f}'.format(omega_rise_pixel_24))
    else:
        print('\n','IF-statement in B_vi FALSE: omega_rise_pixel_24 =  {:.5f}'.format(omega_rise_pixel_24))
    # STEP B - Section (vii)
    lower_int_limit_rise = omega_rise_pixel_24
    print ('\n','OUTPUT from STEP B:','     ','  lower_int_limit_rise = {:.5f}'.format(lower_int_limit_rise))

    # STEP C: Determine the ending integration limit at sun set
    # STEP C - Section (i)
    # Calculate the sine of the sunset angle on a specific pixel using Eq. [13b]
    sin_omega_set_pixel = (a * c + b * np.sqrt(quadratic_function)) / (b**2 + c**2)
    #
    print ('\n','\n','STEP C: Determine the ending integration limit at sun set','\n','\n','STEP C - Section (i)','\n',
        'Before check on sin values within ±1 sin_omega_rise_pixel = {:.5f}'.format(sin_omega_set_pixel))
    #
    if sin_omega_set_pixel < -1.0:
        sin_omega_set_pixel = -1.0
    if sin_omega_set_pixel > 1.0:
        sin_omega_set_pixel = 1.0
    #
    print(' After check on sin values within ±1     sin_omega_set_pixel = {:.5f}'.format(sin_omega_set_pixel))
    print('\n', 'STEP C - Section (ii)')
    #
    # STEP C - Section (ii)
    omega_set_pixel = np.arcsin(sin_omega_set_pixel)
    print (' omega_set_pixel = {:.5f}'.format(omega_set_pixel))

    print('\n', 'STEP C - Section (iii)')
    # STEP C - Section (iii)
    # Calculate cosine of theta using omega_set_pixel using Eq. [14]
    cos_theta_omega_set_pixel = - a + b * np.cos(omega_set_pixel) + c * np.sin(omega_set_pixel)
    print(' Before check on cos_theta_omega_set_pixel in STEP C - Section (iv)     '
          'cos_theta_omega_set_pixel = {:.5f}'.format(cos_theta_omega_set_pixel))
    # STEP C - Section (iv)
    print ('\n','STEP C - Section (iv)')
    print (' IF cos_theta_omega_set_hor <= cos_theta_omega_set_pixel AND cos_theta_omega_set_pixel < 0.001')
    print('                {:.9f}'.format(cos_theta_omega_set_hor),
          '                      {:.9f}'.format(cos_theta_omega_set_pixel), '                          0.001')
    if cos_theta_omega_set_hor <= cos_theta_omega_set_pixel and cos_theta_omega_set_pixel < 0.001:
        omega_set_pixel_24 = omega_set_pixel
        print ('\n','IF-statement is TRUE so that omega_set_pixel_24 = omega_set_pixel')
        print (' first omega_set_pixel_24 = {:.5f}'.format(omega_set_pixel_24))
    else:
        print('\n','IF-statement is FALSE so that a new candidate for omega_set_pixel needs to be selected')
        print ('\n','STEP C - Section (v)')
        print(' The new candidate is called omega_set_pixel_x = np.pi - omega_set_pixel')
    # STEP C - Section (v)
        omega_set_pixel_x = np.pi - omega_set_pixel
        print(' omega_set_pixel_x = {:.5f}'.format(omega_set_pixel_x),
              ' omega_set_pixel = {:.5f}'.format(omega_set_pixel))
        print(' sin_omega_set_pixel_x = {:.5f}'.format(np.sin(omega_set_pixel_x)),
              ' sin_omega_set_pixel = {:.5f}'.format(sin_omega_set_pixel))
    # STEP C - Section (v-a)
        cos_theta_omega_set_pixel_x = - a + b * np.cos(omega_set_pixel_x) + c * np.sin(omega_set_pixel_x)
        print(' cos_theta_omega_set_pixel_x =  {:.5f}'.format(cos_theta_omega_set_pixel_x),'  ',
              ' omega_set_hor =  {:.5f}'.format(omega_set_hor))
        if cos_theta_omega_set_pixel_x > 0.001:
            omega_set_pixel_24 = omega_set_hor
            print(' IF-statement in C_v_a TRUE: omega_set_pixel_24 =  {:.5f}'.format(omega_set_pixel_24))
    # Would it not be better to use instead of omega_set_hor the value for omega_set_pixel_horizontal, i.e. slope = 0?
        # Section C_v_b
        elif cos_theta_omega_set_pixel_x <= 0.001 and omega_set_pixel_x >= omega_set_hor:
            omega_set_pixel_24 = omega_set_hor
            print(' IF-statement in C_v_b TRUE: omega_set_pixel_24 =  {:.5f}'.format(omega_set_pixel_24))
        # Section C_v_c
        elif cos_theta_omega_set_pixel_x <= 0.001 and omega_set_pixel_x < omega_set_hor:
            omega_set_pixel_24 = np.pi - omega_set_pixel
            print(' IF-statement in C_v_c TRUE: omega_set_pixel_24 =  {:.5f}'.format(omega_set_pixel_24))
    # STEP C - Section (vi)
    if omega_set_pixel_24 > omega_set_hor:
        omega_set_pixel_24 = omega_set_hor
        print(' IF-statement in C_vi TRUE: omega_set_pixel_24 =  {:.5f}'.format(omega_set_pixel_24))
    else:
        print('\n', 'IF-statement in C_vi FALSE: omega_set_pixel_24 =  {:.5f}'.format(omega_set_pixel_24))
    # STEP C - Section (vii)
    upper_int_limit_set = omega_set_pixel_24
    print('\n', 'OUTPUT from STEP C:', '     ', '  upper_int_limit_set/omega_set_pixel_24 = {:.5f}'.format(upper_int_limit_set),
          '          ','omega_set_hor = {:.5f}'.format(omega_set_hor))
    print('\n', 'OUTPUT from STEP B:', '     ', '  lower_int_limit_rise/omega_rise_pixel_24 = {:.5f}'.format(lower_int_limit_rise),
          '        ','omega_rise_hor = {:.5f}'.format(omega_rise_hor))

    return omega_set_hor, omega_rise_hor, cos_theta_omega_set_hor, cos_theta_omega_rise_hor, sin_omega_rise_pixel, \
           omega_rise_pixel_24, sin_omega_set_pixel, omega_set_pixel_24, lower_int_limit_rise, upper_int_limit_set,\
           quadratic_function


def calc_solar_topo_image_variables(DOY, LatDeg, SlopeDeg, AspectDeg, local_time, Lz, Lm ):
    """"""
    # This function calculates the cosine and sine for variables that remain constant for a given Landsat image.
    #
    # The solar variables are: LatDeg (latitude), DOY (day of year)
    # The topographical variables are SlopeDeg (slope of a pixel) and AspectDeg (aspect of a pixel)

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
        # read in each raster parameter as a numpy array and continue
        LatDeg = drigo.raster_to_array(LatDeg)
        SlopeDeg = drigo.raster_to_array(SlopeDeg)
        AspectDeg = drigo.raster_to_array(AspectDeg)



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

    if raster_mode:
        print( 'NON RASTER parameters\n')
        print('\n',
              'Output Function calc_solar_topo_image_variables(DOY, LatDeg, SlopeDeg, AspectDeg, local_time, Lz, Lm)')
        print('\n', 'DOY = {:.0f}'.format(DOY),
              '   ', 'local_time = {:.2f}'.format(local_time), '  ', 'solar_time = {:.2f}'.format(solar_time), '  '
              , 'Lz = {:.0f}'.format(Lz), '  ', 'Lm = {:.0f}'.format(Lm))

        print('LAT RAD\n')
        plt.imshow(LatRad)
        plt.show()

        print('SOLAR ANGLES\n')
        print('solar incidence angle theta')
        plt.imshow(solar_incidence_angle_theta)
        plt.show()

        print('solar elevation angle')
        plt.imshow(solar_elevation_angle)
        plt.show()

        print('solar azimuth angle landsat')
        plt.imshow(solar_azimuth_angle_landsat)
        plt.show()

        print('solar azimuth angle')
        plt.imshow(solar_azimuth_angle)
        plt.show()

        print('\n',
              'slope azimuth METRIC: 0 degree for slopes oriented due south, ±180 degrees due north, -90 degrees due east,'
              ' +90 degrees due west', '\n',
              'solar azimuth Landsat: 0 degree for sun in north, 90 degrees sun in east, '
              '180 degrees sun in south, 270 degrees sun in west')

        print('ALLEN 2010 Slope and cos theta \n')

        print('slope degrees')
        plt.imshow(SlopeDeg)
        plt.show()

        print('aspect in degrees')
        plt.imshow(AspectDeg)
        plt.show()

        print(' cosine theta adjusted (ratio)')
        plt.imshow(cos_theta_adj)
        plt.show()

        print('cosine theta unadjusted')
        plt.imshow(cos_theta_unadj)
        plt.show()

        print('cosine of slope')
        plt.imshow(cos_slope)
        plt.show()

        print('coside theta of hirzontal pixel')
        plt.imshow(cos_theta_horizontal_pixel)
        plt.show()

        print('DECLINATIONs and HOUR ANGLES')

        print('Degrees of declination')
        plt.imshow(DeclDeg)
        plt.show()

        print('Declination in Radians')
        plt.imshow(DeclRad)
        plt.show()

        print('hour angle in radians')
        plt.imshow(HourAngleRad)
        plt.show()

        print('Sc')
        plt.imshow(Sc)
        plt.show()

        print('bSc')
        plt.imshow(bSc)
        plt.show()

        print(' PARAMS a, b and c \n')
        plt.imshow(a)
        plt.show()
        plt.imshow(b)
        plt.show()
        plt.imshow(c)
        plt.show()


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
    slope_path = '/Users/dcadol/Desktop/academic_docs_II/m_mountain_metric/Los_Lunas_AOI/dem_clipped/los_lunas_aoi_dem_clipped_slope.tif'
    # RASTER EQUIVALENT TO AspectDeg
    aspect_path = '/Users/dcadol/Desktop/academic_docs_II/m_mountain_metric/Los_Lunas_AOI/dem_clipped/los_lunas_aoi_dem_clipped_aspect.tif'
    # RASTER EQUIVALENT TO LatDeg
    latitude_path = '/Users/dcadol/Desktop/academic_docs_II/m_mountain_metric/Los_Lunas_AOI/dem_clipped/los_lunas_aoi_dem_clipped_lat.tif'

    # todo - make it optional to run a raster mode or a 1-D mode.

    # calculate_integration_limits(DOY, LatDeg, SlopeDeg, AspectDeg, local_time, Lz, Lm, raster_mode=False)

    # # The way the function would be called to run the function in raster mode
    calculate_integration_limits(DOY, latitude_path, slope_path, aspect_path, local_time, Lz, Lm, raster_mode=True)