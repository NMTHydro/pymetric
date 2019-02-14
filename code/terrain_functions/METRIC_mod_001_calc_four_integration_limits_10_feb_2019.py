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
import numpy as np
# ============= standard library imports ========================
#
# Input variables for Eq. [13-1] in Allen et al. (2010) and the similar Eq. [3] in Allen et al. (2006)
# Allen (2010) stands for Allen et al. (2010) METRIC Manual
# Allen (2006) for Allen et al. (2006), Agr and Forest Met. 139: 55-73
# The terrain input variables are: SlopeDeg, AspectDeg, LatDeg
# The solar input variables are: DOY, local_time, Lz and Lm
# DOY = day of year
# AspectRad = aspect in radians: 0 for due south, + or - pi for north, - pi/2 for east, pi/2 for west
# Lz = longitude of the center of the local time zone (degrees west of Greenwich
#    Lz =75, 90, 105, and 120 degrees for Eastern, Central, Rocky Mountain, and Pacific time zones
# Lm = longitude of the center of the satellite image (degrees west of Greenwich).
#
# local_time = local standard clock time for the satellite overpass (hours: 14:30, local t =14.5)
#
# Input variables for Eq. [13-1] in Allen et al. (2010)
#
# DOY = day of year
# AspectRad = aspect in radians: 0 for due south, + or - pi, - pi/2 for east, pi/2 for west
# Lz = longitude of the center of the local time zone (degrees west of Greenwich
#    Lz =75, 90, 105, and 120 degrees for Eastern, Central, Rocky Mountain, and Pacific time zones
# Lm = longitude of the center of the satellite image (degrees west of Greenwich).
#
# local_time = local standard clock time for the satellite overpass (hours: 14:30, local t =14.5)
#
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

DOY = 171
local_time = 12.0
LatDeg = 35.1056
Lz = -105
Lm = -105.0
SlopeDeg = 80
AspectDeg = 180
#
def calc_solar_topo_image_variables(DOY, LatDeg, SlopeDeg, AspectDeg, local_time, Lz, Lm ):
#
# This function calculates the cosine and sine for variables that remain constant for a give Landsat image.
#
# The solar variables are: LatDeg (latitude), DOY (day of year)
# The topographical variables are SlopeDeg (slope of a pixel) and AspectDeg (aspect of a pixel)
#
# The variables a, b, and c are needed for the integration of cos_theta over a day or part of a day
    LatRad = LatDeg * np.pi / 180
    DeclRad = 0.409 * np.sin((2 * np.pi / 365 * DOY) - 1.39)
    DeclDeg = DeclRad * 180 / np.pi
    SlopeRad = SlopeDeg * np.pi / 180
    AspectRad = AspectDeg * np.pi / 180
    cos_lat = np.cos(LatRad)
    sin_lat = np.sin(LatRad)
    cos_decl = np.cos(DeclRad)
    sin_decl = np.sin(DeclRad)
    cos_slope = np.cos(SlopeRad)
    sin_slope = np.sin(SlopeRad)
    cos_aspect = np.cos(AspectRad)
    sin_aspect = np.sin(AspectRad)
    bSc = 2 * np.pi * (DOY - 81) / 364
    Sc = 0.1645 * np.sin(2 * bSc) - 0.1255 * np.cos(bSc) - 0.025 * np.sin(bSc)
    solar_time = (local_time + 0.6667 * (Lz - Lm) + Sc)
    HourAngleRad = np.pi / 12 * ((local_time + 0.6667 * (Lz - Lm) + Sc) - 12)
    cos_hourangle = np.cos(HourAngleRad)
    sin_hourangle = np.sin(HourAngleRad)
    a = sin_decl * cos_lat * sin_slope * cos_aspect - sin_decl * sin_lat * cos_slope
    b = cos_decl * cos_lat * cos_slope + cos_decl * sin_lat * sin_slope * cos_aspect
    c = cos_decl * sin_aspect * sin_slope
    cos_theta_unadj = - a + b * cos_hourangle + c * sin_hourangle
    cos_theta_adj = cos_theta_unadj / cos_slope
    if cos_theta_adj < 0.1:
        cos_theta_adj = 0.1
    if cos_theta_adj > 10.0:
        cos_theta_adj = 10.0
#
# cos_theta_horizontal_pixel at satellite overpass with Eq.[4] in Allen (2006)
# solar incidence angle at satellite overpass on horizontal pixel is solar_incidence_angle_theta
# solar elevation angle at satellite  overpass on horizontal pixel is solar_elevation_angle
# solar azimuth angle at satellite overpass is solar_azimuth_angle with Eq. 11.5 in Campbell and Norman (1998)
#    where solar_azimuth_angle is calculated with respect to due south, increasing in the counter clockwise direction
#    so 90 degrees is east. Afternoon azimuth angles can be calculated by taking 360 degrees minus the solar_azimuth angle
# solar_azimuth_angle_Landsat (also used by NOAA) is measured clockwise starting at 0 in the north.
#
    cos_theta_horizontal_pixel = sin_decl * sin_lat + cos_decl * cos_lat * cos_hourangle
    solar_incidence_angle_theta = np.arccos(cos_theta_horizontal_pixel) * 180 / np.pi
    solar_elevation_angle = 90.0 - solar_incidence_angle_theta
    cos_solar_azimuth_angle = -(sin_decl - cos_theta_horizontal_pixel * sin_lat) / \
                              (cos_lat * np.sin(solar_incidence_angle_theta * np.pi / 180))
    solar_azimuth_angle = np.arccos(cos_solar_azimuth_angle) * 180 / np.pi
    if HourAngleRad <= 0.000:
        solar_azimuth_angle_landsat = 180 - solar_azimuth_angle
    solar_azimuth_angle_landsat = 180 + solar_azimuth_angle

    return LatRad, DeclRad, DeclDeg, SlopeRad, AspectRad, cos_lat, sin_lat, cos_slope, sin_slope, \
    cos_aspect, sin_aspect, sin_decl, cos_decl, bSc, Sc, HourAngleRad, cos_hourangle, sin_hourangle, a, b, c, cos_theta_unadj, \
    cos_theta_adj, cos_theta_horizontal_pixel, solar_incidence_angle_theta, \
           solar_elevation_angle, solar_azimuth_angle, solar_azimuth_angle_landsat, solar_time

LatRad, DeclRad, DeclDeg, SlopeRad, AspectRad, cos_lat, sin_lat, cos_slope, sin_slope, \
    cos_aspect, sin_aspect, sin_decl, cos_decl, bSc, Sc, HourAngleRad, cos_hourangle, sin_hourangle, a, b, c, cos_theta_unadj, \
    cos_theta_adj, cos_theta_horizontal_pixel, solar_incidence_angle_theta, \
           solar_elevation_angle, solar_azimuth_angle, solar_azimuth_angle_landsat, solar_time = \
    calc_solar_topo_image_variables(DOY, LatDeg, SlopeDeg, AspectDeg, local_time, Lz, Lm)

print('\n','Output Function calc_solar_topo_image_variables(DOY, LatDeg, SlopeDeg, AspectDeg, local_time, Lz, Lm)')
print('\n','DOY = {:.0f}'.format(DOY),'    ','LatDeg = {:.2f}'.format(LatDeg),'   ','LatRad = {:.4f}'.format(LatRad),
      '   ','local_time = {:.2f}'.format(local_time), '  ','solar_time = {:.2f}'.format(solar_time), '  '
      ,'Lz = {:.0f}'.format(Lz), '  ','Lm = {:.0f}'.format(Lm))
print ('\n','solar_incidence_angle_theta = {:.1f}'.format(solar_incidence_angle_theta), \
           'solar_elevation_angle = {:.1f}'.format(solar_elevation_angle),
           '    ','solar_azimuth_angle_landsat = {:.1f}'.format(solar_azimuth_angle_landsat),'    ',
           'solar_azimuth_angle = {:.1f}'.format(solar_azimuth_angle))
print ('\n','slope azimuth METRIC: 0 degree for slopes oriented due south, ±180 degrees due north, -90 degrees due east,'
            ' +90 degrees due west','\n','solar azimuth Landsat: 0 degree for sun in north, 90 degrees sun in east, '
                                         '180 degrees sun in south, 270 degrees sun in west')
print ('\n','SlopeDeg = {:.0f}'.format(SlopeDeg),'    ','AspectDeg = {:.0f}'.format(AspectDeg),'    ',
       'cos_theta_adj = {:.4f}'.format(cos_theta_adj),'    ','cos_theta_unadj ={:.4f}'.format(cos_theta_unadj),
       '    ','cos_slope ={:.4f}'.format(cos_slope),'    ',
       'cos_theta_horizontal_pixel ={:.4f}'.format(cos_theta_horizontal_pixel))
print ('\n', 'DeclDeg ={:.4f}'.format(DeclDeg),'    ','DeclRad ={:.4f}'.format(DeclRad),'    ',
       'HourAngleRad = {:.4f}'.format(HourAngleRad),'    ',        'Sc = {:.4f}'.format(Sc),'    ',
       'bSc = {:.4f}'.format(bSc),'    ', 'a = {:.4f}'.format(a),'    ','b = {:.4f}'.format(b),'    ',
       'c = {:.4f}'.format(c))


def calc_two_integration_limits(LatRad, DeclRad, a, b, c):
    # This function is based on Allen (2006)
    # The function is for the northern hemisphere below the arctic circle, i.e. there are no 24 hour days or nights
    # See Appendix A in Allen (2006)
    #
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
    print('\n', '\n', 'Output Function calc_integration_limits(LatRad, DeclRad, a, b, c)', '\n', '\n', 'STEP A')
    print(' omega_rise_hor = {:.5f}'.format(omega_rise_hor), '              ',
          'omega_set_hor = {:.5f}'.format(omega_set_hor))
    print(' cos_theta_omega_rise_hor = {:.5f}'.format(cos_theta_omega_rise_hor), '    ',
          'cos_theta_omega_set_hor = {:.5f}'.format(cos_theta_omega_set_hor))
    #
    # STEP D - Section (i)
    #     In Eqs. [13a, b] the expression (b^2 + C^2 - a^2) under the square root sign must be limited to  > 0 for
    #         numerical stability. Therefore, if (b^2 + C^2 - a^2) is 0 or less, it is set equal to 0.001.
    quadratic_function = b**2 + c**2 - a**2
    print ('\n','STEP D - Section (i)','\n','Before check on positive value','      ','quadratic_function =  {:.6f}'.format(quadratic_function))
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
    #
    print (' Before check on sin values within ±1    sin_omega_rise_pixel = {:.5f}'.format(sin_omega_rise_pixel))
    #
    if sin_omega_rise_pixel < -1.0:
        sin_omega_rise_pixel = -1.0
    if sin_omega_rise_pixel > 1.0:
        sin_omega_rise_pixel = 1.0
    #
    print(' After check on sin values within ±1     sin_omega_rise_pixel = {:.5f}'.format(sin_omega_rise_pixel))
    print('\n', 'STEP B - Section (ii)')
    #
    # STEP B - Section (ii)
    omega_rise_pixel = np.arcsin(sin_omega_rise_pixel)
    print(' omega_rise_pixel = {:.5f}'.format(omega_rise_pixel))

    print('\n', 'STEP B - Section (iii)')
    # STEP B - Section (iii)
    # Calculate cosine of theta using omega_rise_pixel using Eq. [14]
    #cos_theta_omega_rise_pixel_eq14 = - a + b * np.cos(omega_rise_pixel) + c * np.sin(omega_rise_pixel)
    cos_theta_omega_rise_pixel = - a + b * np.cos(omega_rise_pixel) + c * np.sin(omega_rise_pixel)
    print(' Before check on cos_theta_omega_rise_pixel in STEP B - Section (iv)     '
          'cos_theta_omega_rise_pixel = {:.5f}'.format(cos_theta_omega_rise_pixel))
    # STEP B - Section (iv)
    print ('\n','STEP B - Section (iv)')
    print (' IF cos_theta_omega_rise_hor <= cos_theta_omega_rise_pixel AND cos_theta_omega_rise_pixel < 0.001')
    print('                {:.9f}'.format(cos_theta_omega_rise_hor),
          '                      {:.9f}'.format(cos_theta_omega_rise_pixel), '                          0.001')
    if cos_theta_omega_rise_hor <= cos_theta_omega_rise_pixel and cos_theta_omega_rise_pixel < 0.001:
        omega_rise_pixel_24 = omega_rise_pixel
        print ('\n','IF-statement is TRUE so that omega_rise_pixel_24 = omega_rise_pixel')
        print (' first omega_rise_pixel_24 = {:.5f}'.format(omega_rise_pixel_24))
    else:
        print('\n','IF-statement is FALSE so that a new candidate for omega_rise_pixel needs to be selected')
        print ('\n','STEP B - Section (v)')
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

omega_set_hor, omega_rise_hor, cos_theta_omega_set_hor, cos_theta_omega_rise_hor, sin_omega_rise_pixel, \
omega_rise_pixel_24, sin_omega_set_pixel, \
omega_set_pixel_24, lower_int_limit_rise, upper_int_limit_set, quadratic_function  = \
    calc_two_integration_limits(LatRad, DeclRad, a, b, c)
#
#
def calc_two_daytime_integration_limits(omega_rise_pixel_24, omega_set_pixel_24, lower_int_limit_rise,
                                        upper_int_limit_set, sin_slope, sin_lat, sin_decl, sin_aspect, cos_slope,
                                        cos_lat, cos_decl, cos_aspect, a, b, c, quadratic_function):
# STEP D: Additional limits on omega_rise_pixel_24 and omega_set_pixel_24 for numerical stability and
#    twice per day periods of sun
#
# STEP D - Section (i)
# The argument of Eqs. [13a and 13b] cannot be equal to zero for numerical stability. This issue has been taken care of
#   under STEP B and STEP C.

# STEP D - Section (ii)
#
# The value of omega_rise_pixel_24 must be smaller or equal to omega_set_pixel_24 (i.e. sunrise occurs before sunset).
#
    if omega_set_pixel_24 < omega_rise_pixel_24:
        omega_rise_pixel_24 = omega_set_pixel_24
        print ('\n','STEP D - Section (ii)','\n','IF-statement in D_ii is TRUE:','\n','Slope is always shaded because '
                                                         'omega_rise_pixel_24 equals omega_set_pixel_24')
        print ('\n', 'omega_rise_pixel_24 = {:.5f}'.format(omega_rise_pixel_24),'    ',
           'omega_set_pixel_24 = {:.5f}'.format(omega_set_pixel_24))
    else:
        print('\n','STEP D - Section (ii)','\n','Check that omega_rise_pixel_24 is smaller than omega_set_pixel_24','\n',
              'omega_rise_pixel_24 = {:.5f}'.format(omega_rise_pixel_24), '  <  ',
          'omega_set_pixel_24 = {:.5f}'.format(omega_set_pixel_24))
#
# This means that there is no direct solar beam during the day; the slope is always shaded.
#
# STEP D - Section (iii)
# The sine values in Eqs. [13a and 13b] cannot be smaller than -1 or larger than +1. This issue has been taken care of
#   under STEP B and STEP C.
#
# STEP D - Section (iv)
#    Check if possibility exists for two periods of direct beam radiation during the day.
#    The two daytime integration limits are defined as follows:
#       omega_set_during_day_pixel_24 = the time angle when the center of the solar disk disappears the first time
#       omega_rise_during_day_pixel_24 = the time angle when the center of the solar disk reappears over the surface
#
    print('\n','STEP D - Section (iv)','\n','Check IF-statement sin_slope > (sin_lat * cos_decl + cos_lat * sin_decl)')
#    temp = sin_lat * cos_decl + cos_lat * sin_decl
    if sin_slope <= sin_lat * cos_decl + cos_lat * sin_decl:
        print ('\n','sin_slope = {:.5f}'.format(sin_slope),'  <   ',
               'sin_lat * cos_decl + cos_lat * sin_decl = {:.5f}'.format(sin_lat * cos_decl + cos_lat * sin_decl))
        omega_set_during_day_pixel_24 = 0.000000
        omega_rise_during_day_pixel_24 = 0.000000
        print(' There is only one period of direct beam radiation during the day', '\n',
              'No daytime integration limits exist and therefore:','',
              'omega_set_during_day_pixel_24 = {:.5f}'.format(omega_set_during_day_pixel_24),'    ',
              'omega_rise_during_day_pixel_24 = {:.5f}'.format(omega_rise_during_day_pixel_24))
    else:
        print ('\n','sin_slope = {:.5f}'.format(sin_slope),'  >   ',
               'sin_lat * cos_decl + cos_lat * sin_decl = {:.5f}'.format(sin_lat * cos_decl + cos_lat * sin_decl))
        print (' Possibility exists for two periods of direct beam radiation during the day', '\n','\n',
               'STEP D - Section (iv-b, c)')
#
#    STEP D - Section (iv-a)
#    omega_rise_pixel_24 and omega_set_pixel_24 have already been calculated in STEP B and STEP C
#
#    STEP D - Section (iv-b, c)
#    The candidates for intermediate integration limits omega_set_during_day_pixel_24 and omega_rise_during_day_pixel_24
#        are
    sin_A = (a * c + b * np.sqrt(quadratic_function)) / (b**2 + c**2)
    sin_B = (a * c - b * np.sqrt(quadratic_function)) / (b**2 + c**2)
#
#     STEP D - Section(iv - c)
#
    if sin_A < -1.0:
            sin_A = -1.0
    if sin_A > 1.0:
            sin_A = 1.0
    if sin_B < -1.0:
            sin_B = -1.0
    if sin_B > 1.0:
            sin_B = 1.0
    print(' sin_A = {:.5f}'.format(sin_A),'      ','sin_B = {:.5f}'.format(sin_B))
    A = np.arcsin(sin_A)
    B = np.arcsin(sin_B)
    print (' A = {:.5f}'.format(A),'          ','B = {:.5f}'.format(B),'\n','\n',
              'STEP D - Section(iv - d)')
#
#     STEP D - Section(iv - d)
#
    omega_set_during_day_pixel_24 = min(A,B)
    omega_rise_during_day_pixel_24 = max(A,B)
    print('\n', 'omega_set_during_day_pixel_24 = {:.5f}'.format(omega_set_during_day_pixel_24), '      ',
              'omega_rise_during_day_pixel_24 = {:.5f}'.format(omega_rise_during_day_pixel_24))
#
#     STEP D - Section(iv - e)
#
    cos_theta_omega_set_during_day_pixel_24 = - a + b * np.cos(omega_set_during_day_pixel_24) + c \
                                                  * np.sin(omega_set_during_day_pixel_24)
    cos_theta_omega_rise_during_day_pixel_24 = - a + b * np.cos(omega_rise_during_day_pixel_24) + c \
                                                  * np.sin(omega_rise_during_day_pixel_24)
    print('\n', 'cos_theta_omega_set_during_day_pixel_24 = {:.5f}'.format(cos_theta_omega_set_during_day_pixel_24),
              '      ',
              'cos_theta_omega_rise_during_day_pixel_24 = {:.5f}'.format(cos_theta_omega_rise_during_day_pixel_24))
#
    if cos_theta_omega_set_during_day_pixel_24 < -0.001 or cos_theta_omega_set_during_day_pixel_24 > 0.001:
            omega_set_during_day_pixel_24 = - np.pi - omega_set_during_day_pixel_24
    if cos_theta_omega_rise_during_day_pixel_24 < -0.001 or cos_theta_omega_rise_during_day_pixel_24 > 0.001:
            omega_rise_during_day_pixel_24 = np.pi - omega_rise_during_day_pixel_24
    print('\n', 'omega_set_during_day_pixel_24 = {:.5f}'.format(omega_set_during_day_pixel_24),
              '      ',
              'omega_rise_during_day_pixel_24 = {:.5f}'.format(omega_rise_during_day_pixel_24))
    print('\n', 'omega_set_pixel_24 = {:.5f}'.format(omega_set_pixel_24),
              '      ',
              'omega_rise_pixel_24 = {:.5f}'.format(omega_rise_pixel_24))
#
#     STEP D - Section(iv - f and g)
#
#        if omega_set_during_day_pixel_24 < omega_rise_pixel_24 or omega_rise_during_day_pixel_24 > omega_set_pixel_24:
#            print('\n',' There is only one period of direct beam radiation during the day', '\n',' because','\n',
#                  ' omega_set_during_day_pixel_24 < omega_rise_pixel_24 or '
#                  'omega_rise_during_day_pixel_24 < omega_set_pixel_24')
#            print('  omega_set_during_day_pixel_24 = {:.5f}'.format(omega_set_during_day_pixel_24),
#            'omega_rise_pixel_24 = {:.5f}'.format(omega_rise_pixel_24),
#                  'omega_rise_during_day_pixel_24 = {:.5f}'.format(omega_rise_during_day_pixel_24),
#                  ' omega_set_pixel_24 = {:.5f}'.format(omega_set_pixel_24))
#            omega_set_during_day_pixel_24 = 0.000000
#            omega_rise_during_day_pixel_24 = 0.000000
#        else:
#            X = sin_decl * sin_lat * cos_slope * (omega_rise_during_day_pixel_24 - omega_set_during_day_pixel_24)
#            - sin_decl * cos_lat * sin_slope * cos_aspect * (omega_rise_during_day_pixel_24
#                                                             - omega_set_during_day_pixel_24)
#            + cos_decl * cos_lat * cos_slope * (np.sin(omega_rise_during_day_pixel_24)
#                                                - np.sin(omega_set_during_day_pixel_24))
#            + cos_decl * sin_lat * sin_slope * cos_aspect * (np.sin(omega_rise_during_day_pixel_24)
#                                                             - np.sin(omega_set_during_day_pixel_24))
#            - cos_decl * sin_slope * sin_aspect * (np.cos(omega_rise_during_day_pixel_24)
#                                                     - np.cos(omega_set_during_day_pixel_24))
#            if X < 0:
#                print('\n','x = {:.5f}'.format(X),'there are two periods of direct beam radiation during the day')
#            else:
#                print('\n','x = {:.5f}'.format(X),'there is only one period of direct beam radiation during the day')
#
    if omega_set_during_day_pixel_24 < omega_rise_pixel_24:
            omega_set_during_day_pixel_24 = omega_rise_pixel_24
    if omega_rise_during_day_pixel_24 > omega_set_pixel_24:
            omega_rise_during_day_pixel_24 = omega_set_pixel_24


    print('\n', 'omega_set_during_day_pixel_24 = {:.5f}'.format(omega_set_during_day_pixel_24), '      ',
              'omega_rise_during_day_pixel_24 = {:.5f}'.format(omega_rise_during_day_pixel_24))
    print('\n', 'omega_set_pixel_24 = {:.5f}'.format(omega_set_pixel_24), '      ',
              'omega_rise_pixel_24 = {:.5f}'.format(omega_rise_pixel_24))
#
    X = sin_decl * sin_lat * cos_slope * (omega_rise_during_day_pixel_24 - omega_set_during_day_pixel_24) \
        - sin_decl * cos_lat * sin_slope * cos_aspect * (omega_rise_during_day_pixel_24 - omega_set_during_day_pixel_24) \
        + cos_decl * cos_lat * cos_slope * (np.sin(omega_rise_during_day_pixel_24)
                                            - np.sin(omega_set_during_day_pixel_24)) + cos_decl * sin_lat * sin_slope \
        * cos_aspect * (np.sin(omega_rise_during_day_pixel_24) - np.sin(omega_set_during_day_pixel_24)) \
        - cos_decl * sin_slope * sin_aspect * (np.cos(omega_rise_during_day_pixel_24)
                                               - np.cos(omega_set_during_day_pixel_24))
    if X < 0:
        print('\n','x = {:.5f}'.format(X),'there are two periods of direct beam radiation during the day')
    else:
        print('\n','x = {:.5f}'.format(X),'there is only one period of direct beam radiation during the day')

    return lower_int_limit_rise, upper_int_limit_set, omega_rise_pixel_24, omega_set_pixel_24, \
           omega_set_during_day_pixel_24, omega_rise_during_day_pixel_24, X

omega_rise_pixel_24, omega_set_pixel_24, lower_int_limit_rise, upper_int_limit_set, \
omega_set_during_day_pixel_24, omega_rise_during_day_pixel_24, X  \
    = calc_two_daytime_integration_limits(omega_rise_pixel_24, omega_set_pixel_24, lower_int_limit_rise,
                                        upper_int_limit_set, sin_slope, sin_lat, sin_decl, sin_aspect, cos_slope,
                                        cos_lat, cos_decl, cos_aspect, a, b, c, quadratic_function)
#
print ('lower_int_limit_rise = {:.4f}'.format(lower_int_limit_rise))

#
def integral_cos_theta(lower_limit, upper_limit, sin_decl, cos_decl,sin_lat, cos_lat, sin_slope, cos_slope,
                       sin_aspect, cos_aspect):
#   Equation [5] in Allen (2006)
#   lower_limit = lower_limit of integral as hourangle in radisns
#   upper_limit = upper limit of integral as hourangle in radians
#
    int_cos_theta = sin_decl * sin_lat * cos_slope * (upper_limit - lower_limit) \
        - sin_decl * cos_lat * sin_slope * cos_aspect * (upper_limit - lower_limit) \
        + cos_decl * cos_lat * cos_slope * (np.sin(upper_limit) - np.sin(lower_limit)) \
        + cos_decl * sin_lat * sin_slope * cos_aspect * (np.sin(upper_limit) - np.sin(lower_limit)) \
        - cos_decl * sin_slope * sin_aspect * (np.cos(upper_limit) - np.cos(lower_limit))
    print ('\n','lower_limit = {:.5f}'.format(lower_limit))
    print ('\n','upper_limit = {:.5f}'.format(upper_limit))
    print ('\n','upper limit - lower_limit = {:.5f}'.format(upper_limit - lower_limit))
    return int_cos_theta
#
int_cos_theta = integral_cos_theta(omega_rise_hor, omega_set_hor, sin_decl, cos_decl,sin_lat, cos_lat,
                                   sin_slope, cos_slope, sin_aspect, cos_aspect)

Ra = 1367 * ( 1 / (1 + 0.033 * np.cos(DOY * 2 * np.pi / 365))) * int_cos_theta

print ('\n','int_cos_theta = {:.5f}'.format(int_cos_theta))
print ('\n','omega_rise_hor = {:.5f}'.format(omega_rise_hor))
print ('\n','omega_set_hor = {:.5f}'.format(omega_set_hor))
print ('\n','Ra = {:.5f}'.format(Ra))
omega_rise_hor_time_decimal = omega_rise_hor * 180 / np.pi / 15 +12
ihours = int(omega_rise_hor_time_decimal)
omega_rise_hor_time = ihours,(omega_rise_hor_time_decimal - ihours) * 60
print ('%02d:%02d' % omega_rise_hor_time)

omega_set_hor_time_decimal = omega_set_hor * 180 / np.pi / 15 +12
ihours = int(omega_set_hor_time_decimal)
omega_set_hor_time = ihours,(omega_set_hor_time_decimal - ihours) * 60
print ('%02d:%02d' % omega_set_hor_time)

lower_int_limit_rise_time_decimal = lower_int_limit_rise * 180 / np.pi / 15 +12
ihours = int(lower_int_limit_rise_time_decimal)
lower_int_limit_rise_time = ihours,(lower_int_limit_rise_time_decimal - ihours) * 60
print ('%02d:%02d' % lower_int_limit_rise_time)

upper_int_limit_set_time_decimal = upper_int_limit_set * 180 / np.pi / 15 +12
ihours = int(upper_int_limit_set_time_decimal)
upper_int_limit_set_time = ihours,(upper_int_limit_set_time_decimal - ihours) * 60
print ('%02d:%02d' % upper_int_limit_set_time)

omega_set_during_day_pixel_24_time_decimal = omega_set_during_day_pixel_24 * 180 / np.pi / 15 +12
ihours = int(omega_set_during_day_pixel_24_time_decimal)
omega_set_during_day_pixel_24_time = ihours,(omega_set_during_day_pixel_24_time_decimal - ihours) * 60
print ('%02d:%02d' % omega_set_during_day_pixel_24_time)

omega_rise_during_day_pixel_24_time_decimal = omega_rise_during_day_pixel_24 * 180 / np.pi / 15 +12
ihours = int(omega_rise_during_day_pixel_24_time_decimal)
omega_rise_during_day_pixel_24_time = ihours,(omega_rise_during_day_pixel_24_time_decimal - ihours) * 60
print ('%02d:%02d' % omega_rise_during_day_pixel_24_time)
#
#
# Summary of output data for calculation of self-shadowing periods of pixel
print ('\n','Summary of output data for calculation of self-shadowing periods of pixel')
print ('\n','Latitude = {:.4f}'.format(LatDeg),'degrees','    ','DOY = {:.0f}'.format(DOY),'    ','slope = {:.0f}'.format(SlopeDeg),'degrees','     ',
       'aspect = {:.0f}'.format(AspectDeg),'degrees')
print ('\n','omega_rise_hor = {:.4f}'.format(omega_rise_hor),'    ','Daylight starts at','  ',
       '%02d:%02d' % omega_rise_hor_time, '\n','omega_set_hor  =  {:.4f}'.format(omega_set_hor),'    ',
       'Daylight ends at','    ','%02d:%02d' % omega_set_hor_time,)
if X < 0:
    print('\n','X = {:.5f}'.format(X),'at X<0 there are two periods of direct beam radiation during the day')
    print ('\n','lower_int_limit_rise = {:.4f}'.format(lower_int_limit_rise),'    ','Beam radiation starts at','  ',
       '%02d:%02d' % lower_int_limit_rise_time, '\n','upper_int_limit_set  =  {:.4f}'.format(upper_int_limit_set),'    ',
       'Beam radiation ends at','    ','%02d:%02d' % upper_int_limit_set_time,)
    print ('\n','omega_set_during_day_pixel_24   = {:.4f}'.format(omega_set_during_day_pixel_24),'      ',
       'Beam radiation ends during the day at','  ', '%02d:%02d' % omega_set_during_day_pixel_24_time, '\n',
       'omega_rise_during_day_pixel_24  =  {:.4f}'.format(omega_rise_during_day_pixel_24),'      ',
       'Beam radiation starts again at','         ','%02d:%02d' % omega_rise_during_day_pixel_24_time,)
else:
    print('\n','X = {:.5f}'.format(X),'at X=>0 there is only one period of direct beam radiation during the day')
    omega_set_during_day_pixel_24 = 9.0
    omega_rise_during_day_pixel_24 = 9.0
#
print ('\n','the four integration limits for export are')
print('\n', 'lower_int_limit_rise          = {:.4f}'.format(lower_int_limit_rise), '            ',
      ' upper_int_limit_set  =  {:.4f}'.format(upper_int_limit_set))
print(' omega_set_during_day_pixel_24 = {:.4f}'.format( omega_set_during_day_pixel_24), '             ',
      'omega_rise_during_day_pixel_24  =  {:.4f}'.format(omega_rise_during_day_pixel_24))

# Gabe you can ignore everything below this line
#
# minimum for Ra and cos_theta
# sin(d)sin(l)cos(s)-sin(d)cos(l)sin(s)cos(a)
# cos(d)cos(l)cos(s)+cos(d)sin(l)sin(s)cos(a)
# cos(d)sin(s)sin(a)
# Binary "dummy" layer(s)
# sun hour angle w/ 30min int.
# Ra24(W/m2), H.E.
# Image Time
# cos_theta(inst.), H.E.
# (in solar time, ex. 10.45)
# (0.1 is recommended)
# (ex.43.2)
# tempor. Ra 24-test
# DEM
# TauB(inst)
# ea(inst, kPa)
# Kt(inst)
# Rso(inst)
# KB(24)
# Rso(24)pixel
# Rso(inst)flat
# Rso(24)flat
# Crad(max)=1.5 for keeping numerical stability
# Rad. adj. coeff.
# Recommended Kt is "1" for clean air conditions
# P
# TauD(inst)
# KD(24)
# Transmittance
# M001, Pre-Calculations for METRIC 2010 Mountain Model:  --Populated by VBscript 6/11/2013 at 12:31:16 PM
# Copyright (C) 2003-2010 R.G.Allen, M.Tasumi, R.Trezza, J. Kjaersgaard and University of Idaho. All rights reserved.
# Extraterrestrial and clear sky solar radiation for instantaneous and 24-hour periods.   Check SET WINDOW.
# Diffuse+direct Rs
# Eq 19
# Eq 18,14
# Eq 13
# Terrain_albedo
# Default 0.2
# Eq 18, 14
# Eq 19
# ____________________________________________________________________________________________________________________________
# ________________________________________________________________________________________________________________________________
# INPUT
# CALCU-
# LATIONS
# OUTPUT
#
# set cell size for the model
#
# SET CELLSIZE MIN;
#
# set window for the model
#
# SET WINDOW UNION;
#
# set area of interest for the model
#
# SET AOI NONE;
#
# declarations
#
#Integer RASTER n5_slope_degree_p33r37_pr_nad83_cc FILE OLD PUBINPUT NEAREST NEIGHBOR AOI NONE "f:/metric_etrm_jornada_p33r37_2000_2011/commons/slope_degree_p33r37_pr_nad83_cc.img";
#Integer RASTER n6_aspect_degree_p33r37_pr_nad83_cc FILE OLD PUBINPUT NEAREST NEIGHBOR AOI NONE "f:/metric_etrm_jornada_p33r37_2000_2011/commons/aspect_degree_p33r37_pr_nad83_cc.img";
#Float RASTER n15_temp;
#Float RASTER n18_temp;
#Float RASTER n19_temp;
#Float RASTER n33_cos_theta_20110921_p33r37_l5 FILE DELETE_IF_EXISTING PUBOUT USEALL ATHEMATIC FLOAT SINGLE "f:/metric_etrm_jornada_p33r37_2000_2011/example_run_model_001_gmd/cos_theta_20110921_p33r37_l5.img";
#Float RASTER n56_dem_p33r37_pr_nad83_cc FILE OLD PUBINPUT NEAREST NEIGHBOR AOI NONE "f:/metric_etrm_jornada_p33r37_2000_2011/commons/dem_p33r37_pr_nad83_cc.img";
#Float RASTER n58_temp;
#Float RASTER n69_rso_20110921_p33r37_l5 FILE DELETE_IF_EXISTING PUBOUT USEALL ATHEMATIC FLOAT SINGLE "f:/metric_etrm_jornada_p33r37_2000_2011/example_run_model_001_gmd/rso_20110921_p33r37_l5.img";
#Float RASTER n71_temp;
#Float RASTER n74_rso24_20110921_p33r37_l5 FILE DELETE_IF_EXISTING PUBOUT USEALL ATHEMATIC FLOAT SINGLE "f:/metric_etrm_jornada_p33r37_2000_2011/example_run_model_001_gmd/rso24_20110921_p33r37_l5.img";
#Float RASTER n77_rso_flat_20110921_p33r37_l5 FILE DELETE_IF_EXISTING PUBOUT USEALL ATHEMATIC FLOAT SINGLE "f:/metric_etrm_jornada_p33r37_2000_2011/example_run_model_001_gmd/rso_flat_20110921_p33r37_l5.img";
#Float RASTER n80_temp;
#Float RASTER n83_rad_adj_coeff_20110921_p33r37_l5 FILE DELETE_IF_EXISTING PUBOUT USEALL ATHEMATIC FLOAT SINGLE "f:/metric_etrm_jornada_p33r37_2000_2011/example_run_model_001_gmd/rad_adj_coeff_20110921_p33r37_l5.img";
#Float RASTER n88_temp;
#Float RASTER n94_temp;
#Float RASTER n95_temp;
#Float RASTER n98_trans_bb_20110921_p33r37_l5 FILE DELETE_IF_EXISTING PUBOUT USEALL ATHEMATIC FLOAT SINGLE "f:/metric_etrm_jornada_p33r37_2000_2011/example_run_model_001_gmd/trans_bb_20110921_p33r37_l5.img";
#FLOAT TABLE n21_sunhourangles [48];
#FLOAT SCALAR n1_Float;
#FLOAT SCALAR n2_Float;
#FLOAT SCALAR n9_Float;
#FLOAT SCALAR n31_Float;
#FLOAT SCALAR n60_Float;
#FLOAT SCALAR n62_Float;
#FLOAT SCALAR n109_Float;
#
# load scalar n1_Float
#
# n1_Float = 264;
# #
# # load scalar n2_Float
# #
# n2_Float = 33.2;
# #
# # load scalar n9_Float
# #
# n9_Float = 0.1;
# #
# # load table n21_sunhourangles
# #
# n21_sunhourangles = TABLE(-3.076142807, -2.945243113, -2.814343419, -2.683443725, -2.552544031, -2.421644337, -2.290744643, -2.159844949, -2.028945255, -1.898045562, -1.767145868, -1.636246174, -1.50534648, -1.374446786, -1.243547092, -1.112647398, -0.981747704, -0.85084801, -0.719948316, -0.5890486230000001, -0.458148929, -0.327249235, -0.196349541, -0.06544984700000001, 0.06544984700000001, 0.196349541, 0.327249235, 0.458148929, 0.5890486230000001, 0.719948316, 0.85084801, 0.981747704, 1.112647398, 1.243547092, 1.374446786, 1.50534648, 1.636246174, 1.767145868, 1.898045562, 2.028945255, 2.159844949, 2.290744643, 2.421644337, 2.552544031, 2.683443725, 2.814343419, 2.945243113, 3.076142807);
# #
# # load scalar n31_Float
# #
# n31_Float = 10.4;
# #
# # load scalar n60_Float
# #
# n60_Float = 1;
# #
# # load scalar n62_Float
# #
# n62_Float = 1;
# #
# # load scalar n109_Float
# #
# n109_Float = 0.2;
# #
# # function definitions
# #
# print (n109_float)
# n88_temp = 101.3*((293.0-0.0065*n56_dem_p33r37_pr_nad83_cc)/293.0)**5.26
#
# n71_temp = 0.98*exp(-0.00146*$n88_temp/$n62_Float/(sin(0.85+0.3*$n2_Float*pi/180.0*sin(2*pi/365.0*$n1_Float-1.39)-0.42*($n2_Float*pi/180.0)**2))-0.075*((0.14*$n60_Float*$n88_temp+2.1)/(sin(0.85+0.3*$n2_Float*pi/180.0*sin(2*pi/365.0*$n1_Float-1.39)-0.42*($n2_Float*pi/180.0)**2)))**0.4)
# ;
# n95_temp = EITHER (0.35-0.36*$n71_temp) IF ( $n71_temp >= 0.15 ) OR (0.18+0.82*$n71_temp) OTHERWISE
# ;
# n58_temp = 0.98*exp(-0.00146*$n88_temp/$n62_Float/(sin($n2_Float*pi/180.0)*sin(0.409*sin(2*pi/365*$n1_Float-1.39))+cos($n2_Float*pi/180.0)*cos(0.409*sin(2*pi/365*$n1_Float-1.39))*cos(pi/12*($n31_Float-12)))-0.075*((0.14*$n60_Float*$n88_temp+2.1)/(sin($n2_Float*pi/180.0)*sin(0.409*sin(2*pi/365*$n1_Float-1.39))+cos($n2_Float*pi/180.0)*cos(0.409*sin(2*pi/365*$n1_Float-1.39))*cos(pi/12*($n31_Float-12))))**0.4)
# ;
# n94_temp = EITHER (0.35-0.36*$n58_temp) IF ( $n58_temp >= 0.15 ) OR (0.18+0.82*$n58_temp) OTHERWISE
# ;
# n98_trans_bb_20110921_p33r37_l5 = $n58_temp+$n94_temp
# ;
# n80_temp = ($n71_temp+$n95_temp) * 1367.0/PI*(1+0.033*COS(2*PI/365.0*$n1_Float))*(ACOS(-TAN($n2_Float*PI/180.0)*TAN(0.409*SIN(2*PI/365.0*$n1_Float-1.39)))*SIN(0.409*SIN(2*PI/365.0*$n1_Float-1.39))*SIN($n2_Float/180.0*PI)+COS(0.409*SIN(2*PI/365.0*$n1_Float-1.39))*COS($n2_Float/180.0*PI)*SIN(ACOS(-TAN($n2_Float*PI/180.0)*TAN(0.409*SIN(2*PI/365.0*$n1_Float-1.39)))))
# ;
# n77_rso_flat_20110921_p33r37_l5 = ($n58_temp+$n94_temp)*1367.0*(1+0.033*cos($n1_Float*2.0*pi/365.0))*cos(pi/2.0-asin(sin($n2_Float*pi/180)*sin(0.409*sin(2*pi/365.0*$n1_Float-1.39))+cos($n2_Float*pi/180)*cos(0.409*sin(2*pi/365.0*$n1_Float-1.39))*cos(pi/12.0*($n31_Float-12.0))))
# ;
# #define n22_memory Binary($n5_slope_degree_p33r37_pr_nad83_cc*0+1
# \
# )
# #define n25_memory Binary(STACKLAYERS ( $n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory,$n22_memory )
# \
# )
# n19_temp = cos(0.409*sin(2*pi/365.0*$n1_Float-1.39))*sin($n5_slope_degree_p33r37_pr_nad83_cc*pi/180.0)*sin(($n6_aspect_degree_p33r37_pr_nad83_cc-180)*PI/180)
# ;
# n18_temp = cos(0.409*sin(2*pi/365.0*$n1_Float-1.39))*cos($n2_Float*pi/180.0)*cos($n5_slope_degree_p33r37_pr_nad83_cc*pi/180.0)+cos(0.409*sin(2*pi/365.0*$n1_Float-1.39))*sin($n2_Float*pi/180.0)*sin($n5_slope_degree_p33r37_pr_nad83_cc*pi/180.0)*cos(($n6_aspect_degree_p33r37_pr_nad83_cc-180)*PI/180)
# ;
# n15_temp = sin(0.409*sin(2*pi/365.0*$n1_Float-1.39))*sin($n2_Float*pi/180.0)*cos($n5_slope_degree_p33r37_pr_nad83_cc*pi/180.0)-sin(0.409*sin(2*pi/365.0*$n1_Float-1.39))*cos($n2_Float*pi/180.0)*sin($n5_slope_degree_p33r37_pr_nad83_cc*pi/180.0)*cos(($n6_aspect_degree_p33r37_pr_nad83_cc-180)*PI/180)
# ;
# n33_cos_theta_20110921_p33r37_l5 = CONDITIONAL {
#   ( $n15_temp+$n18_temp*cos(pi/12.0*($n31_Float-12.0))+$n19_temp*sin(pi/12.0*($n31_Float-12.0)) <= $n9_Float) ($n9_Float /cos($n5_slope_degree_p33r37_pr_nad83_cc*pi/180.0)) ,
#   ($n31_Float >= 0) (($n15_temp+$n18_temp*cos(pi/12.0*($n31_Float-12))+$n19_temp*sin(pi/12.0*($n31_Float-12)))/cos($n5_slope_degree_p33r37_pr_nad83_cc*pi/180.0))
# }
# ;
# n69_rso_20110921_p33r37_l5 = /*Rs0 on sloped surfaces see Allen et al 2006 Agric For Meteorol. All angles in rad*/
# /*below is the beam Rs0 component eq 12-8*/
# $n58_temp*1367.0*(1+0.033*cos($n1_Float*2.0*pi/365.0))*$n33_cos_theta_20110921_p33r37_l5 +
#
# /*below is the diffuse Rs0 componet, zenith angle is cos(pi/2-eq9 in manual); Allen et al 2006 eq 31, 32*/
# $n94_temp*1367.0*(1+0.033*cos($n1_Float*2.0*pi/365.0))*cos(pi/2.0-asin(sin($n2_Float*pi/180)*sin(0.409*sin(2*pi/365.0*$n1_Float-1.39))+cos($n2_Float*pi/180)*cos(0.409*sin(2*pi/365.0*$n1_Float-1.39))*cos(pi/12.0*($n31_Float-12.0)))) * (0.75 + 0.25*cos($n5_slope_degree_p33r37_pr_nad83_cc*pi/180) - 0.5*($n5_slope_degree_p33r37_pr_nad83_cc*pi/180)/pi) +
#
# /*below is the reflected Rs0 component, Allen et al eq 32, 36 */
# $n77_rso_flat_20110921_p33r37_l5 * $n109_Float * (1-(0.75 + 0.25*cos($n5_slope_degree_p33r37_pr_nad83_cc*pi/180) - 0.5*($n5_slope_degree_p33r37_pr_nad83_cc*pi/180)/pi))
#
#
#
#
# ;
# #define n28_memory Float((STACK SUM (
# \
# CONDITIONAL { ((sin(0.409*sin(2*pi/365.0*$n1_Float-1.39))*sin($n2_Float*pi/180.0) + cos(0.409*sin(2*pi/365.0*$n1_Float-1.39))*cos($n2_Float*pi/180.0)*cos($n21_sunhourangles)) <= 0) (0) , ((($n25_memory*$n15_temp+$n25_memory*$n18_temp*cos($n21_sunhourangles)+$n25_memory*$n19_temp*sin($n21_sunhourangles))) <= 0) (0) , ($n1_Float >= 0) (($n25_memory*$n15_temp+$n25_memory*$n18_temp*cos($n21_sunhourangles)+$n25_memory*$n19_temp*sin($n21_sunhourangles))) }
# \
# ) )
# \
# *0.5*(1+0.033*cos($n1_Float*2*pi/365))*1367.0/24.0
# \
# )
# #define n39_memory Float(EITHER ($n28_memory/cos($n5_slope_degree_p33r37_pr_nad83_cc*pi/180.0))
# \
#   IF ( $n28_memory > $n9_Float*1367.0/pi*(1+0.033*cos($n1_Float*2*pi/365.0))*(acos(-tan($n2_Float*pi/180.0)*tan(0.409*sin(2*pi/365.0*$n1_Float-1.39)))*sin($n2_Float*pi/180.0)*sin(0.409*sin(2*pi/365.0*$n1_Float-1.39))+cos($n2_Float*pi/180.0)*cos(0.409*sin(2*pi/365.0*$n1_Float-1.39))*sin(acos(-tan($n2_Float*pi/180.0)*tan(0.409*sin(2*pi/365.0*$n1_Float-1.39))))) )
# \
#   OR ($n9_Float*1367.0/pi*(1+0.033*cos($n1_Float*2*pi/365.0))*(acos(-tan($n2_Float*pi/180.0)*tan(0.409*sin(2*pi/365.0*$n1_Float-1.39)))*sin($n2_Float*pi/180.0)*sin(0.409*sin(2*pi/365.0*$n1_Float-1.39))+cos($n2_Float*pi/180.0)*cos(0.409*sin(2*pi/365.0*$n1_Float-1.39))*sin(acos(-tan($n2_Float*pi/180.0)*tan(0.409*sin(2*pi/365.0*$n1_Float-1.39))))) ) /cos($n5_slope_degree_p33r37_pr_nad83_cc*pi/180.0)
# \
#   OTHERWISE
# \
#
# \
# )
# n74_rso24_20110921_p33r37_l5 = /*Rs0 on sloped surfaces see Allen et al 2006 Agric For Meteorol. All angles in rad*/
# /*below is the beam Rs0 component*/
# $n71_temp*$n39_memory +
#
# /*below is the diffuse Rs0 component*/
# $n95_temp*1367.0/pi*(1+0.033*cos(2*pi/365.0*$n1_Float))*(acos(-tan($n2_Float*pi/180.0)*tan(0.409*sin(2*pi/365.0*$n1_Float-1.39)))*sin(0.409*sin(2*pi/365.0*$n1_Float-1.39))*sin($n2_Float/180.0*pi)+cos(0.409*sin(2*pi/365.0*$n1_Float-1.39))*cos($n2_Float/180.0*pi)*sin(acos(-tan($n2_Float*pi/180.0)*tan(0.409*sin(2*pi/365.0*$n1_Float-1.39))))) * (0.75 + 0.25*cos($n5_slope_degree_p33r37_pr_nad83_cc*pi/180) - 0.5*($n5_slope_degree_p33r37_pr_nad83_cc*pi/180)/pi) +
#
# /*below is the reflected Rs0 component*/
# $n80_temp * $n109_Float * (1-(0.75 + 0.25*cos($n5_slope_degree_p33r37_pr_nad83_cc*pi/180) - 0.5*($n5_slope_degree_p33r37_pr_nad83_cc*pi/180)/pi))
#
# ;
# n83_rad_adj_coeff_20110921_p33r37_l5 = EITHER (1.5) IF ( ($n77_rso_flat_20110921_p33r37_l5/$n69_rso_20110921_p33r37_l5*$n74_rso24_20110921_p33r37_l5/$n80_temp) > 1.5 ) OR ( $n77_rso_flat_20110921_p33r37_l5/$n69_rso_20110921_p33r37_l5*$n74_rso24_20110921_p33r37_l5/$n80_temp ) OTHERWISE
#