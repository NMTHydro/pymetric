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
from matplotlib import pyplot as plt
# ============= standard library imports ========================
from code.terrain_functions.raster_util import convert_raster_to_array, get_raster_geo_attributes, convert_array_to_raster

# TODO - MAKE SURE THAT YOU USE PYMETRIC-WRITTEN FUNCTIONS TO READ THE RASTERS AS ARRAYS DONT REPEAT FUNCTIONS ALREADY WRITTEN BY DRI


def integral_cos_theta(lower_limit, upper_limit, sin_decl, cos_decl,sin_lat, cos_lat, sin_slope, cos_slope,
                       sin_aspect, cos_aspect, raster_mode=False):
    """"""
    #   Equation [5] in Allen (2006)
    #   lower_limit = lower_limit of integral as hourangle in radisns
    #   upper_limit = upper limit of integral as hourangle in radians

    int_cos_theta = sin_decl * sin_lat * cos_slope * (upper_limit - lower_limit) \
        - sin_decl * cos_lat * sin_slope * cos_aspect * (upper_limit - lower_limit) \
        + cos_decl * cos_lat * cos_slope * (np.sin(upper_limit) - np.sin(lower_limit)) \
        + cos_decl * sin_lat * sin_slope * cos_aspect * (np.sin(upper_limit) - np.sin(lower_limit)) \
        - cos_decl * sin_slope * sin_aspect * (np.cos(upper_limit) - np.cos(lower_limit))

    if not raster_mode:
        print('\n','lower_limit = {:.5f}'.format(lower_limit))
        print('\n','upper_limit = {:.5f}'.format(upper_limit))
        print('\n','upper limit - lower_limit = {:.5f}'.format(upper_limit - lower_limit))

    return int_cos_theta


def calc_two_daytime_integration_limits(omega_rise_pixel_24, omega_set_pixel_24, lower_int_limit_rise,
                                        upper_int_limit_set, sin_slope, sin_lat, sin_decl, sin_aspect, cos_slope,
                                        cos_lat, cos_decl, cos_aspect, a, b, c, quadratic_function, raster_mode=False):

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

    if raster_mode:
        omega_rise_pixel_24[omega_set_pixel_24 < omega_rise_pixel_24] = omega_set_pixel_24[omega_set_pixel_24 < omega_rise_pixel_24]

    else:

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

    # This means that there is no direct solar beam during the day; the slope is always shaded.

    # STEP D - Section (iii)
    # The sine values in Eqs. [13a and 13b] cannot be smaller than -1 or larger than +1. This issue has been taken care of
    #   under STEP B and STEP C.

    # STEP D - Section (iv)
    #    Check if possibility exists for two periods of direct beam radiation during the day.
    #    The two daytime integration limits are defined as follows:
    #       omega_set_during_day_pixel_24 = the time angle when the center of the solar disk disappears the first time
    #       omega_rise_during_day_pixel_24 = the time angle when the center of the solar disk reappears over the surface

    if raster_mode:

        # instantiate the new array
        omega_set_during_day_pixel_24 = np.empty(omega_rise_pixel_24.shape)
        omega_rise_during_day_pixel_24 = np.empty(omega_rise_pixel_24.shape)

        # if sin_slope is less than, there is only one period of direct beam radiation during the day.
        omega_set_during_day_pixel_24[sin_slope <= (sin_lat * cos_decl + cos_lat * sin_decl)] = 0.000000
        omega_rise_during_day_pixel_24[sin_slope <= (sin_lat * cos_decl + cos_lat * sin_decl)] = 0.000000

    else:
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

    if raster_mode:

        sin_A[sin_A < -1.0] = -1.0
        sin_A[sin_A > 1.0] = 1.0
        sin_B[sin_B < -1.0] = -1.0
        sin_B[sin_B > 1.0] = 1.0

        A = np.arcsin(sin_A)
        B = np.arcsin(sin_B)

        #     STEP D - Section(iv - d)
        omega_set_during_day_pixel_24 = np.minimum(A, B)
        omega_rise_during_day_pixel_24 = np.maximum(A, B)

    else:

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


    #     STEP D - Section(iv - e)

    cos_theta_omega_set_during_day_pixel_24 = - a + b * np.cos(omega_set_during_day_pixel_24) + c \
                                              * np.sin(omega_set_during_day_pixel_24)
    cos_theta_omega_rise_during_day_pixel_24 = - a + b * np.cos(omega_rise_during_day_pixel_24) + c \
                                               * np.sin(omega_rise_during_day_pixel_24)

    if raster_mode:
        omega_set_during_day_pixel_24[(cos_theta_omega_set_during_day_pixel_24 < -0.001) |
                                      (cos_theta_omega_set_during_day_pixel_24 > 0.001)] = \
            - np.pi - omega_set_during_day_pixel_24[(cos_theta_omega_set_during_day_pixel_24 < -0.001) |
                                                    (cos_theta_omega_set_during_day_pixel_24 > 0.001)]

        omega_rise_during_day_pixel_24[(cos_theta_omega_rise_during_day_pixel_24 < -0.001) |
                                       (cos_theta_omega_rise_during_day_pixel_24 > 0.001)] = \
            np.pi - omega_rise_during_day_pixel_24[(cos_theta_omega_rise_during_day_pixel_24 < -0.001) |
                                                   (cos_theta_omega_rise_during_day_pixel_24 > 0.001)]

        #     STEP D - Section(iv - f and g) RASTER MODE

        omega_set_during_day_pixel_24[omega_set_during_day_pixel_24 < omega_rise_pixel_24] = \
            omega_rise_pixel_24[omega_set_during_day_pixel_24 < omega_rise_pixel_24]

        omega_rise_during_day_pixel_24[omega_rise_during_day_pixel_24 > omega_set_pixel_24] = \
            omega_set_pixel_24[omega_rise_during_day_pixel_24 > omega_set_pixel_24]

    else:

        print('\n', 'cos_theta_omega_set_during_day_pixel_24 = {:.5f}'.format(cos_theta_omega_set_during_day_pixel_24),
                  '      ',
                  'cos_theta_omega_rise_during_day_pixel_24 = {:.5f}'.format(cos_theta_omega_rise_during_day_pixel_24))

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

        #     STEP D - Section(iv - f and g)

        if omega_set_during_day_pixel_24 < omega_rise_pixel_24:
                omega_set_during_day_pixel_24 = omega_rise_pixel_24
        if omega_rise_during_day_pixel_24 > omega_set_pixel_24:
                omega_rise_during_day_pixel_24 = omega_set_pixel_24



        print('\n', 'omega_set_during_day_pixel_24 = {:.5f}'.format(omega_set_during_day_pixel_24), '      ',
                  'omega_rise_during_day_pixel_24 = {:.5f}'.format(omega_rise_during_day_pixel_24))
        print('\n', 'omega_set_pixel_24 = {:.5f}'.format(omega_set_pixel_24), '      ',
                  'omega_rise_pixel_24 = {:.5f}'.format(omega_rise_pixel_24))


    X = sin_decl * sin_lat * cos_slope * (omega_rise_during_day_pixel_24 - omega_set_during_day_pixel_24) \
        - sin_decl * cos_lat * sin_slope * cos_aspect * (omega_rise_during_day_pixel_24 - omega_set_during_day_pixel_24) \
        + cos_decl * cos_lat * cos_slope * (np.sin(omega_rise_during_day_pixel_24)
                                            - np.sin(omega_set_during_day_pixel_24)) + cos_decl * sin_lat * sin_slope \
        * cos_aspect * (np.sin(omega_rise_during_day_pixel_24) - np.sin(omega_set_during_day_pixel_24)) \
        - cos_decl * sin_slope * sin_aspect * (np.cos(omega_rise_during_day_pixel_24)
                                               - np.cos(omega_set_during_day_pixel_24))

    if raster_mode:
        # Just doing the if/else so that Jan's print statements will still work. No need to modify anything for raster.
        pass
    else:

        if X < 0:
            print('\n','x = {:.5f}'.format(X),'there are two periods of direct beam radiation during the day')
        else:
            print('\n','x = {:.5f}'.format(X),'there is only one period of direct beam radiation during the day')

    return lower_int_limit_rise, upper_int_limit_set, omega_rise_pixel_24, omega_set_pixel_24, \
           omega_set_during_day_pixel_24, omega_rise_during_day_pixel_24, X


# This function applies when Eq [7] from Allen 2006 is met.
def calc_two_integration_limits(LatRad, DeclRad, a, b, c, raster_mode=False, raster_shape=None):
    # This function is based on Allen (2006)
    # The function is for the northern hemisphere below the arctic circle, i.e. there are no 24 hour days or nights
    # See Appendix A in Allen (2006)

    # ========== STEP A ==========
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

    # ========== STEP D - Section (i) ==========
    #     In Eqs. [13a, b] the expression (b^2 + C^2 - a^2) under the square root sign must be limited to  > 0 for
    #         numerical stability. Therefore, if (b^2 + C^2 - a^2) is 0 or less, it is set equal to 0.001.

    quadratic_function = b**2 + c**2 - a**2

    if raster_mode:
        print('\n', 'STEP D - Section (i)', '\n', 'Before check on positive value', '      ',
              'quadratic_function = \n {}'.format(quadratic_function))
        quadratic_function[quadratic_function <= 0.0] = 0.001
        print(' After check on positive value \n '
              'quadratic_function =  {}'.format(quadratic_function[quadratic_function <= 0.0].any()))
    else:
        print('\n', 'STEP D - Section (i)', '\n', 'Before check on positive value', '      ', 'quadratic_function =  {:.6f}'.format(quadratic_function))
        if quadratic_function <= 0.0:
            quadratic_function = 0.001

        print(' After check on positive value', '       ', 'quadratic_function =  {:.6f}'.format(quadratic_function))

    # ========== STEP B ==================
    print('\n', 'STEP B: Determine the beginning integration limit at sun rise',
          '\n','\n', 'STEP B - Section (i)')
    #
    # STEP B: Determine the beginning integration limit at sun rise
    # STEP B - Section (i)
    # Calculate the sine of the sunrise angle on a specific pixel using Eq. [13a]
    sin_omega_rise_pixel = (a * c - b * np.sqrt(quadratic_function)) / (b**2 + c**2)

    if raster_mode:
        print(' Before check on sin values within ±1    sin_omega_rise_pixel = ')

        sin_omega_rise_pixel[sin_omega_rise_pixel < -1.0 ] = -1.0
        sin_omega_rise_pixel[sin_omega_rise_pixel > 1.0] = 1.0
        print(' After check on sin values within ±1 \n '
              'sin_omega_rise_pixel'
              ' = {}'.format(sin_omega_rise_pixel[(sin_omega_rise_pixel < -1.0) & (sin_omega_rise_pixel > 1.0)].any()))

    else:
        print(' Before check on sin values within ±1    sin_omega_rise_pixel = {:.5f}'.format(sin_omega_rise_pixel))

        if sin_omega_rise_pixel < -1.0:
            sin_omega_rise_pixel = -1.0
        if sin_omega_rise_pixel > 1.0:
            sin_omega_rise_pixel = 1.0

        print(' After check on sin values within ±1     sin_omega_rise_pixel = {:.5f}'.format(sin_omega_rise_pixel))


    print('\n', 'STEP B - Section (ii)')

    # STEP B - Section (ii)
    omega_rise_pixel = np.arcsin(sin_omega_rise_pixel)

    if raster_mode:
        print(' omega_rise_pixel = \n {}'.format(omega_rise_pixel))
    else:
        print(' omega_rise_pixel = {:.5f}'.format(omega_rise_pixel))


    print('\n', 'STEP B - Section (iii)')

    # STEP B - Section (iii)
    # Calculate cosine of theta using omega_rise_pixel using Eq. [14]
    #cos_theta_omega_rise_pixel_eq14 = - a + b * np.cos(omega_rise_pixel) + c * np.sin(omega_rise_pixel)
    cos_theta_omega_rise_pixel = - a + b * np.cos(omega_rise_pixel) + c * np.sin(omega_rise_pixel)



    if raster_mode:
        # STEP B - Section (iv) RASTER MODE
        omega_rise_pixel_24 = np.empty(cos_theta_omega_rise_pixel.shape)
        # todo - come back to here to make sure that we get the right candidate is selected

        print('omega rise pixel 24', omega_rise_pixel_24.shape)
        print('cos_theta_omega_rise_hor', cos_theta_omega_rise_hor.shape)
        print('cos_theta_omega_rise_pixel', cos_theta_omega_rise_hor.shape)
        print('omega_rise_pixel', omega_rise_pixel.shape)


        # TODO - THIS REALLY WORKS YEAH?!?!?
        omega_rise_pixel_24[(cos_theta_omega_rise_hor <= cos_theta_omega_rise_pixel) & (cos_theta_omega_rise_hor < 0.001)] = omega_rise_pixel[(cos_theta_omega_rise_hor <= cos_theta_omega_rise_pixel) & (cos_theta_omega_rise_hor < 0.001)]
        print('the omega_rise_pixel array \n', omega_rise_pixel_24)

        print('For the values that have not been set to omega_rise_pixel we need to select new candidates')

        # STEP B - Section (v) RASTER MODE

        omega_rise_pixel_x = -np.pi - omega_rise_pixel

        # STEP B - Section (v-a) RASTER MODE
        cos_theta_omega_rise_pixel_x = - a + b * np.cos(omega_rise_pixel_x) + c * np.sin(omega_rise_pixel_x)
        omega_rise_pixel_24[cos_theta_omega_rise_pixel_x > 0.001] = omega_rise_hor[cos_theta_omega_rise_pixel_x > 0.001]

        # Section B_v_b RASTER MODE

        omega_rise_pixel_24[(cos_theta_omega_rise_pixel_x <= 0.001) & (omega_rise_pixel_x <= omega_rise_hor)] = omega_rise_hor[(cos_theta_omega_rise_pixel_x <= 0.001) & (omega_rise_pixel_x <= omega_rise_hor)]
        omega_rise_pixel_24[(cos_theta_omega_rise_pixel_x <= 0.001) & (omega_rise_pixel_x > omega_rise_hor)] = - np.pi - omega_rise_pixel[(cos_theta_omega_rise_pixel_x <= 0.001) & (omega_rise_pixel_x > omega_rise_hor)]

        # STEP B - Section (vi) RASTER MODE
        omega_rise_pixel_24[omega_rise_pixel_24 < omega_rise_hor] = omega_rise_hor[omega_rise_pixel_24 < omega_rise_hor]
        # STEP B - Section (vii) RASTER MODE
        lower_int_limit_rise = omega_rise_pixel_24

    else:
        # STEP B - Section (iv)
        print(' Before check on cos_theta_omega_rise_pixel in STEP B - Section (iv)     '
              'cos_theta_omega_rise_pixel = {:.5f}'.format(cos_theta_omega_rise_pixel))


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

    if raster_mode:

        print('\n', '\n', 'STEP C: Determine the ending integration limit at sun set', '\n', '\n',
              'STEP C - Section (i)', '\n',
              'Before check on sin values within ±1 sin_omega_rise_pixel ')

        sin_omega_set_pixel[sin_omega_set_pixel < -1.0] = -1.0
        sin_omega_set_pixel[sin_omega_set_pixel > 1.0] = 1.0

    else:

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

    if not raster_mode:
        print(' omega_set_pixel = {:.5f}'.format(omega_set_pixel))

    print('\n', 'STEP C - Section (iii)')
    # STEP C - Section (iii)
    # Calculate cosine of theta using omega_set_pixel using Eq. [14]
    cos_theta_omega_set_pixel = - a + b * np.cos(omega_set_pixel) + c * np.sin(omega_set_pixel)

    if raster_mode:

        # STEP C - Section (iv)
        print('STEP C - Section (iv) RASTER MODE')
        omega_set_pixel_24 = np.empty(cos_theta_omega_set_pixel.shape)
        omega_set_pixel_24[(cos_theta_omega_set_hor <= cos_theta_omega_set_pixel) & (cos_theta_omega_set_pixel < 0.001)] = omega_set_pixel[(cos_theta_omega_set_hor <= cos_theta_omega_set_pixel) & (cos_theta_omega_set_pixel < 0.001)]

        # STEP C - Section (v)
        print('STEP C - Section (v)')
        omega_set_pixel_x = np.pi - omega_set_pixel

        # STEP C - Section (v-a)
        print('STEP C - Section (v-a)')
        cos_theta_omega_set_pixel_x = - a + b * np.cos(omega_set_pixel_x) + c * np.sin(omega_set_pixel_x)
        omega_set_pixel_24[cos_theta_omega_set_pixel_x > 0.001] = omega_set_hor[cos_theta_omega_set_pixel_x > 0.001]

        # Section C_v_b
        print('Section C_v_b')
        omega_set_pixel_24[(cos_theta_omega_set_pixel_x <= 0.001) & (omega_set_pixel_x >= omega_set_hor)] = omega_set_hor[(cos_theta_omega_set_pixel_x <= 0.001) & (omega_set_pixel_x >= omega_set_hor)]

        # Section C_v_c
        print('STEP C - Section (vi)')
        omega_set_pixel_24[(cos_theta_omega_set_pixel_x <= 0.001) & (omega_set_pixel_x < omega_set_hor)] = np.pi - omega_set_pixel[(cos_theta_omega_set_pixel_x <= 0.001) & (omega_set_pixel_x < omega_set_hor)]

        # STEP C - Section (vi)
        print('Section C_v_c')
        omega_set_pixel_24[omega_set_pixel_24 > omega_set_hor] = omega_set_hor[omega_set_pixel_24 > omega_set_hor]

        # STEP C - Section (vii)
        print('STEP C - Section (vii)')
        upper_int_limit_set = omega_set_pixel_24

    else:
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
        #     Section C_v_b
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


def calc_solar_topo_image_variables(DOY, LatDeg, SlopeDeg, AspectDeg, local_time, Lz, Lm, raster_mode=False):
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

    if not raster_mode:
        if cos_theta_adj < 0.1:
            cos_theta_adj = 0.1
        if cos_theta_adj > 10.0:
            cos_theta_adj = 10.0
    if raster_mode:
        # reset values in the array in-place
        cos_theta_adj[cos_theta_adj < 0.1] = 0.1
        cos_theta_adj[cos_theta_adj > 10.0] = 10.0

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
    # this line doesn't change if raster or not but i rewrite it for clarity
    solar_azimuth_angle_landsat = 180 + solar_azimuth_angle

    # RETURN all parameters
    return LatRad, DeclRad, DeclDeg, SlopeRad, AspectRad, cos_lat, sin_lat, cos_slope, sin_slope, \
    cos_aspect, sin_aspect, sin_decl, cos_decl, bSc, Sc, HourAngleRad, cos_hourangle, sin_hourangle, a, b, c, cos_theta_unadj, \
    cos_theta_adj, cos_theta_horizontal_pixel, solar_incidence_angle_theta, \
           solar_elevation_angle, solar_azimuth_angle, solar_azimuth_angle_landsat, solar_time


def calculate_integration_limits(DOY, LatDeg, SlopeDeg, AspectDeg, local_time, Lz, Lm, output_dir=None, raster_mode=False):
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
        # get the geo attributes of a representative file:
        geo = get_raster_geo_attributes(LatDeg)
        print('geo', geo)

        # read in each raster parameter as a numpy array and continue
        LatDeg = convert_raster_to_array(LatDeg)
        SlopeDeg = convert_raster_to_array(SlopeDeg)
        AspectDeg = convert_raster_to_array(AspectDeg)

        # GET the shape of the raster
        raster_shape = LatDeg.shape

        # flatten every array
        LatDeg = LatDeg.flatten()
        SlopeDeg = SlopeDeg.flatten()
        AspectDeg = AspectDeg.flatten()


    LatRad, DeclRad, DeclDeg, SlopeRad, AspectRad, cos_lat, sin_lat, cos_slope, sin_slope, cos_aspect, sin_aspect, \
    sin_decl, cos_decl, bSc, Sc, HourAngleRad, cos_hourangle, sin_hourangle, a, b, c, cos_theta_unadj, cos_theta_adj, \
    cos_theta_horizontal_pixel, solar_incidence_angle_theta, solar_elevation_angle, solar_azimuth_angle, \
    solar_azimuth_angle_landsat, solar_time = \
        calc_solar_topo_image_variables(DOY, LatDeg, SlopeDeg, AspectDeg, local_time, Lz, Lm, raster_mode)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Testing Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if raster_mode:
        solar_topo_raster_vars = [LatRad,  SlopeRad, AspectRad, cos_lat, sin_lat, cos_slope, sin_slope, cos_aspect, sin_aspect, \
      a, b, c, cos_theta_unadj, cos_theta_adj, \
    cos_theta_horizontal_pixel, solar_incidence_angle_theta, solar_elevation_angle, solar_azimuth_angle, \
    solar_azimuth_angle_landsat]

        solar_topo_raster_vars_names = ['LatRad',  'SlopeRad', 'AspectRad', 'cos_lat', 'sin_lat', 'cos_slope', 'sin_slope', 'cos_aspect', 'sin_aspect', \
    'a', 'b', 'c', 'cos_theta_unadj', 'cos_theta_adj', \
    'cos_theta_horizontal_pixel', 'solar_incidence_angle_theta', 'solar_elevation_angle', 'solar_azimuth_angle', \
    'solar_azimuth_angle_landsat']

        print('NON RASTER parameters\n')
        print('\n',
              'Output Function calc_solar_topo_image_variables(DOY, LatDeg, SlopeDeg, AspectDeg, local_time, Lz, Lm)')
        print('\n', 'DOY = {:.0f}'.format(DOY),
              '   ', 'local_time = {:.2f}'.format(local_time), '  ', 'solar_time = {:.2f}'.format(solar_time), '  '
              , 'Lz = {:.0f}'.format(Lz), '  ', 'Lm = {:.0f}\n'.format(Lm))

        for name, var in zip(solar_topo_raster_vars_names, solar_topo_raster_vars):
            print('{} is a raster\n'.format(name))
            var = np.reshape(var, raster_shape)
            print(var, '\n')

            # plt.imshow(i)
            # plt.show()

        print(
            ' \n more non-raster params \n DeclRad -> {}, DeclDeg -> {} sin_decl -> {}, cos_decl -> {}, bSc -> {}  Sc -> {} \n  HourAngleRad {}, '
            'cos_hourangle {} sin_hourangle {}'.format(DeclRad, DeclDeg, sin_decl, cos_decl, bSc, Sc,
                                                       HourAngleRad, cos_hourangle, sin_hourangle))

        print('\n',
              'slope azimuth METRIC: 0 degree for slopes oriented due south, ±180 degrees due north, -90 degrees due east,'
              ' +90 degrees due west', '\n',
              'solar azimuth Landsat: 0 degree for sun in north, 90 degrees sun in east, '
              '180 degrees sun in south, 270 degrees sun in west')

        print('ALLEN 2010 Slope and cos theta \n')

    else:
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
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # RUN calc_two_integration_limits function
    if raster_mode:

        omega_set_hor, omega_rise_hor, cos_theta_omega_set_hor, cos_theta_omega_rise_hor, sin_omega_rise_pixel, \
        omega_rise_pixel_24, sin_omega_set_pixel, \
        omega_set_pixel_24, lower_int_limit_rise, upper_int_limit_set, quadratic_function = calc_two_integration_limits(LatRad, DeclRad, a, b, c, raster_mode, raster_shape=raster_shape)

    else:
        omega_set_hor, omega_rise_hor, cos_theta_omega_set_hor, cos_theta_omega_rise_hor, sin_omega_rise_pixel, \
        omega_rise_pixel_24, sin_omega_set_pixel, \
        omega_set_pixel_24, lower_int_limit_rise, upper_int_limit_set, quadratic_function = \
            calc_two_integration_limits(LatRad, DeclRad, a, b, c, raster_mode)

    # RUN calc_two_daytime_integration_limits function
    if raster_mode:
        omega_rise_pixel_24, omega_set_pixel_24, lower_int_limit_rise, upper_int_limit_set, \
        omega_set_during_day_pixel_24, omega_rise_during_day_pixel_24, X \
            = calc_two_daytime_integration_limits(omega_rise_pixel_24, omega_set_pixel_24, lower_int_limit_rise,
                                                  upper_int_limit_set, sin_slope, sin_lat, sin_decl, sin_aspect,
                                                  cos_slope,
                                                  cos_lat, cos_decl, cos_aspect, a, b, c, quadratic_function, raster_mode=True)

    else:
        omega_rise_pixel_24, omega_set_pixel_24, lower_int_limit_rise, upper_int_limit_set, \
        omega_set_during_day_pixel_24, omega_rise_during_day_pixel_24, X \
            = calc_two_daytime_integration_limits(omega_rise_pixel_24, omega_set_pixel_24, lower_int_limit_rise,
                                                  upper_int_limit_set, sin_slope, sin_lat, sin_decl, sin_aspect,
                                                  cos_slope,
                                                  cos_lat, cos_decl, cos_aspect, a, b, c, quadratic_function)

    # RUN integral_cos_theta
    if raster_mode:
        int_cos_theta = integral_cos_theta(omega_rise_hor, omega_set_hor, sin_decl, cos_decl, sin_lat, cos_lat,
                                           sin_slope, cos_slope, sin_aspect, cos_aspect, raster_mode)
    else:
        int_cos_theta = integral_cos_theta(omega_rise_hor, omega_set_hor, sin_decl, cos_decl, sin_lat, cos_lat,
                                           sin_slope, cos_slope, sin_aspect, cos_aspect)

    # Radiation from Allen 2006 Eq [1]
    Ra = 1367 * (1 / (1 + 0.033 * np.cos(DOY * 2 * np.pi / 365))) * int_cos_theta

    if not raster_mode:
        print('\n', 'int_cos_theta = {:.5f}'.format(int_cos_theta))
        print('\n', 'omega_rise_hor = {:.5f}'.format(omega_rise_hor))
        print('\n', 'omega_set_hor = {:.5f}'.format(omega_set_hor))
        print('\n', 'Ra = {:.5f}'.format(Ra))


    if raster_mode:
        # for the raster, better now to leave the time decimal as a decimal for conversion later on...

        omega_rise_hor_time_decimal = omega_rise_hor * 180 / np.pi / 15 + 12
        omega_set_hor_time_decimal = omega_set_hor * 180 / np.pi / 15 + 12

        lower_int_limit_rise_time_decimal = lower_int_limit_rise * 180 / np.pi / 15 + 12
        upper_int_limit_set_time_decimal = upper_int_limit_set * 180 / np.pi / 15 + 12

        omega_set_during_day_pixel_24_time_decimal = omega_set_during_day_pixel_24 * 180 / np.pi / 15 + 12
        omega_rise_during_day_pixel_24_time_decimal = omega_rise_during_day_pixel_24 * 180 / np.pi / 15 + 12

    else:
        omega_rise_hor_time_decimal = omega_rise_hor * 180 / np.pi / 15 + 12
        print('the decimal', omega_rise_hor_time_decimal)
        ihours = int(omega_rise_hor_time_decimal)
        # returns omega_rise_hor_time as a tuple of hours and seconds...
        omega_rise_hor_time = ihours, (omega_rise_hor_time_decimal - ihours) * 60
        print('%02d:%02d' % omega_rise_hor_time)

        omega_set_hor_time_decimal = omega_set_hor * 180 / np.pi / 15 + 12
        ihours = int(omega_set_hor_time_decimal)
        omega_set_hor_time = ihours, (omega_set_hor_time_decimal - ihours) * 60
        print('%02d:%02d' % omega_set_hor_time)

        lower_int_limit_rise_time_decimal = lower_int_limit_rise * 180 / np.pi / 15 + 12
        ihours = int(lower_int_limit_rise_time_decimal)
        lower_int_limit_rise_time = ihours, (lower_int_limit_rise_time_decimal - ihours) * 60
        print('%02d:%02d' % lower_int_limit_rise_time)

        upper_int_limit_set_time_decimal = upper_int_limit_set * 180 / np.pi / 15 + 12
        ihours = int(upper_int_limit_set_time_decimal)
        upper_int_limit_set_time = ihours, (upper_int_limit_set_time_decimal - ihours) * 60
        print('%02d:%02d' % upper_int_limit_set_time)

        omega_set_during_day_pixel_24_time_decimal = omega_set_during_day_pixel_24 * 180 / np.pi / 15 + 12
        ihours = int(omega_set_during_day_pixel_24_time_decimal)
        omega_set_during_day_pixel_24_time = ihours, (omega_set_during_day_pixel_24_time_decimal - ihours) * 60
        print('%02d:%02d' % omega_set_during_day_pixel_24_time)

        omega_rise_during_day_pixel_24_time_decimal = omega_rise_during_day_pixel_24 * 180 / np.pi / 15 + 12
        ihours = int(omega_rise_during_day_pixel_24_time_decimal)
        omega_rise_during_day_pixel_24_time = ihours, (omega_rise_during_day_pixel_24_time_decimal - ihours) * 60
        print('%02d:%02d' % omega_rise_during_day_pixel_24_time)


    # Summary of output data for calculation of self-shadowing periods of pixel
    if raster_mode:
        omega_set_during_day_pixel_24[X < 0] = 9.0
        omega_rise_during_day_pixel_24[X < 0] = 9.0

    else:
        print('\n', 'Summary of output data for calculation of self-shadowing periods of pixel')
        print('\n', 'Latitude = {:.4f}'.format(LatDeg), 'degrees', '    ', 'DOY = {:.0f}'.format(DOY), '    ',
              'slope = {:.0f}'.format(SlopeDeg), 'degrees', '     ',
              'aspect = {:.0f}'.format(AspectDeg), 'degrees')
        print('\n', 'omega_rise_hor = {:.4f}'.format(omega_rise_hor), '    ', 'Daylight starts at', '  ',
              '%02d:%02d' % omega_rise_hor_time, '\n', 'omega_set_hor  =  {:.4f}'.format(omega_set_hor), '    ',
              'Daylight ends at', '    ', '%02d:%02d' % omega_set_hor_time, )
        if X < 0:
            print('\n', 'X = {:.5f}'.format(X), 'at X<0 there are two periods of direct beam radiation during the day')
            print('\n', 'lower_int_limit_rise = {:.4f}'.format(lower_int_limit_rise), '    ', 'Beam radiation starts at',
                  '  ',
                  '%02d:%02d' % lower_int_limit_rise_time, '\n',
                  'upper_int_limit_set  =  {:.4f}'.format(upper_int_limit_set), '    ',
                  'Beam radiation ends at', '    ', '%02d:%02d' % upper_int_limit_set_time, )
            print('\n', 'omega_set_during_day_pixel_24   = {:.4f}'.format(omega_set_during_day_pixel_24), '      ',
                  'Beam radiation ends during the day at', '  ', '%02d:%02d' % omega_set_during_day_pixel_24_time, '\n',
                  'omega_rise_during_day_pixel_24  =  {:.4f}'.format(omega_rise_during_day_pixel_24), '      ',
                  'Beam radiation starts again at', '         ', '%02d:%02d' % omega_rise_during_day_pixel_24_time, )
        else:
            print('\n', 'X = {:.5f}'.format(X), 'at X=>0 there is only one period of direct beam radiation during the day')
            omega_set_during_day_pixel_24 = 9.0
            omega_rise_during_day_pixel_24 = 9.0
        #
        print('\n', 'the four integration limits for export are')
        print('\n', 'lower_int_limit_rise          = {:.4f}'.format(lower_int_limit_rise), '            ',
              ' upper_int_limit_set  =  {:.4f}'.format(upper_int_limit_set))
        print(' omega_set_during_day_pixel_24 = {:.4f}'.format(omega_set_during_day_pixel_24), '             ',
              'omega_rise_during_day_pixel_24  =  {:.4f}'.format(omega_rise_during_day_pixel_24))

    # =========== 2D RASTER OUTPUT =================
    if raster_mode:
        # Now we output the images as rasters.
        # These are the images that will be output as rasters for Reasearch and Development of future scripts.
        # Jan, add more here if I have left any out.
        vars_to_output = {'omega_rise_time_horizontal': omega_rise_hor_time_decimal,
                          'omega_set_time_horizonatl': omega_set_hor_time_decimal,
                          'lower_limit_rise_time': lower_int_limit_rise_time_decimal,
                          'upper_limit_rise_time': upper_int_limit_set_time_decimal,
                          'omega_set_daytime': omega_set_during_day_pixel_24_time_decimal,
                          'omega_rise_daytime': omega_rise_during_day_pixel_24_time_decimal,
                          'int_cos_theta': int_cos_theta,
                          }

        for k, v in vars_to_output.items():
            filename = '{}.tif'.format(k)
            outpath = os.path.join(output_dir, filename)
            arr = np.reshape(v, raster_shape)
            convert_array_to_raster(outpath, arr, geo)



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

    # =============== RASTER MODE INPUTS (MSEC COMP) =======================

    # We may not need the DEM.
    dem_path = '/Users/dcadol/Desktop/academic_docs_II/m_mountain_metric/Los_Lunas_AOI/dem_clipped/los_lunas_aoi_dem_clipped.tif'
    # RASTER EQUIVALENT TO SlopeDeg
    slope_path = '/Users/dcadol/Desktop/academic_docs_II/m_mountain_metric/Los_Lunas_AOI/dem_clipped/los_lunas_aoi_dem_clipped_slope.tif'
    # RASTER EQUIVALENT TO AspectDeg
    aspect_path = '/Users/dcadol/Desktop/academic_docs_II/m_mountain_metric/Los_Lunas_AOI/dem_clipped/los_lunas_aoi_dem_clipped_aspect.tif'
    # RASTER EQUIVALENT TO LatDeg
    latitude_path = '/Users/dcadol/Desktop/academic_docs_II/m_mountain_metric/Los_Lunas_AOI/dem_clipped/los_lunas_aoi_dem_clipped_lat.tif'

    # OUTPUT
    output_dir = '/Users/dcadol/Desktop/academic_docs_II/m_mountain_metric/Los_Lunas_AOI/script_test_output'

    # We make it optional to run a raster mode or a 1-D mode.
    raster_mode = True

    if raster_mode:
        # # The way the function would be called to run the function in raster mode
        calculate_integration_limits(DOY, latitude_path, slope_path, aspect_path, local_time, Lz, Lm, output_dir, raster_mode)
    else:
        calculate_integration_limits(DOY, LatDeg, SlopeDeg, AspectDeg, local_time, Lz, Lm, raster_mode)