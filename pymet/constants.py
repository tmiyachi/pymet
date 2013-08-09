# coding: utf-8
u"""
=========================================
定数モジュール (:mod:`pymet.constants`)
=========================================

物理・数学定数


数学定数
--------

============= =======================================================
``pi``        PI
============= =======================================================

地球科学定数
------------

================== ==================================================
``earth_radius``   地球半径     [m]
``earth_daysec``   地球自転周期 [sec]
``earth_omega``    地球回転角速度 [rad/sec]
``earth_gravity``  重力加速度 [m s^-2]
================== ==================================================

.. autosummary::

   earth_f
   earth_beta

気体に関する定数
-----------------

================== ============================================================================================
``air_rd``         乾燥空気の気体定数 :math:`R_{d}` [J K^-1 Kg^-1]
``air_rv``         水蒸気の気体定数   :math:`R_{v}` [J K^-1 Kg^-1]
``air_cp``         乾燥空気の定圧比熱 :math:`C_{p}` [J K^-1 Kg^-1]
``air_cv``         乾燥空気の定積比熱 :math:`C_{v}` [J K^-1 Kg^-1]
``air_gamma``      比熱比 :math:`\gamma=C_{p}/C_{v}`
``air_kappa``      :math:`\kappa=R_{d}/C_{p}`
``air_eps``        :math:`\epsilon=R_{d}/R_{v}`
================== ============================================================================================

-------
"""

import math as _math
import numpy as _np
#constant

#math_constants
pi  = _math.pi
d2r = _np.deg2rad   # degree => radian
r2d = _np.rad2deg   # radian => degree

#Geophysical constants
earth_radius = 6.371e+6
earth_daysec = 86400.0
earth_omega  = 2.0*pi/earth_daysec
earth_gravity = 9.80665

#air_constants
air_rd = 287.0
air_rv = 461.0
air_cp = 1004.0
air_cv = 717.0
air_kappa = air_rd/air_cp
air_gamma = air_cp/air_cv
air_eps = air_rd/air_rv

#functions for geophysical constants
def earth_f(lat):
    u"""
    コリオリパラメータを緯度から計算する。

    :Parameters:
     **lat** : array_like
       コリオリパタメータを求める緯度(度)

    :Returns:    
     **f** : float or array of floats
       コリオリパラメータ

    .. note::
     Computes ``f = 2 * earth_omega*sin(lat)``

    **Examples**
     >>> from pymet.constants.mod_const import earth_f
     >>> earth_f(__np.array([-30., 0., 30.))
     array()
    """
    return 2.*earth_omega*_np.sin(pi/180.*lat)

def earth_beta(lat):
    u"""
    ベータパラメータを緯度から計算する。

    :Parameters:
     **lat** : array_like
          latitude (degrees) to be converted

    :Returns:    
     **beta** : float or array of floats
           beta parameter

    **Notes**
     Computes ``beta = df/dy = 2*earth_omega/earth_radius*cos(lat)``

    **Examples**
     >>> from pymet.constants.mod_const import earth_beta
     >>> earth_beta(__np.array([-30., 0., 30.))
     array()
    """
    return 2.*earth_omega/earth_radius*_np.cos(pi/180.*lat)

def dlon2dx(dlon, clat):
    u"""
    経度差から局所直交系での実距離を求める。
    """
    return d2r(dlon)*earth_radius*_np.cos(d2r(clat))

def dx2dlon(dx, clat):
    u"""
    距離から経度差を求める。
    """
    return r2d(dx/earth_radius/_np.cos(d2r(clat)))

def dlat2dy(dlat):
    u"""
    緯度差から局所直交系での実距離を求める。
    """
    return d2r(dlat)*earth_radius

def dy2dlat(dy):
    u"""
    距離から緯度差を求める。
    """
    return r2d(dy/earth_radius)
