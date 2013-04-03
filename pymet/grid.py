# coding: utf-8
u"""
===============================================
格子点データ計算モジュール (:mod:`pymet.grid`)
===============================================

差分計算
========

.. autosummary::

   dvardx
   dvardy
   dvardp
   d2vardx2
   d2vardy2
   dvardvar
   div
   rot

積分
====

.. autosummary::

   vint
   
内挿
====

.. autosummary::

   vinterp


-----------------------   
"""

import numpy as np
import constants as constants
import scipy.interpolate
import _internal

NA=np.newaxis
a0=constants.earth_radius
g = constants.earth_gravity
PI = constants.pi
d2r=PI/180.

__all__ = ['dvardx', 'dvardy', 'dvardp', 'd2vardx2', 'd2vardy2', 'div', 'rot', 'dvardvar',
           'vint',
           'vinterp']

#=== 微分と差分 ====================================================================================

def dvardx(var, lon, xdim, cyclic=True):
    u"""
    経度方向の微分を中央差分で計算。

    :Arguments:
     **var** : ndarray
       微分を計算する領域の格子点の値
     **lon** : array_like
       経度
     **xdim**: int
       経度度次元のインデックス。len(var.shape[xdim]) == len(lon)でなければならない。
     **cyclic** : bool, optional
       東西の境界を周期境界として扱うかどうか。False の場合は周期境界を用いずに前方、後方差分
       で計算する。デフォルトは True。

       
    :Returns:
     **result** : ndarray
       varと同じ形状。

       
    **Examples**
    
    >>> var.shape
    (24, 73, 144)
    >>> lon = np.arange(0, 360, 2.5)
    >>> result = dvardx(var, lon, 2, cyclic=True)
    >>> result.shape
    (24, 73, 144)

    領域が全球でない場合。
    
    >>> var.shape
    (24, 73, 72)
    >>> lon = np.arange(0, 180, 2.5)
    >>> result = dvardx(var, lon, 2, cyclic=False)
    >>> result.shape
    (24, 73, 72)
    """ 
    var = np.array(var)
    ndim = var.ndim
    var = np.rollaxis(var,xdim,ndim)
    if cyclic:
        dvar = np.concatenate(\
                    ((var[...,1]-var[...,-1])[...,NA],\
                    (var[...,2:]-var[...,:-2]),\
                    (var[...,0]-var[...,-2])[...,NA]), axis=-1)
        dx   = PI/180. * \
               np.r_[(lon[1]+360-lon[-1]),\
                    (lon[2:]-lon[:-2]   ),\
                    (lon[0]+360-lon[-2] )]
    else:
        dvar = np.concatenate(\
                    ((var[...,1]-var[...,0])[...,NA],\
                    (var[...,2:]-var[...,:-2]),\
                    (var[...,-1]-var[...,-2])[...,NA]), axis=-1)
        dx   = PI/180. * \
               np.r_[(lon[1]-lon[0] ),\
                    (lon[2:]-lon[:-2]   ),\
                    (lon[-1]-lon[-2])]
    out = dvar/dx/a0
    out = np.rollaxis(out,ndim-1,xdim)

    return out

def dvardy(var, lat, ydim):
    u"""
    緯度方向の微分を中央差分で計算。南北端は前方、後方差分

    :Arguments:
     **var** : ndarray
       微分を計算する領域の格子点の値
     **lat** : array_like
       緯度
     **ydim**: int
       緯度次元のインデックス。len(var.shape[ydim]) == len(lat)でなければならない。

       
    :Returns:
     **result** : ndarray
       varと同じ形状。

       
    **Examples**
    
    >>> var.shape
    (24, 73, 144)
    >>> lat = np.arange(-90, 90.1, 2.5)
    >>> result = dvardy(var, lat, 1)
    >>> result.shape
    (24, 73, 144)
    """ 
    var = np.array(var)
    ndim = var.ndim
    var = np.rollaxis(var,ydim,ndim)

    dvar = np.concatenate(\
              [(var[...,1] -var[...,0]  )[...,NA],\
               (var[...,2:]-var[...,:-2]),\
               (var[...,-1]-var[...,-2] )[...,NA]], axis=-1) 
    dy   = PI/180. * \
           np.r_[(lat[1]-lat[0]),\
                 (lat[2:]-lat[:-2]),\
                 (lat[-1]-lat[-2])]

    out = dvar/dy/a0
    out = np.rollaxis(out,ndim-1,ydim)

    return out

def dvardp(var, lev, zdim, punit=100.):
    u"""
    鉛直方向の微分をlog(p)の中央差分で計算。上下端は前方、後方差分

    :Arguments:
     **var**   : ndarray
       微分を計算する領域の格子点の値
     **lev**   : 1d-array
       等圧面のレベル
     **zdim**  : int 
       鉛直次元のインデックス。len(var.shape[zdim]) == len(lev)でなければならない。
     **punit** : float, optional
       圧力座標の単位(Pa)。デフォルトは100、すなわちhPaとして扱う。

       
    :Returns:
     **dvardp** : ndarray
       varと同じ形状をもつ。

       
    **Examples**
    
    >>>
    >>>
    """ 
    var = np.array(var)
    ndim = var.ndim
    lev = lev * punit
        
    #roll lat dim axis to last
    var = np.rollaxis(var,zdim,ndim)
    dvar = np.concatenate(\
              [ (var[...,1] -var[...,0]  )[...,NA],\
                (var[...,2:]-var[...,:-2]),\
                (var[...,-1]-var[...,-2] )[...,NA] ], axis=-1) 
    dp   = np.r_[ np.log(lev[1]/lev[0])*lev[0] ,\
                  np.log(lev[2:]/lev[:-2])*lev[1:-1],\
                  np.log(lev[-1]/lev[-2])*lev[-1] ]
    
    out = dvar/dp
    #reroll lat dim axis to original dim
    out = np.rollaxis(out,ndim-1,zdim)

    return out

def d2vardx2(var, lon, xdim, cyclic=True):
    u"""
    経度方向の2階微分を中央差分で計算。

    :Arguments:
     **var** : ndarray
       微分を計算する領域の格子点の値
     **lon** : array_like
       経度
     **xdim** : int 
       経度次元のインデックス。len(var.shape[xdim]) == len(lon)でなければならない。
     **cyclic** : bool, optional
       経度方向にはサイクリックとするかどうか。デフォルトではTrue。Falseの場合は東西端はゼロとする。


    :Returns:
     **result** : ndarray
       varと同じ形状をもつ。


    **Examples**
    
    >>>
    >>>
    """ 
    var = np.array(var)
    ndim = var.ndim

    #roll lon dim axis to last
    var = np.rollaxis(var,xdim,ndim)

    if cyclic:
        dvar = np.concatenate(\
               ((var[...,1]-2*var[...,0]+var[...,-1])[...,NA],\
                (var[...,2:]-2*var[...,1:-1]+var[...,:-2]),\
                (var[...,0]-2*var[...,-1]+var[...,-2])[...,NA]),\
                axis=-1)
        dx   = PI/180. * \
               np.r_[(lon[1]+360-lon[-1]),\
                    (lon[2:]-lon[:-2]   ),\
                    (lon[0]+360-lon[-2] )]
    else: #edge is zero
        dvar = np.concatenate(\
               ((var[...,0]-var[...,0])[...,NA],\
                (var[...,2:]-2*var[...,1:-1]+var[...,:-2]),\
                (var[...,0]-var[...,0])[...,NA]),\
                axis=-1)
        dx   = PI/180. * \
               np.r_[(lon[1]-lon[0] ),\
                    (lon[2:]-lon[:-2]   ),\
                    (lon[-1]-lon[-2])]
    out = 4.*dvar/dx/dx/a0/a0
    #reroll lon dim axis to original dim
    out = np.rollaxis(out,ndim-1,xdim)

    return out

def d2vardy2(var, lat, ydim):
    u"""
    緯度方向の2階微分を中央差分で計算。南北端はゼロとする。

    .. todo:: 南北端をマスクするようなオプションの追加。 
    
    :Arguments:
     **var**  : ndarray
       微分を計算する領域の格子点の値
     **lat**  : array_like
       緯度
     **ydim** : int 
       緯度次元のインデックス。len(var.shape[ydim]) == len(lat)でなければならない。

       
    :Returns:
     **result** : ndarray
       varと同じ形状をもつ。

       
    **Examples**
    
    >>>
    >>>
    """ 
    var = np.array(var)
    ndim = var.ndim

    #roll lat dim axis to last
    var = np.rollaxis(var,ydim,ndim)

    #edge is zero
    dvar = np.concatenate(\
              [(var[...,0] -var[...,0]  )[...,NA],\
               (var[...,2:]-2*var[...,1:-1]+var[...,:-2]),\
               (var[...,0]-var[...,0] )[...,NA]], axis=-1) 
    dy   = PI/180. * \
           np.r_[(lat[1]-lat[0]),\
                 (lat[2:]-lat[:-2]),\
                 (lat[-1]-lat[-2])]

    out = 4.*dvar/dy/dy/a0/a0
    #reroll lat dim axis to original dim
    out = np.rollaxis(out,ndim-1,ydim)

    return out

def div(u, v, lon, lat, xdim, ydim, cyclic=True):
    u"""
    水平発散を計算する。

    :Arguments:
     **u, v** : ndarray
     
     **lon, lat** : array_like
       緯度と経度
     **xdim, ydim** : int
       緯度、経度の軸
     **cyclic** : bool, optional
       経度微分の際に周境界とするかどうか。デフォルトはTrue

    
    :Returns:
     **div** : ndarray
       u, v と同じ形状のndarray

       
    **Notes**
    
    :math:`\partial u/\partial x+\partial v/\partial y`


    **Examples**

    >>> from pymet.grid import div
    >>> 
    """
    u, v = np.array(u), np.array(v)
    ndim = u.ndim

    out = dvardx(u,lon,xdim,cyclic=cyclic) + dvardy(v,lat,ydim)

    return out

def rot(u, v, lon, lat, xdim, ydim, cyclic=True):
    u"""
    回転を計算
    """
    u, v = np.array(u), np.array(v)
    ndim = u.ndim

    out = dvardx(v,lon,xdim,cyclic=cyclic) - dvardy(u,lat,ydim)

    return out
    
def dvardvar(var1, var2, dim, cyclic=True):
    u"""
    d(var1)/d(var2) を axis=dim に沿って差分を計算
    """
    var1, var2 = np.array(var1), np.array(var2)
    ndim = var1.ndim
        
    #roll dim axis to last
    var1 = np.rollaxis(var1,dim,ndim)
    var2 = np.rollaxis(var2,dim,ndim)

    dvar1 = np.concatenate(\
              [ (var1[...,1] -var1[...,0]  )[...,NA],\
                (var1[...,2:]-var1[...,:-2]),\
                (var1[...,-1]-var1[...,-2] )[...,NA] ], axis=-1) 
    dvar2 = np.concatenate(\
              [ (var2[...,1] -var2[...,0]  )[...,NA],\
                (var2[...,2:]-var2[...,:-2]),\
                (var2[...,-1]-var2[...,-2] )[...,NA] ], axis=-1) 

    out = dvar1/dvar2
    #reroll lat dim axis to original dim
    out = np.rollaxis(out,ndim-1,dim)

    return out

#=============================================================================

def vint(var, bottom, top, lev, zdim, punit=100.):
    u"""
    質量重み付き鉛直積分。
    vint = 
    
    :Arguments:
     **var** : array_like
       データ。       
     **bottom, top** : float
       積分の下端、上端。      
     **lev** : 1d-array
       等圧面のレベル
     **zdim** : int
       鉛直次元のインデックス。
     **punit** : float, optional
       levの単位(Pa)

       
    :Returns:
     **result** : array_like
       入力配列から鉛直次元を除いた形状と同じ。

       
    **Examples**

    >>>
    >>>
    """
    var = np.array(var)
    ndim = var.ndim

    lev = np.array(lev)[(lev <= bottom)&(lev >= top)]
    lev_m = np.r_[bottom,(lev[1:] + lev[:-1])/2.,top]
    dp = lev_m[:-1] - lev_m[1:]

    #roll lat dim axis to last
    var = np.rollaxis(var,zdim,ndim)
    out = var[...,(lev <= bottom)&(lev >= top)] * dp / g * punit
    out = out.sum(axis=-1)

    return out

def vinterp(var, oldz, newz, zdim, logintrp=True, bounds_error=True):
    u"""
    指定した鉛直レベルへの線形内挿。
    
    :Arguments:
     **var** : array_like
       内挿するデータ。
     **old_z** : array_like
       元データの鉛直レベル。
     **new_z** : array_like
       内挿するレベル。
     **zdim** : int
       鉛直次元のインデックス。
     **logintrp** : bool, optional
       log(z)で線形内挿するかどうか。デフォルトはTrue
     **bounds_error** : bool, optional


    :Returns:
     **result** : array_like
       内挿後のデータ。result.shape[zdim] == len(new_z)になる。

       
    **Examples**

    >>> from pymet.grid import vinterp
    >>>
        
    """
    var = np.array(var)
    ndim = var.ndim    

    new_z = np.array(newz)
    old_z = np.array(oldz)
    if logintrp:
        old_z = np.log(old_z)
        new_z = np.log(new_z)
    old_zn = var.shape[zdim]
    new_zn = len(new_z)
    
    #roll z dim axis to last    
    var = np.rollaxis(var,zdim,ndim)
    old_shape = var.shape
    new_shape = list(old_shape)
    new_shape[-1] = new_zn
    var = var.reshape(-1,old_zn)
    if old_z.ndim == ndim:
        out = np.empty((var.shape[0],new_zn))
        old_z = np.rollaxis(old_z,zdim,ndim).reshape(-1,old_zn)
        out = _internal.linear_interp(var, old_z, new_z)
    elif old_z.ndim == 1:
        f = scipy.interpolate.interp1d(old_z,var,kind='linear',bounds_error=bounds_error)
        out = f(new_z)

    #reroll lon dim axis to original dim
    out = out.reshape(new_shape)
    out = np.rollaxis(out,ndim-1,zdim)

    return out



def distance(lon1,lon2,lat1,lat2):
    u"""
    2点間の距離を計算する。

    :Arguments:
     **lon1, lon2** : float or ndarray of floats
        2地点の経度(degrees)
     **lat1, lat2** : float or ndarray of floats
        2地点の緯度(degrees)

        
    :Returns:
     **result** : floats or ndarray of floats
        (lon1, lat1), (lon2, lat2) 間の距離 (km)


    **Examples**

    >>> from pymet.grid import distance
    >>> distance(100., 180., 30., 40.)
    
    """
    return a0 * np.arccos(np.sin(lat1*d2r)*np.sin(lat2*d2r)
                        +np.cos(lat1*d2r)*np.cos(lat2*d2r)*np.cos((lon2-lon1)*d2r))

