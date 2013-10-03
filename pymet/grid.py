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
   div
   rot
   grad
   skgrad
   laplacian

積分
====

.. autosummary::

   vint
   
内挿
====

.. autosummary::

   vinterp

その他
======

.. autosummary::
   distance
   
-----------------------   
"""

import numpy as np
import constants as constants
import scipy.interpolate
import _internal
import tools

NA=np.newaxis
a0=constants.earth_radius
g = constants.earth_gravity
PI = constants.pi
d2r=PI/180.

__all__ = ['dvardx', 'dvardy', 'dvardp', 'd2vardx2', 'd2vardy2', 'div', 'rot', 'grad', 'skgrad', 'laplacian', #'dvardvar',
           'vint', 'vmean',
           'vinterp',
           'distance']

#=== 微分と差分 ====================================================================================

def dvardx(var, lon, lat, xdim, ydim, cyclic=True, sphere=True):
    ur"""
    経度方向のx微分を中央差分で計算。

    :Arguments:
     **var** : ndarray
       微分を計算する領域の格子点の値
     **lon** : array_like
       経度, もしくはx座標
     **lat** : array_like
       緯度, もしくはy座標
     **xdim, ydim**: int
       入力配列の経度、緯度次元のインデックス
     **cyclic** : bool, optional
       東西の境界を周期境界として扱うかどうか。False の場合は周期境界を用いずに前方、後方差分
       で計算する。デフォルトは True。
     **sphere** : bool, optional
       球面緯度経度座標かどうか。デフォルトはTrue。Falseにすると直交座標として扱う。
       
    :Returns:
     **result** : ndarray
       varと同じ形状。

    .. note::
       球面緯度経度座標系におけるx微分は、次のように定義される。

       .. math:: \frac{\partial \Phi}{\partial x} = \frac{1}{a\cos\phi}\frac{\partial \Phi}{\partial \lambda}

       ここで、aは地球半径。これを中央差分で次のように計算する。
       
       .. math:: \left( \frac{\partial \Phi}{\partial x} \right)_{i,j}
                  = \frac{1}{a\cos\phi_{j}}\frac{\Phi_{i+1,j} - \Phi_{i-1,j}}{\lambda_{i+1} - \lambda_{i-1}}

       cyclic=Falseの場合は、両端は前方、後方差分を用いて、

       .. math:: \left( \frac{\partial \Phi}{\partial x} \right)_{0,j}
               = \frac{1}{a\cos\phi_{j}}\frac{\Phi_{1,j} - \Phi_{0,j}}{\lambda_{1} - \lambda_{0}}
                 \hspace{3em}
                 \left( \frac{\partial \Phi}{\partial x} \right)_{N-1,j}
               = \frac{1}{a\cos\phi_{j}}\frac{\Phi_{N-1,j} - \Phi_{N-2,j}}{\lambda_{N-1} - \lambda_{N-2}}

       で計算する。sphere=Falseの場合は、
       
       .. math:: \left( \frac{\partial \Phi}{\partial x} \right)_{i,j}
                  = \frac{\Phi_{i+1,j} - \Phi_{i-1,j}}{x_{i+1} - x_{i-1}}


               
    **Examples**    
     >>> var.shape
     (24, 73, 144)
     >>> lon = np.arange(0, 360, 2.5)
     >>> lat = np.arange(-90, 90.1, 2.5)
     >>> result = dvardx(var, lon, lat, 2, 1, cyclic=True)
     >>> result.shape
     (24, 73, 144)

     領域が全球でない場合。
    
     >>> var.shape
     (24, 73, 72)
     >>> lon = np.arange(0, 180, 2.5)
     >>> lat = np.arange(-90, 90.1, 2.5)
     >>> result = dvardx(var, lon, lat, 2, 1, cyclic=False)
     >>> result.shape
     (24, 73, 72)
    """ 
    var = np.array(var)
    ndim = var.ndim
    var = np.rollaxis(var,xdim,ndim)
    if cyclic and sphere:
        dvar = np.concatenate(\
                    ((var[...,1]-var[...,-1])[...,NA],\
                    (var[...,2:]-var[...,:-2]),\
                    (var[...,0]-var[...,-2])[...,NA]), axis=-1)
        dx   = np.r_[(lon[1]+360-lon[-1]),\
                    (lon[2:]-lon[:-2]   ),\
                    (lon[0]+360-lon[-2] )]
    else:
        dvar = np.concatenate(\
                    ((var[...,1]-var[...,0])[...,NA],\
                    (var[...,2:]-var[...,:-2]),\
                    (var[...,-1]-var[...,-2])[...,NA]), axis=-1)
        dx   = np.r_[(lon[1]-lon[0] ),\
                    (lon[2:]-lon[:-2]   ),\
                    (lon[-1]-lon[-2])]
                    
    dvar = np.rollaxis(dvar,ndim-1,xdim)
    if sphere:
        dx = a0*PI/180.*tools.expand(dx,ndim,xdim) * tools.expand(np.cos(lat*d2r),ndim,ydim)
    else:
        dx = tools.expand(dx,ndim,xdim)
    out = dvar/dx
    
    return out

def dvardy(var, lat, ydim, sphere=True):
    ur"""
    緯度方向のy微分を中央差分で計算。南北端は前方、後方差分

    :Arguments:
     **var** : ndarray
       微分を計算する領域の格子点の値
     **lat** : array_like
       緯度、もしくはy座標
     **ydim**: int
       緯度次元のインデックス。len(var.shape[ydim]) == len(lat)でなければならない。       
     **sphere** : bool, optional
       球面緯度経度座標かどうか。デフォルトはTrue。Falseにすると直交座標として扱う。

    :Returns:
     **result** : ndarray
       varと同じ形状。

    .. note::
       球面緯度経度座標系におけるy微分は、次のように定義される。

       .. math:: \frac{\partial \Phi}{\partial y} = \frac{1}{a}\frac{\partial \Phi}{\partial \phi}

       ここで、aは地球半径。これを中央差分で次のように計算する。
       
       .. math:: \left( \frac{\partial \Phi}{\partial y} \right)_{i,j}
                  = \frac{1}{a}\frac{\Phi_{i,j+1} - \Phi_{i,j-1}}{\phi_{j+1} - \phi_{j-1}}

       南北端は前方、後方差分を用いて、

       .. math:: \left( \frac{\partial \Phi}{\partial y} \right)_{i,0}
                 = \frac{1}{a}\frac{\Phi_{i,1} - \Phi_{i,0}}{\phi_{1} - \phi_{0}}
                   \hspace{3em}
                   \left( \frac{\partial \Phi}{\partial y} \right)_{i,N-1}
                 = \frac{1}{a}\frac{\Phi_{i,N-1} - \Phi_{i,N-2}}{\phi_{N-1} - \phi_{N-2}}

       で計算する。
       
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
    dy   = np.r_[(lat[1]-lat[0]),\
                 (lat[2:]-lat[:-2]),\
                 (lat[-1]-lat[-2])]

    if sphere:
        dy = a0*PI/180.*dy
    out = dvar/dy
    out = np.rollaxis(out,ndim-1,ydim)

    return out

def dvardp(var, lev, zdim, punit=100.):
    ur"""
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

    .. note::
       気圧座標系におけるp微分は、次のように書き換えられる。

       .. math:: \frac{\partial \Phi}{\partial p} = \frac{1}{p}\frac{\partial \Phi}{\partial \ln{p}}

       したがって、中央差分では、
       
       .. math:: \left( \frac{\partial \Phi}{\partial p} \right)_{k}
                  = \frac{1}{p_{k}}\frac{\Phi_{k+1} - \Phi_{k-1}}{\ln{p_{k+1}} - \ln{p_{k-1}}}

       とかける。上下端は前方、後方差分を用いて、

       .. math:: \left( \frac{\partial \Phi}{\partial p} \right)_{0}
                     = \frac{1}{p_{0}}\frac{\Phi_{1} - \Phi_{0}}{\ln{p_{1}} - \ln{p_{0}}}
                 \hspace{3em}
                 \left( \frac{\partial \Phi}{\partial p} \right)_{N-1}
                     = \frac{1}{p_{N-1}}\frac{\Phi_{N-1} - \Phi_{N-2}}{\ln{p_{N-1}} - \ln{p_{N-2}}}

       で計算する。
       
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

# not compleate----------------------------------------------------------------------------------------------------

def d2vardx2(var, lon, lat, xdim, ydim, cyclic=True, sphere=True):
    ur"""
    経度方向の2階x微分を中央差分で計算。

    :Arguments:
     **var** : ndarray
       微分を計算する領域の格子点の値
     **lon** : array_like
       経度、もしくはx座標
     **lat** : array_like
       緯度、もしくはy座標
     **xdim, ydim** : int 
       経度、緯度次元のインデックス。
     **cyclic** : bool, optional
       経度方向にはサイクリックとするかどうか。デフォルトではTrue。Falseの場合は東西端はゼロとする。
     **sphere** : bool, optional
       球面緯度経度座標かどうか。デフォルトはTrue。Falseにすると直交座標として扱う。    

    :Returns:
     **result** : ndarray
       varと同じ形状をもつ。

    .. note::
       球面緯度経度座標系におけるx2回微分は、次のように定義される。

       .. math:: \frac{\partial^2 \Phi}{\partial x^2} = \frac{1}{a^2\cos^2\phi}\frac{\partial^2 \Phi}{\partial \lambda^2}

       ここで、aは地球半径。これを中央差分で次のように計算する。
       
       .. math:: \left( \frac{\partial^2 \Phi}{\partial x^2} \right)_{i,j}
                  = \frac{4}{a^2\cos^2\phi_{j}}\frac{\Phi_{i+1,j} - 2\Phi_{i,j} + \Phi_{i-1,j}}{\lambda_{i+1} - \lambda_{i-1}}

       cyclic=Falseの場合は、両端はゼロとする。

    **Examples**    
     >>>
     >>>
    """ 
    var = np.array(var)
    ndim = var.ndim

    #roll lon dim axis to last
    var = np.rollaxis(var,xdim,ndim)

    if cyclic and sphere:
        dvar = np.concatenate(\
               ((var[...,1]-2*var[...,0]+var[...,-1])[...,NA],\
                (var[...,2:]-2*var[...,1:-1]+var[...,:-2]),\
                (var[...,0]-2*var[...,-1]+var[...,-2])[...,NA]),\
                axis=-1)
        dx   = np.r_[(lon[1]+360-lon[-1]),\
                    (lon[2:]-lon[:-2]   ),\
                    (lon[0]+360-lon[-2] )]
    else: #edge is zero
        dvar = np.concatenate(\
               ((var[...,0]-var[...,0])[...,NA],\
                (var[...,2:]-2*var[...,1:-1]+var[...,:-2]),\
                (var[...,0]-var[...,0])[...,NA]),\
                axis=-1)
        dx   = np.r_[(lon[1]-lon[0] ),\
                    (lon[2:]-lon[:-2]   ),\
                    (lon[-1]-lon[-2])]

    dvar = np.rollaxis(dvar,ndim-1,xdim)
    if sphere:
        dx2   = a0**2 * (PI/180.)**2 * tools.expand(dx**2,ndim,xdim) * tools.expand(np.cos(lat*d2r)**2,ndim,ydim)
    else:
        dx2   = tools.expand(dx**2,ndim,xdim)
    out = 4.*dvar/dx2
    #reroll lon dim axis to original dim
    out = np.rollaxis(out,ndim-1,xdim)

    return out

def d2vardy2(var, lat, ydim, sphere=True):
    ur"""
    緯度方向の2階微分を中央差分で計算。南北端はゼロとする。

    .. todo:: 南北端をマスクするようなオプションの追加。 
    
    :Arguments:
     **var**  : ndarray
       微分を計算する領域の格子点の値
     **lat**  : array_like
       緯度
     **ydim** : int 
       緯度次元のインデックス。len(var.shape[ydim]) == len(lat)でなければならない。
     **sphere** : bool, optional
       球面緯度経度座標かどうか。デフォルトはTrue。Falseにすると直交座標として扱う。
       
    :Returns:
     **result** : ndarray
       varと同じ形状をもつ。

    .. note::
       球面緯度経度座標系におけるy2回微分は、次のように定義される。

       .. math:: \frac{\partial^2 \Phi}{\partial y^2} = \frac{1}{a^2}\frac{\partial^2 \Phi}{\partial \phi^2}

       ここで、aは地球半径。これを中央差分で次のように計算する。
       
       .. math:: \left( \frac{\partial^2 \Phi}{\partial y^2} \right)_{i,j}
                  = \frac{4}{a^2}\frac{\Phi_{i,j+1} - 2\Phi_{i,j} + \Phi_{i,j-1}}{\lambda_{j+1} - \lambda_{j-1}}

       cyclic=Falseの場合は、両端はゼロとする。
       
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
    dy   = np.r_[(lat[1]-lat[0]),\
                 (lat[2:]-lat[:-2]),\
                 (lat[-1]-lat[-2])]

    if sphere:
        dy2 = a0**2 * dy**2
    else:
        dy2 = dy**2
    out = 4.*dvar/dy2
    #reroll lat dim axis to original dim
    out = np.rollaxis(out,ndim-1,ydim)

    return out

def div(u, v, lon, lat, xdim, ydim, cyclic=True, sphere=True):
    ur"""
    水平発散を計算する。

    :Arguments:
     **u, v** : ndarray
       ベクトルの東西、南北成分。
     **lon, lat** : array_like
       緯度と経度
     **xdim, ydim** : int
       緯度、経度の軸
     **cyclic** : bool, optional
       経度微分の際に周境界とするかどうか。デフォルトはTrue
     **sphere** : bool, optional
       球面緯度経度座標かどうか。デフォルトはTrue。Falseにすると直交座標として扱う。
    
    :Returns:
     **div** : ndarray
       u, v と同じ形状のndarray
       
    .. note::
       球面緯度経度座標系における水平発散は次のように定義される。

       .. math:: \mbox{\boldmath $\nabla$}\cdot\mbox{\boldmath $v$}
                    = \frac{1}{a\cos\phi} \left[ \frac{\partial u}{\partial \lambda} +\frac{\partial (v\cos\phi)}{\partial \phi} \right] 
                    = \frac{1}{a\cos\phi}\frac{\partial u}{\partial \lambda} + \frac{1}{a}\frac{\partial v}{\partial \phi} - \frac{v\tan\phi}{a}

       ここで、aは地球半径。これを中央差分で次のように計算する。
       
       .. math:: (\mbox{\boldmath $\nabla$}\cdot\mbox{\boldmath $v$})_{i,j}
                    = \left(\frac{1}{a\cos\phi}\frac{\partial u}{\partial \lambda}\right)_{i,j}
                        + \left(\frac{1}{a}\frac{\partial v}{\partial \phi}\right)_{i,j}
                        - \frac{v_{i,j}\tan\phi_{j}}{a}

       緯度、経度微分の項は、:py:func:`dvardx`、 :py:func:`dvardy` を内部で呼ぶ。境界の扱いもこれらに準ずる。
    
    **Examples**
     >>> from pymet.grid import div
     >>> 
    """
    u, v = np.array(u), np.array(v)
    ndim = u.ndim
    
    out = dvardx(u,lon,lat,xdim,ydim,cyclic=cyclic,sphere=sphere) + dvardy(v,lat,ydim,sphere=sphere)
    if sphere:
       out = out - v*tools.expand(np.tan(lat*d2r),ndim,ydim)/a0

    out = np.rollaxis(out, ydim, 0)
    out[0,...]  = 0.
    out[-1,...] = 0.
    out = np.rollaxis(out, 0, ydim+1)
    
    return out

def rot(u, v, lon, lat, xdim, ydim, cyclic=True, sphere=True):
    ur"""
    回転の鉛直成分を計算する。

    :Arguments:
     **u, v** : ndarray
       ベクトルの東西、南北成分。
     **lon, lat** : array_like
       緯度と経度
     **xdim, ydim** : int
       緯度、経度の軸
     **cyclic** : bool, optional
       経度微分の際に周境界とするかどうか。デフォルトはTrue
     **sphere** : bool, optional
       球面緯度経度座標かどうか。デフォルトはTrue。Falseにすると直交座標として扱う。
    
    :Returns:
     **div** : ndarray
       u, v と同じ形状のndarray
       
    .. note::
       球面緯度経度座標系における回転の鉛直成分は次のように定義される。

       .. math:: \mbox{\boldmath $k$}\cdot(\mbox{\boldmath $\nabla$}\times\boldmath{v})
                    &= \frac{1}{a\cos\phi} \left[ \frac{\partial v}{\partial \lambda} - \frac{\partial (u\cos\phi)}{\partial \phi} \right] \\
                    &= \frac{1}{a\cos\phi}\frac{\partial v}{\partial \lambda} - \frac{1}{a}\frac{\partial u}{\partial \phi} + \frac{u\tan\phi}{a}

       ここで、aは地球半径。これを中央差分で次のように計算する。
       
       .. math:: [\mbox{\boldmath $k$}\cdot(\mbox{\boldmath $\nabla$}\times\boldmath{v})]_{i,j}
                    = \left(\frac{1}{a\cos\phi}\frac{\partial v}{\partial \lambda}\right)_{i,j}
                        - \left(\frac{1}{a}\frac{\partial u}{\partial \phi}\right)_{i,j}
                        + \frac{u_{i,j}\tan\phi_{j}}{a}

       緯度、経度微分の項は、:py:func:`dvardx`、 :py:func:`dvardy` を内部で呼ぶ。境界の扱いもこれらに準ずる。
    
    """
    u, v = np.array(u), np.array(v)
    ndim = u.ndim

    out = dvardx(v,lon,lat,xdim,ydim,cyclic=cyclic,sphere=sphere) - dvardy(u,lat,ydim,sphere=sphere)
    if sphere:
       out = out + u*tools.expand(np.tan(lat*d2r),ndim,ydim)/a0
    
    out = np.rollaxis(out, ydim, 0)
    out[0,...]  = 0.
    out[-1,...] = 0.
    out = np.rollaxis(out, 0, ydim+1)
    
    return out

def laplacian(var, lon, lat, xdim, ydim, cyclic=True, sphere=True):
    ur"""
    球面上での2次元ラプラシアンを計算する。

    :Arguments:
     **var** : ndarray
       スカラー場
     **lon, lat** : array_like
       緯度と経度
     **xdim, ydim** : int
       緯度、経度の軸
     **cyclic** : bool, optional
       経度微分の際に周境界とするかどうか。デフォルトはTrue
     **sphere** : bool, optional
       球面緯度経度座標かどうか。デフォルトはTrue。Falseにすると直交座標として扱う。
    
    :Returns:
     **out** : ndarray
       varと同じ形状のndarray
       
    .. note::
       球面緯度経度座標系における2次元ラプラシアンは次のように定義される。

       .. math:: \nabla_{h}^2\Phi
                    &= \frac{1}{a^2\cos^2\phi} \left[ \frac{\partial^2 \Phi}{\partial \lambda^2} + \cos\phi\frac{\partial}{\partial \phi}
                       \left( \cos\phi \frac{\partial \Phi}{\partial \phi} \right) \right] \\
                    &= \frac{1}{a^2\cos^2\phi}\frac{\partial \Phi}{\partial \lambda} + \frac{1}{a^2}\frac{\partial^2 \Phi}{\partial \phi^2}
                       - \frac{\tan\phi}{a^2}\frac{\partial \Phi}{\partial \phi}

       ここで、aは地球半径。これを中央差分で次のように計算する。
       
       .. math:: (\nabla_{h}^2\Phi)_{i,j}
                    = \left(\frac{1}{a^2\cos^2\phi}\frac{\partial^2 \Phi}{\partial \lambda^2}\right)_{i,j}
                        - \left(\frac{1}{a^2}\frac{\partial^2 \Phi}{\partial \phi^2}\right)_{i,j}
                        - \frac{\tan\phi_{j}}{a} \left( \frac{1}{a}\frac{\partial \Phi}{\partial \phi}\right)_{i,j}

       緯度、経度微分の項は、:py:func:`d2vardx2`、 :py:func:`d2vardy2` を内部で呼ぶ。境界の扱いもこれらに準ずる。        
    """
    var = np.asarray(var)
    ndim = var.ndim

    if sphere:
        out = d2vardx2(var, lon, lat, xdim, ydim, cyclic=cyclic, sphere=sphere) + d2vardy2(var, lat, ydim, sphere=sphere) - tools.expand(np.tan(lat*d2r),ndim,ydim)*dvardy(var, lat, ydim)/a0
    else:
        out = d2vardx2(var, lon, lat, xdim, ydim, cyclic=cyclic, sphere=sphere) + d2vardy2(var, lat, ydim, sphere=sphere)
        
    return out

def grad(var, lon, lat, xdim, ydim, cyclic=True, sphere=True):
    ur"""
    水平勾配を計算する
    :Arguments:
     **var** : ndarray
       スカラー場
     **lon, lat** : array_like
       緯度と経度、もしくはx座標とy座標
     **xdim, ydim** : int
       緯度、経度の軸
     **cyclic** : bool, optional
       経度微分の際に周境界とするかどうか。デフォルトはTrue
     **sphere** : bool, optional
       球面緯度経度座標かどうか。デフォルトはTrue。Falseにすると直交座標として扱う。
    
    :Returns:
     **outu, outv** : ndarray
       varと同じ形状のndarray

    .. note::
    
       球面緯度経度座標系における勾配は、次のように定義される。

       .. math:: \mbox{\boldmath $\nabla$}\Phi = \frac{\mbox{\boldmath $i$}}{a\cos\phi}\frac{\partial \Phi}{\partial \lambda}
                                                  + \frac{\mbox{\boldmath $j$}}{a}\frac{\partial \Phi}{\partial \phi}                                               

       ここで、aは地球半径。
       
    """
    var = np.asarray(var)

    outu = dvardx(var,lon,lat,xdim,ydim,cyclic=True, sphere=sphere)
    outv = dvardy(var,lat,ydim, sphere=sphere)    
    
    return outu, outv

def skgrad(var, lon, lat, xdim, ydim, cyclic=True, sphere=True):
    ur"""
    skew gradient（流線関数からベクトルを求める）を計算する
    :Arguments:
     **var** : ndarray
       スカラー場
     **lon, lat** : array_like
       緯度と経度、もしくはx座標とy座標
     **xdim, ydim** : int
       緯度、経度の軸
     **cyclic** : bool, optional
       経度微分の際に周境界とするかどうか。デフォルトはTrue
     **sphere** : bool, optional
       球面緯度経度座標かどうか。デフォルトはTrue。Falseにすると直交座標として扱う。
    
    :Returns:
     **outu, outv** : ndarray
       varと同じ形状のndarray

    .. note::
    
       球面緯度経度座標系におけるskew gradientは、次のように定義される。

       .. math:: \mbox{\boldmath $k$}\cdot\mbox{\boldmath $\nabla$}\Phi = -\frac{\mbox{\boldmath $i$}}{a}\frac{\partial \Phi}{\partial \phi}
                                                  + \frac{\mbox{\boldmath $j$}}{a\cos\phi}\frac{\partial \Phi}{\partial \lambda}                                               


       ここで、aは地球半径。
       
    """
    var = np.asarray(var)
 
    outu = -dvardy(var,lat,ydim, sphere=sphere)       
    outv = dvardx(var,lon,lat,xdim,ydim,cyclic=True, sphere=sphere)
    
    return outu, outv

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

def vint(var, bottom, top, lev, zdim, punit=100.):
    ur"""
    質量重み付き鉛直積分。
    
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

    .. note::
       p面座標系で与えられるデータの鉛直積分は次式で定義される。
       
       .. math:: vint = \frac{1}{g}\int_{bottom}^{top} \Phi dp
       
    **Examples**

    >>>
    >>>
    """
    var = np.ma.asarray(var)
    lev = np.asarray(lev)
    ndim = var.ndim

    lev = lev[(lev <= bottom)&(lev >= top)]
    lev_m = np.r_[bottom,(lev[1:] + lev[:-1])/2.,top]
    dp = lev_m[:-1] - lev_m[1:]

    #roll lat dim axis to last
    var = tools.mrollaxis(var,zdim,ndim)
    out = var[...,(lev <= bottom)&(lev >= top)] * dp / g * punit
    if bottom > top:
        out = out.sum(axis=-1)
    else:
        out = -out.sum(axis=-1)
    return out

def vmean(var, bottom, top, lev, zdim, punit=100.):
    ur"""
    質量重み付き鉛直平均。
    
    :Arguments:
     **var** : array_like
       データ。       
     **bottom, top** : float
       平均の下端、上端。      
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
    var = np.ma.asarray(var)
    lev = np.asarray(lev)
    ndim = var.ndim

    lev = lev[(lev <= bottom)&(lev >= top)]
    lev_m = np.r_[bottom,(lev[1:] + lev[:-1])/2.,top]
    dp = lev_m[:-1] - lev_m[1:]

    #roll lat dim axis to last
    var = tools.mrollaxis(var,zdim,ndim)
    out = var[...,(lev <= bottom)&(lev >= top)] * dp 
    out = out.sum(axis=-1)/(dp.sum())

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
    nx1, ny1 = np.size(lon1), np.size(lat1)
    nx2, ny2 = np.size(lon2), np.size(lat2)

    if nx1 == 1 and nx2 == 1 :
        lon2 = np.asarray(lon2)
        lat2 = np.asarray(lat2)
    elif nx2 == 1 and ny2 == 1 :
        lon1 = np.asarray(lon1)
        lat1 = np.asarray(lat1)
    elif nx1 > 1 and ny1 > 1 and nx2 > 1 and ny2 > 1:
        lon1 = np.asarray(lon1)
        lat1 = np.asarray(lat1)
        lon2 = np.asarray(lon2)
        lat2 = np.asarray(lat2)

    return a0 * np.arccos(np.sin(lat1*d2r)*np.sin(lat2*d2r)
                        +np.cos(lat1*d2r)*np.cos(lat2*d2r)*np.cos((lon2-lon1)*d2r))


