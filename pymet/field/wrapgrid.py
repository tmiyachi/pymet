# coding: utf-8
import pymet.grid
import numpy as np
from core import *
import pymet.tools as tools

__all__ = ['dvardx', 'dvardy', 'dvardp', 'div', 'rot', 'd2vardx2', 'd2vardy2',
           'vint', 'dvardt', 'vmean']

def dvardx(field, cyclic=True)
    u"""
    経度方向の微分を中央差分で計算。

    :Arguments:
     **field** : McField object

     **cyclic** : bool, optional
       東西の境界を周期境界として扱うかどうか。False の場合は周期境界を用いずに
       前方、後方差分で計算する。デフォルトは True。
       
    :Returns:
     **result** : McField object

    .. seealso::

     .. autosummary::
        :nosignatures:
     
        pymet.grid.dvardx
    """
    if not isinstance(field, McField):
        raise TypeError, "field must be McField instance"
    grid = field.grid.copy()
    data = np.ma.getdata(field, subok=False)
    mask = np.ma.getmask(field)

    result = pymet.grid.dvardx(data, grid.lon, grid.lat, grid.xdim, grid.ydim, cyclic=True, sphere=grid.sphere)
    if np.size(result)<2:
        return result
    return McField(result, name=field.name, grid=grid, mask=mask)

def dvardy(field):
    u"""
    緯度方向の微分を中央差分で計算。

    :Arguments:
     **field** : McField object

    :Returns:
     **result** : McField object

    .. seealso::

     .. autosummary::
        :nosignatures:
     
        pymet.grid.dvardy
    """
    if not isinstance(field, McField):
        raise TypeError, "field must be McField instance"
    grid = field.grid.copy()
    data = np.ma.getdata(field, subok=False)
    mask = np.ma.getmask(field)

    result = pymet.grid.dvardy(data, grid.lat, grid.ydim, sphere=grid.sphere)
    if np.size(result)<2:
        return result
    return McField(result, name=field.name, grid=grid, mask=mask)

def dvardp(field):
    u"""
    鉛直方向の微分をlog(p)の中央差分で計算。

    :Arguments:
     **field** : McField object
       
    :Returns:
     **result** : McField object

    .. seealso::

     .. autosummary::
        :nosignatures:
     
        pymet.grid.dvardp
    """
    if not isinstance(field, McField):
        raise TypeError, "field must be McField instance"
    grid = field.grid.copy()
    data = np.ma.getdata(field, subok=False)
    mask = np.ma.getmask(field)

    result = pymet.grid.dvardp(data, grid.lev, grid.zdim, grid.punit)
    if np.size(result) < 2:
        return result
    return McField(result, name=field.name, grid=grid, mask=mask)

def d2vardx2(field, cyclic=True):
    u"""
    経度方向の2階微分を中央差分で計算。

    :Arguments:
     **field** : McField object

     **cyclic** : bool, optional
       東西の境界を周期境界として扱うかどうか。False の場合は周期境界を用いずに
       前方、後方差分で計算する。デフォルトは True。
       
    :Returns:
     **result** : McField object

    .. seealso::

     .. autosummary::
        :nosignatures:
     
        pymet.grid.d2vardx2
    """
    if not isinstance(field, McField):
        raise TypeError, "field must be McField instance"
    grid = field.grid.copy()
    data = np.ma.getdata(field, subok=False)
    mask = np.ma.getmask(field)

    result = pymet.grid.d2vardx2(data, grid.lon, grid.lat, grid.xdim, grid.ydim, cyclic=True)
    if np.size(result) < 2:
        return result
    return McField(result, name=field.name, grid=grid, mask=mask)

def d2vardy2(field):
    u"""
    緯度方向の微分を中央差分で計算。

    :Arguments:
     **field** : McField object
       
    :Returns:
     **result** : McField object

    .. seealso::

     .. autosummary::
        :nosignatures:
     
        pymet.grid.d2vardy2
    """
    if not isinstance(field, McField):
        raise TypeError, "field must be McField instance"
    grid = field.grid.copy()
    data = np.ma.getdata(field, subok=False)
    mask = np.ma.getmask(field)

    result = grid.d2vardy2(data, grid.lat, grid.ydim)
    if np.size(result) < 2:
        return result
    return McField(result, name=field.name, grid=grid, mask=mask)

def div(ufield, vfield, cyclic=True):
    u"""
    水平発散を計算する。

    :Arguments:
     **ufield, vfield** : McField object
      ベクトルの東西、南北成分。

     **cyclic** : bool, optional
       東西の境界を周期境界として扱うかどうか。False の場合は周期境界を用いずに
       前方、後方差分で計算する。デフォルトは True。
       
    :Returns:
     **result** : McField object

    .. seealso::

     .. autosummary::
        :nosignatures:
     
        pymet.grid.div
    """
    if not isinstance(ufield, McField) or not isinstance(vfield, McField):
        raise TypeError, "input must be McField instance"
    grid = ufield.grid.copy()
    u = np.ma.getdata(ufield, subok=False)
    v = np.ma.getdata(vfield, subok=False)
    mask = np.ma.getmask(ufield) | np.ma.getmask(vfield)

    result = pymet.grid.div(u, v, grid.lon, grid.lat, grid.xdim, grid.ydim, cyclic=True, sphere=grid.sphere)
    if np.size(result) < 2:
        return result
    return McField(result, name='div', grid=grid, mask=mask)

def rot(ufield, vfield, cyclic=True):
    u"""
    回転の鉛直成分を計算する。

    :Arguments:
     **ufield, vfield** : McField object
       ベクトルの東西、鉛直成分。
     **cyclic** : bool, optional
       東西の境界を周期境界として扱うかどうか。False の場合は周期境界を用いずに
       前方、後方差分で計算する。デフォルトは True。
       
    :Returns:
     **result** : McField object

    .. seealso::

     .. autosummary::
        :nosignatures:
     
        pymet.grid.rot
    """
    if not isinstance(ufield, McField) or not isinstance(vfield, McField):
        raise TypeError, "input must be McField instance"
    grid = ufield.grid.copy()
    u = np.ma.getdata(ufield, subok=False)
    v = np.ma.getdata(vfield, subok=False)
    mask = np.ma.getmask(ufield) | np.ma.getmask(vfield)

    result = pymet.grid.rot(u, v, grid.lon, grid.lat, grid.xdim, grid.ydim, cyclic=True, sphere=grid.sphere)
    if np.size(result) < 2:
        return result
    return McField(result, name='rot', grid=grid, mask=mask)

def vint(field, bottom, top):
    u"""
    質量重み付き鉛直積分。

    :Argument:
     **field** : McField object
      入力データ。
     **bottom** : float
      下端
     **top** : float
      上端

    :Returns:
     **result** : McField object     
    """
    if not isinstance(field, McField):
        raise TypeError, "input must be McField instance"
    var = np.ma.asarray(field)
    grid = field.grid.copy()
#    var = np.ma.getdata(field, subok=False)
#    mask = np.ma.getmask(field)

    result = pymet.grid.vint(var, bottom, top, lev=grid.lev, zdim=grid.zdim, punit=grid.punit)

    if np.size(result) < 2:
        return result
    grid.name = field.name + '_vint'
    grid.lev = None
    return McField(result, name=grid.name, grid=grid)

def vmean(field, bottom, top):
    u"""
    質量重み付き鉛直平均。

    :Argument:
     **field** : McField object
      入力データ。
     **bottom** : float
      下端
     **top** : float
      上端

    :Returns:
     **result** : McField object     
    """
    if not isinstance(field, McField):
        raise TypeError, "input must be McField instance"
    var = np.ma.asarray(field)
    grid = field.grid.copy()

    result = pymet.grid.vmean(var, bottom, top, lev=grid.lev, zdim=grid.zdim, punit=grid.punit)

    if np.size(result) < 2:
        return result
    grid.name = field.name + '_vmean'
    grid.lev = None
    return McField(result, name=grid.name, grid=grid)

def dvardt(field, bound='mask'):
    grid = field.grid.copy()
    var = np.ma.getdata(field, subok=False)
    mask = np.ma.getmask(field)
    ndim = np.ndim(var)
    tdim = grid.tdim
    
    dvar = np.diff(var,n=2,axis=tdim)
    dt   = map(lambda td: td.total_seconds(), grid.time[2:]-grid.time[:-2])

    if bound == 'mask':
       out = np.ma.empty(var.shape, dtype=field.dtype)
       out = tools.mrollaxis(out, tdim, 0)
       out[0,...]    = np.ma.masked
       out[1:-1,...] = dvar/tools.expand(dt, ndim, tdim)
       out[-1,...]   = np.ma.masked
       out = tools.mrollaxis(out, 0, tdim+1)
    elif bound == 'valid':
        grid.time = grid.time[1:-1]
        out = np.ma.asarray(dvar/dt)
    return McField(out, name=grid.name, grid=grid, mask=out.mask | mask)
       
        
    
    
    
