# coding: utf-8
import pymet.grid as grid
from core import *

__all__ = ['dvardx', 'dvardy', 'dvardp', 'div', 'rot', 'd2vardx2', 'd2vardy2']

def dvardx(field, cyclic=True):
    u"""
    経度方向の微分を中央差分で計算。

    :Arguments:
     **field** : McField object

     **cyclic** : bool, optional
       東西の境界を周期境界として扱うかどうか。False の場合は周期境界を用いずに
       前方、後方差分で計算する。デフォルトは True。
       
    :Returns:
     **result** : McField object

    **See Also**

    """
    if not isinstance(field, McField):
        raise TypeError, "field must be McField instance"
    grid = field.grid.copy()
    data = np.ma.getdata(field, subok=False)
    mask = np.ma.getmask(field)

    result = grid.dvardx(data, grid.lon, grid.xdim, cyclic=True)
    if not result.ndim:
        return result
    return McField(result, name=field.name, grid=grid, mask=mask)

def dvardy(field):
    if not isinstance(field, McField):
        raise TypeError, "field must be McField instance"
    grid = field.grid.copy()
    data = np.ma.getdata(field, subok=False)
    mask = np.ma.getmask(field)

    result = grid.dvardy(data, grid.lat, grid.ydim)
    if not result.ndim:
        return result
    return McField(result, name=field.name, grid=grid, mask=mask)

def dvardp(field):
    if not isinstance(field, McField):
        raise TypeError, "field must be McField instance"
    grid = field.grid.copy()
    data = np.ma.getdata(field, subok=False)
    mask = np.ma.getmask(field)

    result = grid.dvardp(data, grid.lev, grid.zdim, grid.punit)
    if not result.ndim:
        return result
    return McField(result, name=field.name, grid=grid, mask=mask)

def d2vardx2(field, cyclic=True):
    if not isinstance(field, McField):
        raise TypeError, "field must be McField instance"
    grid = field.grid.copy()
    data = np.ma.getdata(field, subok=False)
    mask = np.ma.getmask(field)

    result = grid.d2vardx2(data, grid.lon, grid.xdim, cyclic=True)
    if not result.ndim:
        return result
    return McField(result, name=field.name, grid=grid, mask=mask)

def d2vardy2(field):
    if not isinstance(field, McField):
        raise TypeError, "field must be McField instance"
    grid = field.grid.copy()
    data = np.ma.getdata(field, subok=False)
    mask = np.ma.getmask(field)

    result = grid.d2vardy2(data, grid.lat, grid.ydim)
    if not result.ndim:
        return result
    return McField(result, name=field.name, grid=grid, mask=mask)

def div(ufield, vfield, cyclic=True):
    if not isinstance(ufield, McField) or not isinstance(vfield, McField):
        raise TypeError, "input must be McField instance"
    grid = ufield.grid.copy()
    u = np.ma.getdata(ufield, subok=False)
    v = np.ma.getdata(vfield, subok=False)
    mask = np.ma.getmask(ufield) | np.ma.getmask(vfield)

    result = grid.div(u, v, grid.lon, grid.lat, grid.xdim, grid.ydim, cyclic=True)
    if not result.ndim:
        return result
    return McField(result, name='div', grid=grid, mask=mask)

def rot(ufield, vfield, cyclic=True):
    if not isinstance(ufield, McField) or not isinstance(vfield, McField):
        raise TypeError, "input must be McField instance"
    grid = ufield.grid.copy()
    u = np.ma.getdata(ufield, subok=False)
    v = np.ma.getdata(vfield, subok=False)
    mask = np.ma.getmask(ufield) | np.ma.getmask(vfield)

    result = grid.rot(u, v, grid.lon, grid.lat, grid.xdim, grid.ydim, cyclic=True)
    if not result.ndim:
        return result
    return McField(result, name='rot', grid=grid, mask=mask)
