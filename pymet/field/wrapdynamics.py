# coding: utf-8
#--------------------------------------------------------
# -- dynamicsへのラッパー
#--------------------------------------------------------
import pymet.dynamics as dynamics
import numpy as np
from core import *

__all__ = ['pottemp', 'ertelpv', 'stability', 'tnflux2d', 'tnflux3d', 'absvrt', 'rhmd']

def pottemp(tfield, p0=1000000.):
    u"""
    温位を気温から計算する。

    :Arguments: 
     **tfield** : McField object
      気温 [K]
     **p0** : float, optional
      標準気圧。デフォルトは1000hPa。

    .. seealso::
    
       .. autosummary::
          :nosignatures:

          pymet.dynamics.pottemp          
     
    """
    if not isinstance(tfield, McField):
        raise TypeError, "input must be McField instance"
    grid = tfield.grid.copy()
    t = np.ma.getdata(tfield, subok=False)
    mask = np.ma.getmask(tfield) 

    result = dynamics.pottemp(t, grid.lev, grid.zdim, punit=grid.punit, p0=p0)

    return McField(result, name='theta', grid=grid, mask=mask)

def ertelpv(ufield, vfield, tfield, cyclic=True):
    u"""
    エルテルのポテンシャル渦度を計算する。

    :Arguments:
     **ufield, vfield** : McField object
      東西風、南北風 [m/s]
     **tfield**       : McField object
      気温 [K]
     **cyclic** : bool, optional
      東西境界を周期境界で扱うか。デフォルトはTrue
      
    :Returns:
      **out** : McField object

    .. seealso::
    
       .. autosummary::
          :nosignatures:

          pymet.dynamics.ertelpv          
          
    """
    if not isinstance(ufield, McField) or not isinstance(vfield, McField) or not isinstance(tfield, McField):
        raise TypeError, "input must be McField instance"
    grid = ufield.grid.copy()
    uwnd = np.ma.getdata(ufield, subok=False)
    vwnd = np.ma.getdata(vfield, subok=False)
    temp = np.ma.getdata(tfield, subok=False)
    mask = np.ma.getmask(ufield) | np.ma.getmask(vfield) | np.ma.getmask(tfield) 

    result = dynamics.ertelpv(uwnd, vwnd, temp, grid.lon, grid.lat, grid.lev,
                              grid.xdim, grid.ydim, grid.zdim,
                              cyclic=cyclic, punit=grid.punit, sphere=grid.sphere)

    return McField(result, name='ertelpv', grid=grid, mask=mask)

def stability(tfield):
    u"""
    p座標系での静的安定度(Brunt-Vaisala振動数)を計算する。

    :Arguments:
     **tfield** : McField object
      気温 [K]
      
    :Returns:
      **out** : McField object

    .. seealso::
    
       .. autosummary::
          :nosignatures:

          pymet.dynamics.stability
          
    """
    if not isinstance(tfield, McField):
        raise TypeError, "input must be McField instance"
    grid = ufield.grid.copy()
    temp = np.ma.getdata(tfield, subok=False)
    mask = np.ma.getmask(tfield)

    result = dynamics.stability(temp, grid.lev, grid.zdim, punit=grid.punit)

    return McField(result, name='stability', grid=grid, mask=mask)

def tnflux2d(Ufield, Vfield, strmfield, cyclic=True, limit=100.):
    u"""
    Takaya & Nakamura (2001) の波活動度フラックスの水平成分をp面上で計算する。

    :Arguments:
     **Ufield, Vfield** : McField object
      気候値東西風、南北風 [m/s]
     **strmfield** : McField object
      流線関数偏差 [m^2/s]
     **cyclic** : bool, optional
      東西境界を周期境界で扱うか。デフォルトはTrue。
     **limit** : float, optional
      デフォルトは100.
     
    :Returns:
     **tnx, tny** : McField object
      フラックスの東西、南北成分

    .. seealso::
    
       .. autosummary::
          :nosignatures:

          pymet.dynamics.tnflux2d    
     """     

    if not isinstance(Ufield, McField) or not isinstance(Vfield, McField) or not isinstance(strmfield, McField):
        raise TypeError, "input must be McField instance"
    grid = Ufield.grid.copy()
    U = np.ma.getdata(Ufield, subok=False)
    V = np.ma.getdata(Vfield, subok=False)
    strm = np.ma.getdata(strmfield, subok=False)
    mask = np.ma.getmask(Ufield) | np.ma.getmask(Vfield) | np.ma.getmask(strmfield) 

    tnx, tny = dynamics.tnflux2d(U, V, strm, grid.lon, grid.lat, grid.xdim, grid.ydim,
                               cyclic=cyclic, limit=limit)
    tnx = McField(tnx, name='tnflux_x', grid=grid.copy(), mask=mask)
    tny = McField(tny, name='tnflux_y', grid=grid.copy(), mask=mask)
    return tnx, tny

def tnflux3d(Ufield, Vfield, Tfield, strmfield, cyclic=True, limit=100.):
    u"""
    Takaya & Nakamura (2001) の波活動度フラックスをp面上で計算する。

    :Arguments:
     **Ufield, Vfield** : McField object
      気候値東西風、南北風 [m/s]
     **Tfield**    : McField object
      気候値気温 [K]
     **strmfield** : McField object
      流線関数偏差 [m^2/s]
     **cyclic** : bool, optional
      東西境界を周期境界で扱うか。デフォルトはTrue。
     **limit** : float, optional
      デフォルトは100.
     
    :Returns:
     **tnx, tny, tnz** : McField object
      フラックスの東西、南北、鉛直成分

    .. seealso::
    
       .. autosummary::
          :nosignatures:

          pymet.dynamics.tnflux3d    

    """
    if not isinstance(Ufield, McField) or not isinstance(Vfield, McField) or not isinstance(strmfield, McField) or not isinstance(Tfield, McField):
        raise TypeError, "input must be McField instance"
    grid = Ufield.grid.copy()
    U = np.ma.getdata(Ufield, subok=False)
    V = np.ma.getdata(Vfield, subok=False)
    T = np.ma.getdata(Tfield, subok=False)
    strm = np.ma.getdata(strmfield, subok=False)
    mask = np.ma.getmask(Ufield) | np.ma.getmask(Vfield) | np.ma.getmask(Tfield) | np.ma.getmask(strmfield) 

    tnx, tny, tnz = dynamics.tnflux3d(U, V, T, strm, grid.lon, grid.lat, grid.lev, grid.xdim, grid.ydim,
                                      grid.zdim, cyclic=cyclic, limit=limit, punit=grid.punit)
    tnx = McField(tnx, name='tnflux_x', grid=grid.copy(), mask=mask)
    tny = McField(tny, name='tnflux_y', grid=grid.copy(), mask=mask)
    tnz = McField(tnz, name='tnflux_z', grid=grid.copy(), mask=mask)
    return tnx, tny, tnz

def absvrt(ufield, vfield, cyclic=True):
    u"""
    """
    if not isinstance(ufield, McField) or not isinstance(vfield, McField):
        raise TypeError, "input must be McField instance"
    grid = ufield.grid.copy()
    uwnd = np.ma.getdata(ufield, subok=False)
    vwnd = np.ma.getdata(vfield, subok=False)
    mask = np.ma.getmask(ufield) | np.ma.getmask(vfield)

    result = dynamics.absvrt(uwnd, vwnd, grid.lon, grid.lat,
                             grid.xdim, grid.ydim, cyclic=cyclic, sphere=grid.sphere)

    return McField(result, name='absvrt', grid=grid, mask=mask)

def rhmd(tfield, qvalfield, qtyp='q'):
    u"""
    """
    if not isinstance(tfield, McField) or not isinstance(qvalfield, McField):
        raise TypeError, "input must be McField instance"
    grid = tfield.grid.copy()
    temp = np.ma.getdata(tfield, subok=False)
    qval = np.ma.getdata(qvalfield, subok=False)
    mask = np.ma.getmask(tfield) | np.ma.getmask(qvalfield)

    result = dynamics.rhmd(temp, qval, grid.lev, grid.zdim, punit=grid.punit, qtyp=qtyp)

    return McField(result, name='rh', grid=grid, mask=mask)

