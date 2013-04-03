# coding: utf-8
#--------------------------------------------------------
# -- dynamicsへのラッパー
#--------------------------------------------------------
import pymet.dynamics as dynamics
from core import *

__all__ = ['pottemp', 'ertelpv', 'tnflux2d']

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

def ertelpv(ufield, vfield, tfield, cyclic=True, punit=100.):
    if not isinstance(ufield, McField) or not isinstance(vfield, McField) or not isinstance(tfield, McField):
        raise TypeError, "input must be McField instance"
    grid = ufield.grid.copy()
    uwnd = np.ma.getdata(ufield, subok=False)
    vwnd = np.ma.getdata(vfield, subok=False)
    temp = np.ma.getdata(tfield, subok=False)
    mask = np.ma.getmask(ufield) | np.ma.getmask(vfield) | np.ma.getmask(tfield) 

    result = dynamics.ertelpv(uwnd, vwnd, temp, grid.lon, grid.lat, grid.lev,
                              grid.xdim, grid.ydim, grid.zdim,
                              cyclic=cyclic, punit=grid.punit)

    return McField(result, name='ertelpv', grid=grid, mask=mask)

def tnflux2d(Ufield, Vfield, strmfield, cyclic=True, limit=100.):
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
    tny = McField(result, name='tnflux_y', grid=grid.copy(), mask=mask)
    return tnx, tny

