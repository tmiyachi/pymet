# coding:utf-8
"""
"""
from grads import GaCore
from pymet.field import McField, McGrid
import numpy as np
from datetime import datetime
import os, os.path
import locale

__all__ = ['GradsIO']

class GradsIO:
    u"""
    pygradsを用いてMcFieldデータを読み込むためのクラス。

    :Parameters:
      **Echo** : bool, optional
         gradsの出力を標準出力に表示するかどうか。デフォルトはFalse。

    .. note::


    **Examples**
     >>> gaio = pymet.io.GradsIO()

    **Attrubutes**

    ======= ====================================
    ga
    ======= ====================================

    **Methods**
    .. currentmodule:: pymet.io.gradsio.GradsIO
    
    ..  autosummary::

        open
        command
        setdim
        get
        __init__

    """
    def __init__(self, Echo=False):
        locale.setlocale(locale.LC_ALL,'en_US')
        ga = GaCore('grads -b', Echo=Echo)
        self.ga = ga
    def open(self, fname, Quiet=True):
        if not os.path.isfile(fname):
            raise IOError, "cannot find {0}".format(fname)
        root, ext = os.path.splitext(fname)
        print ext, fname
        if ext.upper() == '.NC' or ext.upper() == '.NETCDF':
            ftype='SDF'
        else:
            ftype='default'
        self.ga.open(fname, ftype=ftype, Quiet=Quiet)
        coords = self.ga.coords()
        if coords.lon[-1] - coords.lon[0] == 360.:            
            self.ga('set lon %f %f' % (coords.lon[0], coords.lon[-2]))

    def command(self,command_string):
        self.ga(command_string)
    def setdim(self, time=None, lon=None, lat=None, lev=None, ens=None):
        ga = self.ga
        if hasattr(lon, '__iter__'):
            ga('set lon %f %f' % lon)
        elif lon==None:
            coords = ga.coords()
            if coords.lon[-1] - coords.lon[0] == 360.:            
                ga('set lon %f %f' % (coords.lon[0], coords.lon[-2]))
        else:
            ga('set lon %f' % lon)
        if hasattr(lat, '__iter__'):
            ga('set lat %f %f' % lat)
        elif lat!=None:
            ga('set lat %f' % lat)
        if hasattr(lev, '__iter__'):
            ga('set lev %f %f' % lev)
        elif lev!=None:
            ga('set lev %f' % lev)
        if hasattr(time, '__iter__'):
            ga('set time %s %s' % (time[0].strftime('%Hz%d%b%Y'), time[1].strftime('%Hz%d%b%Y')))
        elif time != None:
            ga('set time %s' % time.strftime('%Hz%d%b%Y'))
        if hasattr(ens, '__iter__'):
            ga('set e %d %d' % ens)
        elif type(ens) == type(int(1)):
            ga('set e %d' % ens)
        elif ens == None:
            ga('set e 1 %d' % ga.query('file',Quiet=True).ne)
            
    def get(self, var):
        """
        """
        ga = self.ga
        info = ga.coords()
        out  = np.asarray(ga.eval(var), dtype=np.float32)
        out  = np.ma.array(out, mask=(out==info.undef))

        out  = out.reshape(info.shape)

        lon  = np.asarray(info.lon)
        lat  = np.asarray(info.lat)
        lev  = np.asarray(info.lev)
        tyme = [datetime.strptime(time, '%HZ%d%b%Y') for time in info.time]
        dims = info.dims
        print lon
        # create McGrid
        grid = McGrid(name=var, lon=lon, lat=lat, lev=lev, tyme=tyme, dims=dims)
        field = McField(out, name=var, grid=grid, mask=out.mask)
        return field

        
