# coding:utf-8
import os, os.path
import netCDF4
from pymet.field.core import McField, McGrid
import time
import numpy as np
from datetime import datetime

__all__ = ['NetcdfIO', 'NetcdfWrite']

class NetcdfIO(object):
    u"""
    netCDF形式のファイルを読み込むためのクラス。

    **Examples**
     >>> ncio = pymet.io.NetcdfIO()
     >>> ncio.open('air.nc')
     >>> ncio.setdim(lon=(0,180), lat=(-90,90), lev=500, time=datetime(2009,10,13,12))
     >>> air = ncio.get('air')
     >>> ncio.allclose()

    **Attributes**

     ========= =========================================================
     ncs       現在開いているnetCDF4.Datasetクラスのリスト
     fn        現在開いているファイル数
     fnames    現在開いているファイル名
     vars      現在開いているファイルに含まれる変数名
     ========= =========================================================

    **Methods**
     .. autosummary::

        open
        close
        allclose
        setdim
        get
    """
    def __init__(self):
        self.fn = 0
        self.ncs = []
        self.vars = []
        self.fnames = []
        self.dimkwargs = {}
        
    def open(self, fname):
        u"""
        netCDFファイルを開く。

        :Arguments:
         **fname** : str
          netCDFファイルのパス。
        """
        if not os.path.isfile(fname):
            raise IOError, "cannot find file '{0}'".format(fname)
        nc = netCDF4.Dataset(fname, 'r')
        self.fn += 1
        self.ncs.append(nc)        
        self.vars.append(nc.variables.keys())
        self.fnames.append(fname)        
        
    def close(self, fid=None):
        u"""
        netCDFファイルを閉じる。

        :Arguments:
         **fid** : int, optional
          閉じるファイルのid。現在開いているファイルのうち古い順にfid番目を閉じる。指定しない場合は最も新しいファイルを閉じる。
        """
        if fid:
            if fid > self.fn:
                raise ValueError, "file id '{0}' is larger than number of opened files".format(fid)
        else:
            fid = self.fn
        nc = self.ncs.pop(fid-1)
        self.fnames.pop(fid-1)
        self.vars.pop(fid-1)
        nc.close()
        self.fn += -1

    def allclose(self):
        u"""
        全ての開いているファイルを閉じる。
        """
        for nc in self.ncs:
            nc.close()
        self.ncs = []
        self.fnames = []
        self.vars = []
        self.fn = 0
        self.dimkwargs = {}
        
    def setdim(self, **kwargs):
        u"""
        次元の範囲を設定する。

        :Arguments:
         **lon,lat,lev,time,ens** :
          経度、緯度、鉛直、時間、アンサンブル次元の範囲。
          指定の仕方は :py:func:`McGrid.gridmask` に準ずる。
        """
        self.dimkwargs.update(kwargs)

    def get(self, var, fid=None, **kwargs):
        u"""
        値をMcFieldとして取得する。

        :Arguments:
         **var** : str
          変数名。
         **fid** : int, optional
          ファイルid。変数を取得するファイルidを指定する。指定しない場合は、変数名を含む最も新しく開いたファイルから取得する。
         **lon,lat,lev,time,ens** :
          経度、緯度、鉛直、時間、アンサンブル次元の範囲。
          指定の仕方は :py:func:`McGrid.gridmask` に準ずる。
          :py:func:`setdim`で指定した範囲がある場合は上書きされる。
        """
        if fid:
            if fid > self.fn:
                raise ValueError, "file id '{0}' is larger than number of opened files".format(fid)
            if self.vars[fid-1].count(var) < 1:
                raise ValueError, "Cannot find variable {0} in opened files".format(gavar)
        else:
            var_in_fids = [var in fvars for fvars in self.vars]            
            if var_in_fids.count(True) == 0:
                raise ValueError, "Cannot find variable '{0}' in opened files".format(var)
            elif var_in_fids.count(True) > 1: 
                print "Warnig, multiple variables '{0}' are fond in opened files".format(var)
            fid = var_in_fids.index(True,-1) + 1

        nc = self.ncs[fid-1]

        ncvar = nc.variables[var]
        ncdims = ncvar.dimensions

        grid = McGrid()
        for i, dim in enumerate(ncdims):
            dimval = nc.variables[dim]
            if dimval.axis == 'X':
                grid.lon = dimval[:]
            elif dimval.axis == 'Y':
                grid.lat = dimval[:]
            elif dimval.axis == 'Z':
                grid.lev = dimval[:]
            elif dimval.axis == 'T':
                calendar = getattr(dimval, 'calendar', 'standard')
                times = netCDF4.num2date(dimval[:],units=dimval.units,calendar=calendar)
                grid.time = times
            elif dimval.axis == 'E':
                grid.ens = dimval[:]
                        
        # 指定した範囲のマスク配列をつくる
        dimkwargs = self.dimkwargs.copy()
        dimkwargs.update(kwargs)
        dimkwargs = dict((k,v) for k,v in dimkwargs.items() if k in grid.dims)
        gridmask = grid.gridmask(**dimkwargs)
        
        # 範囲の値を取得する
        grid = grid.getgrid(**dimkwargs)        
        data = np.ma.asarray(ncvar[:])[gridmask].squeeze()

        return McField(data, name=var, grid=grid, mask=data.mask)

class NetcdfWrite(object):
    u"""
    netCDF ver.4 形式のファイルを作成するためのクラス。

    :Argument:
     **fname** : str
      ファイル名のパス。
     **makedirs** : bool, optional
      ファイル名のパスに存在ないディレクトリが含まれるときに作成する場合はTrue。
      デフォルトはFalseで、存在しない場合は例外を発生させる。
    """
    def __init__(self, fname, makedirs=False):
        self.open(fname, makedirs)
        
    def __del__(self):
        if hasattr(self, 'nc') and isinstance(self.nc, netCDF4.Dataset):
            nc.close()
        
    def open(self, fname, makedirs=False):
        u"""
        作成するnetcdfファイルを開き、Datasetオブジェクトを作成する。
        
        :Argument:
         **fname** : str
          ファイル名のパス。
         **makedirs** : bool, optional
          ファイル名のパスに存在ないディレクトリが含まれるときに作成する場合はTrue。
          デフォルトはFalseで、存在しない場合は例外を発生させる。
        """
        if hasattr(self,'nc'):
            raise IOError, "netCDF4.Dataset instance is still created"
        dirname = os.path.dirname(fname)
        if not os.path.isdir(dirname) and len(dirname)>0:
            if makedirs:
                os.makedirs(dirname)
            else:
                raise IOError, "'{0}' does not exsit".format(dirname)
        nc = netCDF4.Dataset(fname, 'w', fprmat='NETCDF4')
        self.nc = nc
    def close(self):
        u"""
        開いているDatasetオブジェクトを閉じる。
        """
        self.nc.close()
        self.nc = None
        
    def setglobalattr(self, **kwargs):
        u"""
        グローバル属性を設定する。        
        """
        nc = self.nc
        history = kwargs.get('history', None)
        if history:
            history += '\n'
        else:
            history = ''
        history += '{0} created by pymet.io.netcdfio'.format(time.ctime())
        kwargs['history'] = history
        for key, val in kwargs.items():
            nc.setncattr(key, val)
      
    def createdim(self, dimname, dimval):
        u"""
        次元を作成し、次元変数を設定する。

        :Arguments:
         **dimname** : str
          次元名
         **dimval** : array_like
          次元の値のリスト

        .. note::
         lon,lat,lev,time,ensは基本次元としてデフォルト値を設定する。
        """
        nc = self.nc
        if dimname in nc.dimensions.keys():
            raise ValueError, "dimension '{0}' still exists".format(dimname)
        dim = nc.createDimension(dimname, size=len(dimval))

        if dimname == 'lon':
            longitudes = nc.createVariable('lon','f4',('lon',))
            longitudes.setncattr('long_name', 'Longitude')
            longitudes.setncattr('units', 'degrees_east')
            longitudes.setncattr('standart_name', 'longitude')
            longitudes.setncattr('axis', 'X')
            longitudes[:] = np.float32(dimval)     
        elif dimname == 'lat':
            latitudes = nc.createVariable('lat','f4',('lat',))            
            latitudes.setncattr('long_name', 'Latitude')
            if dimval[0] > dimval[-1]:
                latitudes.setncattr('units', 'degrees_north')
            else:
                latitudes.setncattr('units', 'degrees_south')                
            latitudes.setncattr('standart_name', 'latitude')
            latitudes.setncattr('axis', 'Y')
            latitudes[:] = np.float32(dimval)                        
        elif dimname == 'lev':
            levels = nc.createVariable('lev','f4',('lev',))
            levels.setncattr('long_name', 'Level')
            levels.setncattr('units', 'millibar')
            levels.setncattr('axis', 'Z')
            levels.setncattr('positive', 'down')
            levels[:] = np.float32(dimval)            
        elif dimname == 'time':
            times = nc.createVariable('time','f8',('time',))
            times.setncattr('units','hours since 0001-01-01 00:00:00.0')
            times.setncattr('calendar', 'standard')
            times.setncattr('long_name', 'Time')
            times.setncattr('axis', 'T')
            times[:] = np.float64(netCDF4.date2num(dimval, units=times.units,
                                  calendar=times.calendar))
        elif dimname == 'ens':
            ensembles = nc.createVariable('ens','i4',('ens',))
            ensembles.setncattr('long_name', 'Ensemble member')
            ensembles.setncattr('axis', 'E') 
            
            
    def createvar(self, field, dims=None, **kwargs):
        if isinstance(field, McField):
            dims = dims or field.grid.dims            
        if dims == None:
            raise ValueError, "dims is required if input field is not McField object"
        
        nc = self.nc
        grid = field.grid
        data = np.ma.asarray(field)
        
        name = kwargs.get('name', field.name)
        long_name = kwargs.get('long_name', field.name)
        units = kwargs.get('units', '')
        dtype = kwargs.get('dtype', 'i2')

        for dimname in dims:
            dimval = getattr(grid,dimname)
            try:
                self.createdim(dimname, dimval)
            except ValueError:
                continue
            
        var = nc.createVariable(name, dtype, tuple(grid.dims))
        if dtype == 'i2':
            ## 最大値、最小値が2バイト整数の範囲に入るようにスケール
            maxval = data.max()
            minval = data.min()
            add_offset = (maxval+minval)/2.
            scale_factor = (maxval-minval)/(32767-(-32767))
            missing_value = np.int16(-32768)
            var.setncattr('add_offset', add_offset)
            var.setncattr('scale_factor', scale_factor)
            var.setncattr('missing_value', missing_value)
            var.set_auto_maskandscale(True)
        else:
            missing_value = data.fill_value
            var.set_auto_maskandscale(True)
            
        var.setncattr('long_name', long_name)
        var.setncattr('units', units)
        var[:] = data
            
