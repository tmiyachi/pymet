# coding:utf-8
#------------------------------------------------------------------------------        
#---  緯度経度ラベルのためのticker
#------------------------------------------------------------------------------
from mpl_toolkits.basemap import Basemap
import matplotlib.ticker
import matplotlib.dates
import numpy as np
import pymet.tools as tools

__all__ = ['BasemapXaxisLocator','BasemapYaxisLocator',
           'BasemapXaxisFormatter','BasemapYaxisFormatter',
           'XaxisFormatter','YaxisFormatter',
           'DateFormatter', 'ScaleFormatter']


class BasemapXaxisLocator(matplotlib.ticker.FixedLocator):
    u"""
    経度目盛を等間隔でふるためのLocator。
    """
    def __init__(self,baseMap,lonlocs):
        self.baseMap=baseMap
        self.lonlocs = np.asarray(lonlocs)
        self.latlocs = np.ones(len(lonlocs))*baseMap.llcrnrlat
        
    def __call__(self):
        xlocs, ylocs = self.baseMap(self.lonlocs,self.latlocs)
        return xlocs
    
class BasemapYaxisLocator(matplotlib.ticker.FixedLocator):
    u"""
    緯度目盛りを等間隔でふるためのLocator。
    """
    def __init__(self,baseMap,latlocs):
        self.baseMap=baseMap
        self.lonlocs = np.ones(len(latlocs))*baseMap.llcrnrlon
        self.latlocs = np.asarray(latlocs)
    def __call__(self):
        xlocs, ylocs = self.baseMap(self.lonlocs,self.latlocs)
        return ylocs    

class BasemapXaxisFormatter(matplotlib.ticker.Formatter):
    u"""
    Basemapの経度ラベルためのFormatter。
    """
    def __init__(self, baseMap):
        self.baseMap=baseMap
    def __call__(self,x ,pos=None):
        lon,lat = self.baseMap(x,0,inverse=True)
        return tools.lon2txt(lon)
    
class BasemapYaxisFormatter(matplotlib.ticker.Formatter):
    u"""
    Basemapの緯度ラベルのためのFormatter。
    """
    def __init__(self,baseMap):
        self.baseMap=baseMap
    def __call__(self,y,pos=None):
        lon,lat = self.baseMap(0,y,inverse=True)
        return tools.lat2txt(lat)

class XaxisFormatter(matplotlib.ticker.Formatter):
    u"""
    経度ラベルのためのFormatter
    """
    def __call__(self, lon, pos=None):
        return tools.lon2txt(lon)

class YaxisFormatter(matplotlib.ticker.Formatter):
    u"""
    緯度ラベルのためのFormatter
    """
    def __call__(self, lat, pos=None):
        return tools.lat2txt(lat)
    
    
class DateFormatter(matplotlib.ticker.Formatter):
    u"""
    時間軸ラベルのためのFormatter。localeの設定に関係なく%bは英語大文字の略称になる。

    :Arguments:
     **fmt** : str, optional
      デフォルトは、'%HZ%d%b%Y'
    """
    def __init__(self, fmt='%HZ%d%b%Y', tz=None):
        self.fmt = fmt
        self.tz  = tz
    def __call__(self, t, pos=1):
        d = matplotlib.dates.num2date(t, self.tz) 
        return tools.d2s(d, fmt=self.fmt)

class DegreeAxisFormatter(matplotlib.ticker.Formatter):
    def __init__(self, fmt='%g'):
        self.fmt = fmt
    def __call__(self, deg, pos=None):
        fmt = self.fmt
        labstr =  u'%s\N{DEGREE SIGN}'%fmt
        return labstr%deg

class ScaleFormatter(matplotlib.ticker.Formatter):
    u"""
    経度ラベルのためのFormatter
    """
    def __init__(self, scale=1.):
        self.scale = scale
    def __call__(self, d, pos=None):
        return d/scale.
    
    
