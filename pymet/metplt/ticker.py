# coding:utf-8
#------------------------------------------------------------------------------        
#---  緯度経度ラベルのためのticker
#------------------------------------------------------------------------------
from mpl_toolkits.basemap import Basemap
import matplotlib.ticker
import numpy as np

__all__ = ['BasemapXaxisFormatter','BasemapYaxisFormatter',
           'BasemapXaxisLocator','BasemapYaxisLocator',
           'lon2txt','lat2txt']


class BasemapXaxisFormatter(matplotlib.ticker.Formatter):
    u"""
    経度ラベルためのFormatter。
    """
    def __init__(self, baseMap):
        self.baseMap=baseMap
    def __call__(self,x ,pos=1):
        lon,lat = self.baseMap(x,0,inverse=True)
        return lon2txt(lon,pos=pos)
    
class BasemapYaxisFormatter(matplotlib.ticker.Formatter):
    u"""
    緯度ラベルのためのFormatter。
    """
    def __init__(self,baseMap):
        self.baseMap=baseMap
    def __call__(self,y,pos=1):
        lon,lat = self.baseMap(0,y,inverse=True)
        return lat2txt(lat,pos=pos)

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

def lon2txt(lon, pos=None):
    u"""
    経度の値を文字列に変換する。

    0度からの相対経度を東経、西経の文字列(30\N{DEGREE SIGN}E,
    140\N{DEGREE SIGN}Wなど)に変換する。

    :Argumets:
     **lon** : int
      経度(degrees)。0度からの相対経度で表す。
     **pos** : optional
     
    :Returns:
     **lonlab** : str
      経度のラベル。
      
    **Examples**
     >>> lon2txt(135)
     '135\N{DEGREE SIGN}E'
     >>> lon2txt(-30)
     '30\N{DEGREE SIGN}W'
     >>> lon2txt(250)
     '110\N{DEGREE SIGN}W'
    """
    fmt = '%g'
    lon = (lon+360) % 360
    if lon>180:
        lonlabstr = u'%s\N{DEGREE SIGN}W'%fmt
        lonlab = lonlabstr%abs(lon-360)
    elif lon<180 and lon != 0:
        lonlabstr = u'%s\N{DEGREE SIGN}E'%fmt
        lonlab = lonlabstr%lon
    else:
        lonlabstr = u'%s\N{DEGREE SIGN}'%fmt
        lonlab = lonlabstr%lon
    return lonlab

def lat2txt(lat, pos=None):
    u"""
    緯度の値を文字列に変換する。

    緯度を北緯、南緯の文字列(30\N{DEGREE SIGN}N,60\N{DEGREE SIGN}Sなど)
    に変換する。

    :Argumets:
     **lon** : int
      緯度(degrees)。南緯はマイナス、北緯はプラス。
     **pos** : optional
     
    :Returns:
     **lonlab** : str
      緯度のラベル。
      
    **Examples**
     >>> lat2txt(60)
     '60\N{DEGREE SIGN}N'
     >>> lat2txt(-30)
     '30\N{DEGREE SIGN}S'
    """
    fmt = '%g'
    if lat<0:
        latlabstr = u'%s\N{DEGREE SIGN}S'%fmt
        latlab = latlabstr%abs(lat)
    elif lat>0:
        latlabstr = u'%s\N{DEGREE SIGN}N'%fmt
        latlab = latlabstr%lat
    else:
        latlabstr = u'%s\N{DEGREE SIGN}'%fmt
        latlab = latlabstr%lat
    return latlab
