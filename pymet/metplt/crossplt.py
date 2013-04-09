# coding:utf-8
#------------------------------------------------------------------------------
# 鉛直断面図のための関数
#------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

import ticker

__all__ = ['crossplot','crosscontour','crosscontourf']

def crosscontourf(xy,z,data,zlog=True,zinvert=True,xlabel=None,**kwargs):
    u"""
    東西・南北-鉛直断面の塗りつぶしコンタープロット。

    :Arguments:
     **xy** : array_like
      東西、南北座標または緯度・経度。
     **z** : array_like
      鉛直座標。
     **data** : 2darray
      プロットするデータ。
     **zlog** : bool, optional
      縦軸を対数座標にするかどうか。デフォルトはTrue。
     **zinvert** : bool, optional
      縦軸の向きを逆(減少方向)にするかどうか。デフォルトはTrue。
     **xlabel** : {'lon', 'lat'}, optional
      横軸を経度表示にする場合はlon,緯度表示にする場合はlatを指定する。

    :Returns:
     **CF**

    **Examples**
    
    """
    ax = kwargs.get('ax',plt.gca())
    if zinvert:
        ax.set_ylim(z.max(),z.min())
    if zlog: 
        ax.set_yscale('log')
        subs = [1,2,3,4,5,6,7,8,9]
        loc = matplotlib.ticker.LogLocator(base=10.,subs=subs)
        fmt = matplotlib.ticker.FormatStrFormatter("%g")
        ax.yaxis.set_major_locator(loc)
        ax.yaxis.set_major_formatter(fmt)
    if xlabel=='lon':
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticker.lon2txt))
    elif xlabel=='lat':
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticker.lat2txt))
    return ax.contourf(xy,z,data,**kwargs)

def crosscontour(xy,z,data,zlog=True,zinvert=True,xlabel=None,cint=None,**kwargs):
    u"""
    東西・南北-鉛直断面のコンタープロット。

    :Arguments:
     **xy** : array_like
      東西、南北座標または緯度・経度。
     **z** : array_like
      鉛直座標。
     **data** : 2darray
      プロットするデータ。
     **zlog** : bool, optional
      縦軸を対数座標にするかどうか。デフォルトはTrue。
     **zinvert** : bool, optional
      縦軸の向きを逆(減少方向)にするかどうか。デフォルトはTrue。
     **xlabel** : {'lon', 'lat'}, optional
      横軸を経度表示にする場合はlon,緯度表示にする場合はlatを指定する
     **cint** : float, optional
      コンター間隔。

    :Returns:
     **CR**

    **Examples**
    
    """
    ax = kwargs.get('ax',plt.gca())
    kwargs.setdefault('colors','k')
    if zinvert:
        ax.set_ylim(z.max(),z.min())
    if zlog: 
        ax.set_yscale('log')
        z = np.asarray(z)
        subs = [1,2,3,4,5,6,7,8,9]
        loc = matplotlib.ticker.LogLocator(base=10.,subs=subs)
        fmt = matplotlib.ticker.FormatStrFormatter("%g")
        ax.yaxis.set_major_locator(loc)
        ax.yaxis.set_major_formatter(fmt)
    if xlabel=='lon':
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticker.lon2txt))
    elif xlabel=='lat':
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticker.lat2txt))
    if cint != None:
        kwargs.setdefault('locator',matplotlib.ticker.MultipleLocator(cint))
        ax.text(1,-0.1,'contour interval = %g' % cint, 
                transform=ax.transAxes,ha='right')#,fontsize=10)
    return ax.contour(xy,z,data,**kwargs)

def crossplot(xy,z,zlog=True,zinvert=True,**kwargs):
    u"""
    鉛直断面プロット。

    :Arguments:
     **xy** : array_like
      x or y座標。
     **z** : array_like
      鉛直座標。
     **zlog** : bool, optional
      縦軸を対数座標にするかどうか。デフォルトはTrue。
     **zinvert** : bool, optional
      縦軸の向きを逆(減少方向)にするかどうか。デフォルトはTrue。
     **xlabel** : {'lon', 'lat'}, optional
      横軸を経度表示にする場合はlon,緯度表示にする場合はlatを指定する

    :Returns:


    **Examples**

    
    """
    ax = kwargs.get('ax',plt.gca())
    if zinvert:
        ax.set_ylim(z.max(),z.min())
    if zlog: 
        ax.set_yscale('log')
        z = np.asarray(z)
        subs = [1,2,3,4,5,6,7,8,9]
        loc = matplotlib.ticker.LogLocator(base=10.,subs=subs)
        fmt = matplotlib.ticker.FormatStrFormatter("%g")
        ax.yaxis.set_major_locator(loc)
        ax.yaxis.set_major_formatter(fmt)
    if xlabel=='lon':
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticker.lon2txt))
    elif xlabel=='lat':
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticker.lat2txt))

    return ax.plot(xy, z, **kwargs)

