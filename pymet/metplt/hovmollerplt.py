# coding:utf-8
#------------------------------------------------------------
#--- Hovmollerプロット
#------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib.dates
import matplotlib.ticker
from datetime import datetime
import numpy as np
import locale

import ticker

__all__ = ['hovplot','hovcontour','hovcontourf']

def hovplot(xy,time,**kwargs):
    u"""
    ホフメラー図のプロット。

    :Arguments:
     **xy** : array_like
      x or y or 緯度経度座標。
     **time** : array_like of datetime objects
      時間軸。datetimeオブジェクトのarray。
     **xlabel** : {'lon', 'lat'}, optional
      横軸を経度表示にする場合はlon,緯度表示にする場合はlatを指定する
     **fmt** : str, optional
      時間軸ラベルの表示形式。デフォルトは%Hz%d%b%Y。

    :Returns:
     **CR**

    **Examples**      
    """
    locale.setlocale(locale.LC_ALL,'en_US')
    ax = kwargs.get('ax',plt.gca())
    xlabel = kwargs.pop('xlabel',None)
    ax.set_ylim(time.max(),time.min())
    fmt = kwargs.pop('fmt','%Hz%d%b%Y')
    if isinstance(time[0],datetime) :
        if fmt==None:
            ax.yaxis.set_major_formatter(matplotlib.ticker.NullFormatter())
        else:
            ax.yaxis.set_major_formatter(matplotlib.dates.DateFormatter(fmt))
    if xlabel=='lon':
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticker.lon2txt))
    elif xlabel=='lat':
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticker.lat2txt))
    return ax.plot(xy,time,*args,**kwargs)


def hovcontour(xy,time,data,*args,**kwargs):
    u"""
    ホフメラー図のコンタープロット。

    :Arguments:
     **xy** : array_like
      x or y or 緯度経度座標。
     **time** : array_like of datetime objects
      時間軸。datetimeオブジェクトのarray。
     **data** : 2darray
      プロットするデータ。
     **xlabel** : {'lon', 'lat'}, optional
      横軸を経度表示にする場合はlon,緯度表示にする場合はlatを指定する
     **fmt** : str, optional
      時間軸ラベルの表示形式。デフォルトは%Hz%d%b%Y。
     **cint** : float, optional
      コンター間隔。

    :Returns:
     **CR**

    **Examples**
      
    """
    locale.setlocale(locale.LC_ALL,'en_US')
    ax = kwargs.get('ax',plt.gca())
    xlabel = kwargs.pop('xlabel',None)
    cint =  kwargs.pop('cint',None)
    kwargs.setdefault('colors','k')
    ax.set_ylim(time.max(),time.min())
    fmt = kwargs.pop('fmt','%Hz%d%b%Y')
    if isinstance(time[0],datetime):
        ax.yaxis.set_major_formatter(matplotlib.dates.DateFormatter(fmt))
    if xlabel=='lon':
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticker.lon2txt))
    elif xlabel=='lat':
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticker.lat2txt))
    if cint != None:
        kwargs.setdefault('locator',matplotlib.ticker.MultipleLocator(cint))
        ax.text(1,-0.15,'contour interval = %g' % cint, 
                transform=ax.transAxes,ha='right')#,fontsize=10)
    return ax.contour(xy,time,data,*args,**kwargs)

def hovcontourf(xy,time,data,*args,**kwargs):
    u"""
    ホフメラー図の塗りつぶしコンタープロット。

    :Arguments:
     **xy** : array_like
      x or y or 緯度経度座標。
     **time** : array_like of datetime objects
      時間軸。datetimeオブジェクトのarray。
     **data** : 2darray
      プロットするデータ。
     **xlabel** : {'lon', 'lat'}, optional
      横軸を経度表示にする場合はlon,緯度表示にする場合はlatを指定する
     **fmt** : str, optional
      時間軸ラベルの表示形式。デフォルトは%Hz%d%b%Y。
     **cint** : float, optional
      コンター間隔。

    :Returns:
     **CF**

    **Examples**
      
    """
    ax = kwargs.get('ax',plt.gca())
    xlabel = kwargs.pop('xlabel',None)
    ax.set_ylim(time.max(),time.min())
    fmt = kwargs.pop('fmt','%Hz%d%b%Y')
    xlint = kwargs.pop('xlint',None)
    if isinstance(time[0],datetime):
        ax.yaxis.set_major_formatter(matplotlib.dates.DateFormatter(fmt))
    if xlabel=='lon':
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticker.lon2txt))
    elif xlabel=='lat':
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticker.lat2txt))
    if xlint!=None:
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(xlint))
    return ax.contourf(xy,time,data,*args,**kwargs)

