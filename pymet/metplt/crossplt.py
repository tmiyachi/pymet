# coding:utf-8
#------------------------------------------------------------------------------
# 鉛直断面図のための関数
#------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

import ticker

__all__ = ['crossplot','crosscontour','crosscontourf']


def crossplot(xy, z, **kwargs):
    u"""
    鉛直断面プロット。

    :Arguments:
     **xy** : array_like
      x or y座標。
     **z** : array_like
      鉛直座標。
    :Returns:
     
    **Keyword**
     =========== ========= ======================================================
     Value       Default   Description
     =========== ========= ======================================================
     zinvert     True      縦軸(鉛直レベル)を降順にするかどうか
     zlog        True      縦軸(鉛直レベル)を対数座標にするかどうか
     xylab       None      横軸を経度表示にする場合は'lon',
                           緯度表示にする場合は'lat'を指定する
     =========== ========= ======================================================

    **Examples**

    
    """
    ax = kwargs.get('ax',plt.gca())
    
    if kwargs.pop('zinvert', True):
        ax.set_ylim(z.max(),z.min())
        
    if kwargs.pop('zlog', True): 
        ax.set_yscale('log')
        z = np.asarray(z)
        subs = [1,2,3,4,5,6,7,8,9]
        loc = matplotlib.ticker.LogLocator(base=10.,subs=subs)
        fmt = matplotlib.ticker.FormatStrFormatter("%g")
        ax.yaxis.set_major_locator(loc)
        ax.yaxis.set_major_formatter(fmt)
        
    xylabel = kwargs.pop('xylabel', None)
    if xlabel=='lon':
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticker.lon2txt))
    elif xlabel=='lat':
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticker.lat2txt))

    return ax.plot(xy, z, **kwargs)


def crosscontour(xy, z, data, **kwargs):
    u"""
    東西・南北-鉛直断面のコンタープロット。

    :Arguments:
     **xy** : array_like
      東西、南北座標または緯度・経度。
     **z** : array_like
      鉛直座標。
     **data** : 2darray
      プロットするデータ。
    :Returns:
     **CR**

    **Keyword**
     =========== ========= ======================================================
     Value       Default   Description
     =========== ========= ======================================================
     zinvert     True      縦軸(鉛直レベル)を降順にするかどうか
     zlog        True      縦軸(鉛直レベル)を対数座標にするかどうか
     colors      'k'       コンターの色。
     cint        None      コンター間隔。指定しない場合は自動で決定される。
     cintlab     True      図の右下にcontour interval= 'cint unit'の形式で
                           コンター間隔の説明を表示。cintを指定した場合に有効。
     labunit     ''        cinttext=Trueの場合に表示する単位。
     labxtext     0        コンター間隔を表す文字列の位置
     labytext    -15       コンター間隔を表す文字列の位置
     xylab       None      横軸を経度表示にする場合は'lon',
                           緯度表示にする場合は'lat'を指定する
     =========== ========= ======================================================

    **Examples**

    """
    ax = kwargs.get('ax',plt.gca())
    kwargs.setdefault('colors','k')
    
    if kwargs.pop('zinvert', True):
        ax.set_ylim(z.max(),z.min())
        
    if kwargs.pop('zlog', True): 
        ax.set_yscale('log')
        z = np.asarray(z)
        subs = [1,2,3,4,5,6,7,8,9]
        loc = matplotlib.ticker.LogLocator(base=10.,subs=subs)
        fmt = matplotlib.ticker.FormatStrFormatter("%g")
        ax.yaxis.set_major_locator(loc)
        ax.yaxis.set_major_formatter(fmt)
        
    xylabel = kwargs.pop('xylabel', None)
    if xylabel=='lon':
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticker.lon2txt))
    elif xylabel=='lat':
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticker.lat2txt))

    cint = kwargs.pop('cint', None)        
    if cint != None:
        kwargs.setdefault('locator',matplotlib.ticker.MultipleLocator(cint))
        if kwargs.pop('cintlab',True):
            labxtext = kwargs.pop('labxtext', 0)
            labytext = kwargs.pop('labytext', -15)
            labunit = kwargs.pop('labunit', '')
            ax.annotate('contour interval = %g %s' % (cint, labunit),
                        xy=(1, 0), xycoords='axes fraction',
                        xytext=(labxtext, labytext), textcoords='offset points',
                        ha='right',va='top')
    
    return ax.contour(xy,z,data,**kwargs)

def crosscontourf(xy, z, data, **kwargs):
    u"""
    東西・南北-鉛直断面の塗りつぶしコンタープロット。

    :Arguments:
     **xy** : array_like
      東西、南北座標または緯度・経度。
     **z** : array_like
      鉛直座標。
     **data** : 2darray
      プロットするデータ。
    :Returns:
     **CF**

    **Keyword**
     =========== ========= ======================================================
     Value       Default   Description
     =========== ========= ======================================================
     zinvert     True      縦軸(鉛直レベル)を降順にするかどうか
     zlog        True      縦軸(鉛直レベル)を対数座標にするかどうか
     colors      'k'       コンターの色。
     cint        None      コンター間隔。指定しない場合は自動で決定される。
     xylab       None      横軸を経度表示にする場合は'lon',
                           緯度表示にする場合は'lat'を指定する
     =========== ========= ======================================================

    **Examples**
     .. plot:: ../examples/crosscontour.py              
    """
    ax = kwargs.get('ax',plt.gca())
    if kwargs.pop('zinvert', True):
        ax.set_ylim(z.max(),z.min())
        
    if kwargs.pop('zlog', True): 
        ax.set_yscale('log')
        subs = [1,2,3,4,5,6,7,8,9]
        loc = matplotlib.ticker.LogLocator(base=10.,subs=subs)
        fmt = matplotlib.ticker.FormatStrFormatter("%g")
        ax.yaxis.set_major_locator(loc)
        ax.yaxis.set_major_formatter(fmt)
        
    cint = kwargs.pop('cint', None)
    if cint != None:
        kwargs.setdefault('locator',matplotlib.ticker.MultipleLocator(cint))

    xylabel = kwargs.pop('xylabel', None)
    if xylabel=='lon':
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticker.lon2txt))
    elif xylabel=='lat':
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticker.lat2txt))
    return ax.contourf(xy,z,data,**kwargs)

