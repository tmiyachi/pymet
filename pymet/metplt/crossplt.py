# coding:utf-8
#------------------------------------------------------------------------------
# 鉛直断面図のための関数
#------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pymet.tools as tools
import ticker
from matplotlib.font_manager import FontProperties
import matplotlib
rcParams = matplotlib.rcParams

__all__ = ['crossplot','crosscontour','crosscontourf',
           'crossquiver', 'crossquiverkey']

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
        
    xylab = kwargs.pop('xylab', None)
    if xylab=='lon':
        ax.xaxis.set_major_formatter(ticker.XaxisFormatter())
    elif xylab=='lat':
        ax.xaxis.set_major_formatter(ticker.YaxisFormatter())

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
     xylab       None      横軸を経度表示にする場合は'lon',
                           緯度表示にする場合は'lat'を指定する
     hoffset     0         コンター間隔を表す文字列の位置のオフセット
     voffset     auto      コンター間隔を表す文字列の位置のオフセット
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
        
    xylab = kwargs.pop('xylab', None)
    if xylab=='lon':
        ax.xaxis.set_major_formatter(ticker.XaxisFormatter())
    elif xylab=='lat':
        ax.xaxis.set_major_formatter(ticker.YaxisFormatter())

    cint = kwargs.pop('cint', None)        
    if cint != None:
        kwargs.setdefault('locator',matplotlib.ticker.MultipleLocator(cint))
        if kwargs.pop('cintlab',True): 
            hoffset = kwargs.pop('hoffset', 0)
            voffset = -FontProperties(size=rcParams['font.size']).get_size_in_points() - rcParams['xtick.major.pad'] - rcParams['xtick.major.size']
            voffset = kwargs.pop('voffset', voffset)                
            labunit = kwargs.pop('labunit', '')
            ax.annotate('contour interval = %g %s' % (cint, labunit),
                        xy=(1, 0), xycoords='axes fraction',
                        xytext=(hoffset, voffset), textcoords='offset points',
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

     basemapと共通のキーワード(デフォルトを独自に設定しているもの)
    
     ========== ======== ======================================================
     Value      Default  Description
     ========== ======== ======================================================
     exntend    'both'
     ========== ======== ======================================================
     
    **Examples**
     .. plot:: ../examples/crosscontourf.py
     
    """
    ax = kwargs.get('ax',plt.gca())
    kwargs.setdefault('extend', 'both')    
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

    xylab = kwargs.pop('xylab', None)
    if xylab=='lon':
        ax.xaxis.set_major_formatter(ticker.XaxisFormatter())
    elif xylab=='lat':
        ax.xaxis.set_major_formatter(ticker.YaxisFormatter())

    return ax.contourf(xy,z,data,**kwargs)

def crossquiver(xy, z, uv, w, **kwargs):
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
     zscale      最大値比  
     xylab       None      横軸を経度表示にする場合は'lon',
                           緯度表示にする場合は'lat'を指定する
     xyskip      1         データの間引き間隔。data[::xyskip,::zskip]
     zskip       1         データの間引き間隔。data[::xyskip,::zskip]     
     =========== ========= ======================================================

    **Examples**
     .. plot:: ../examples/crossquiver.py              
    """
    ax = kwargs.get('ax',plt.gca())
    # 縦軸を降順にする
    if kwargs.pop('zinvert', True):
        ax.set_ylim(z.max(),z.min()*0.91)
        w = -w
    # grid arrayをつくる
    if xy.ndim<2 and z.ndim<2:
        xy, z = np.meshgrid(xy, z)
        
    # 鉛直成分の水平成分に対するスケールの決定
    zscale = kwargs.pop('zscale', None)
    if zscale==None:
        if np.all(w == 0) or np.all(uv ==0):
            zscale = 1.
        else:
            # 絶対値最大値の比を小数点0桁で丸め込み
            zscale = tools.roundoff(np.abs(uv).max()/np.abs(w).max(), digit=1)
    w = w*zscale
    
    if kwargs.pop('zlog', True): 
        ax.set_yscale('log')
        subs = [1,2,3,4,5,6,7,8,9]
        loc = matplotlib.ticker.LogLocator(base=10.,subs=subs)
        fmt = matplotlib.ticker.FormatStrFormatter("%g")
        ax.yaxis.set_major_locator(loc)
        ax.yaxis.set_major_formatter(fmt)
        
    xylab = kwargs.pop('xylab', None)
    if xylab=='lon':
        ax.xaxis.set_major_formatter(ticker.XaxisFormatter())
    elif xylab=='lat':
        ax.xaxis.set_major_formatter(ticker.YaxisFormatter())
    ax.set_xlim(xy[0,0], xy[0,-1])

    xyskip = kwargs.pop('xyskip',1)
    zskip = kwargs.pop('zskip',1)

    ## scaleは指定しない場合は、quiverによって決まる
    # スケールを決める
    data_scale = kwargs.pop('data_scale', None)
    dot_scale  = kwargs.pop('dot_scale', 30.)   
    if not kwargs.has_key('scale_units'):
        kwargs.setdefault('scale_units', 'width')
        ax_width = ax.bbox.width                     # 図の横幅(dot)
        # データの最大値を有効数字1桁で丸めてデータの基準長とする
        if data_scale == None:
            data_scale = tools.roundoff(np.hypot(uv[::zskip,::xyskip], w[::zskip,::xyskip]).max(), digit=1)
            scale = data_scale * ax_width / dot_scale           # (test) dot_scale(dot)がdata_scaleになるようにする
            kwargs.setdefault('scale', scale)
    
    kwargs.setdefault('headwidth',15)
    kwargs.setdefault('headaxislength',0)
    kwargs.setdefault('headlength',10)
    kwargs.setdefault('linewidth',0.5)
    kwargs.setdefault('width',0.001)

    QV =  ax.quiver(xy[::zskip,::xyskip], z[::zskip,::xyskip],
                    uv[::zskip,::xyskip], w[::zskip,::xyskip], **kwargs)

    # quiverkeyのために基準長さを残しておく
    if data_scale!=None:
        QV.data_scale = data_scale
    QV.zscale = zscale

    return QV

def crossquiverkey(QV, **kwargs):
    u"""
    鉛直ベクトルプロットのrefference arrowを描く

    :Arguments:
     **QV**
     
    **Keyword**
     独自キーワード
         
     =============== =========== ======================================================
     Value           Default     Description
     =============== =========== ======================================================
     data_scale
     zscale          QV.zscale   水平成分に対する鉛直成分スケール比
     dot_scale       30
     uvunit                      水平成分のreferrence arrowにつける単位
     wunit                       鉛直成分のreferrence arrowにつける単位
     =============== =========== ======================================================

     デフォルトを独自に設定しているキーワード
         
     =============== ======= ======================================================
     Value           Default Description
     =============== ======= ======================================================
     
     =============== ======= ======================================================
    """
    from matplotlib.text import OffsetFrom
    ax = kwargs.pop('ax', plt.gca())
    uvunit = kwargs.pop('uvunit', '')
    wunit = kwargs.pop('wunit', '')
    dot_scale = kwargs.pop('dot_scale', 30)        

    ## Refference Allowの長さ (データ単位)    
    data_scale = kwargs.pop('data_scale', None)
    zscale = kwargs.pop('zscale', None)
    if zscale==None: zscale = getattr(QV, 'zscale', 1)
    if data_scale==None:
        data_scale = getattr(QV, 'data_scale', tools.roundoff(np.hypot(QV.U, QV.V).max(), digit=1))

    offset_from = OffsetFrom(ax, (0.95,1))
    uvlabel = '%.1f %s' % (data_scale, uvunit)
    wlabel = '%.1f %s' % (data_scale/zscale, wunit)
    ax.annotate('', xy=(0,5), xycoords=offset_from,
                xytext=(-dot_scale,dot_scale+5), textcoords=offset_from,
                ha='right', va='bottom',
                arrowprops=dict(arrowstyle='<->',
                            connectionstyle='angle,rad=0'))
    ax.annotate(uvlabel, xy=(1,1), xycoords='axes fraction',
                xytext=(3,5), textcoords='offset points',
                ha='right', va='bottom')
    ax.annotate(wlabel, xy=(0,5), xycoords=offset_from,
                xytext=(-dot_scale*0.8,dot_scale*0.8), textcoords=offset_from,
                ha='left', va='bottom')
    

