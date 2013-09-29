# coding:utf-8
u"""
"""
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pymet.tools as tools
import ticker
from matplotlib.font_manager import FontProperties
import matplotlib
import scipy.ndimage
rcParams = matplotlib.rcParams

__all__ = ['MyBasemap']

class MyBasemap(Basemap):
    u"""
    Basemapを使いやすくするためのラッパー。

    **Keeyward**
     **projection**
      地図の種類。
      
      ======== ============================================================
      Value    Description
      ======== ============================================================
      cyl      正距円筒図法(Cylindrical Equidistant)。同緯度経度で表示。
      merc     メルカトル図法(Mercator)。角度が正しい図法。
      lcc      ランベルト正角円錐図法(Lambert Conformal Conic)。
      npstere  北極ステレオ投影図法(Polar Stereographic Projection)。
      spstere  南極ステレオ投影図法(South Stereographic Projection)。
      ======== ============================================================

     **地図の範囲**
      ============= =================================================================================
      Keyword       Description
      ============= =================================================================================
      lon           地図の経度の範囲。デフォルトは(-180, 180)。
      lat           地図の緯度の範囲。デフォルトは(-90, 90)。npstere,spstereでは(lat[0],90),
                    (-90,lat[-1])になる。
      lon_0         地図の中心経度。cyl,mercでは無視される。デフォルトは(lon[0]+lon[1])/2。
      lat_0         地図の中心緯度。cyl,mercでは無視される。デフォルトは(lat[0]+lat[1])/2。
      boundinglat   npstere,spstereのときの外周緯度。デフォルトは赤道。
      ============= =================================================================================
      
     **軸、枠**  
      ======== =================================================================================
      Keyword  Description
      ======== =================================================================================
      xlint    経度メモリの間隔(degrees)。
      ylint    緯度メモリの間隔(degrees)。
      grid     緯線,経線を引くかどうか。デフォルトはFalse。
      label    緯度経度のラベルをふるかどうか。デフォルトはTrue
      ======== =================================================================================
      
     **その他**
      ======== =========================================================================
      Keyword  Description
      ======== =========================================================================
      ax       地図を作成するaxesオブジェクトを指定する場合。指定しない場合はplt.gca()。
      ======== ========================================================================= 
     
    **Methods**
     .. currentmodule:: pymet.plot.mybasemap.MyBasemap

     .. autosummary::

        __init__
        fillcontinents
        set_xlint
        set_ylint
        plot
        contour
        contourf
        colorbar
        quiver
        quiverkey


    **Examples**
     >>> m = MyBasemap(lon=(0,180), lat=(0,90), projection='cyl', xlint=30, ylin=30)
     >>>

    .. seealso::
       see more deatail <http://matplotlib.org/basemap/>


        
    """
    def __init__(self,lon=None,lat=None,xlint=None,ylint=None,nocoast=False,
                 grid=False,**keys):
        if lon != None:
            slon = keys.setdefault('llcrnrlon', lon[0])
            elon = keys.setdefault('urcrnrlon', lon[1])
        else:
            slon = keys.setdefault('llcrnrlon', 0.)
            elon = keys.setdefault('urcrnrlon', 360.)
        if lat != None:
            slat = keys.setdefault('llcrnrlat', lat[0])
            elat = keys.setdefault('urcrnrlat', lat[1])                        
        else:
            slat = -90.
            elat = 90.            
        projection = keys.setdefault('projection','cyl')
        label = keys.pop('label', True)
        # 地図に合わせたデフォルト値
        if projection == 'npstere':
            boundinglat = keys.setdefault('boundinglat',max(slat, 0.))
            keys.setdefault('lon_0',0.5*(slon+elon))
            keys.setdefault('suppress_ticks',True)
            keys.setdefault('round',True)
        elif projection == 'spstere':
            keys.setdefault('boundinglat',min(elat, 0.))
            keys.setdefault('lon_0',0.5*(slon+elon))
            keys.setdefault('suppress_ticks',True)
            keys.setdefault('round',True)
        elif projection == 'lcc':
            keys.setdefault('lon_0',0.5*(slon+elon))
#            keys.setdefault('lat_0',0.5*(slat+elat))
            keys.setdefault('suppress_ticks',True)
        if not grid: keys.setdefault('suppress_ticks',False)

        # draw map
        Basemap.__init__(self,**keys)
        self.drawmapboundary()
        if not nocoast: self.drawcoastlines(linewidth=0.5)
        ax = self.ax or self._check_ax()
        if projection == 'npstere' or projection == 'lcc' or grid:
            if xlint is None: xlint=60
            if ylint is None: ylint=20
            if label:
                self.drawmeridians(np.arange(0,360,xlint),labels=[1,1,1,1],linewidth=0.5)
                self.drawparallels(np.arange(0,90.1,ylint),labels=[1,0,0,0],linewidth=0.5)
            else:
                self.drawmeridians(np.arange(0,360,xlint),labels=[0,0,0,0],linewidth=0.5)
                self.drawparallels(np.arange(0,90.1,ylint),labels=[0,0,0,0],linewidth=0.5)
        else:
            if label:
                ax.xaxis.set_major_formatter(ticker.BasemapXaxisFormatter(self))   
                ax.yaxis.set_major_formatter(ticker.BasemapYaxisFormatter(self))
            else:
                ax.xaxis.set_major_formatter(matplotlib.ticker.NullFormatter())
                ax.yaxis.set_major_formatter(matplotlib.ticker.NullFormatter())                   
            if xlint != None:
                self.set_xlint(xlint)
            if ylint != None:
                self.set_ylint(ylint)
            ax.tick_params(direction='out', top='off', right='off')
            
            
    def fillcontinents(self,color='lightgray',zorder=0,**keys):
        u"""
        大陸に色をつける。

        :Arguments:
         **color** : str, optional
          塗りつぶす色。デフォルトはlightgray。
         **zorder** : int, optional
          Artistsオブジェクトを描く順番。小さい番号から描画される。

        **Exmaples**
         zorderをcontourやcontourfで指定したものより大きくすると、大陸マスクのように描画できる。
          >>> m.contourf(lon, lat, data, zorder=0)
          >>> m.fillcontinents(zorder=1)
        
         .. plot:: ../examples/fillcontinents.py
        """
        Basemap.fillcontinents(self,color=color,zorder=zorder,**keys)

    def set_xlint(self,xlint,sx=0.):
        u"""
        x軸目盛りの間隔(degrees)を指定する。

        :Arguments:
         **xlint**: float
          目盛り間隔(degrees)。
         **sx** : float, optional
          スタート位置。デフォルトはゼロ。
        """
        ax = self.ax or self._check_ax()
        xlocs = np.arange(sx,sx+360.0001,xlint)
        ax.xaxis.set_major_locator(ticker.BasemapXaxisLocator(self,xlocs))

    def set_ylint(self,ylint,sy=-90.):
        u"""
        y軸目盛りの間隔(degrees)を指定する。

        :Arguments:
         **ylint**: float
          目盛り間隔(degrees)。
         **sy** : float, optional
          目盛りのスタート位置。デフォルトは-90。          
        """
        ax = self.ax or self._check_ax()
        ylocs = np.arange(sy, 90.001, ylint)
        ax.yaxis.set_major_locator(ticker.BasemapYaxisLocator(self,ylocs))

    def plot(self,lon,lat,*args,**kwargs):
        u"""
        通常のプロット。

        :Arguments:
         **lon** : array_like
          経度
         **lat** : array_like
          緯度
        """
        x, y = self(lon,lat)

        return Basemap.plot(self, x, y, *args, **kwargs)

    def contour(self, lon, lat, data, *args, **kwargs):
        u"""
        コンターをプロットする。

        :Arguments:
         **lon** : array_like
          緯度。
         **lat** : array_like
          経度。
         **data** : 2darray
          プロットするデータ。
         :Returns:
          **CR** :

        **Keyword**
         指定可能なキーワード(一部)
         
         =========== ========= ======================================================
         Value       Default   Description
         =========== ========= ======================================================
         colors      'k'       コンターの色。
         skip        1         データの表示間隔。data[::skip,::skip]が使われる。
         cint        None      コンター間隔。指定しない場合は自動で決定される。
         cintlab     True      図の右下にcontour interval= 'cint unit'の形式で
                               コンター間隔のラベルを表示。cintを指定した場合に有効。
         labunit               cinttext=Trueの場合に表示する単位。
         hoffset     0         コンター間隔を表す文字列の位置のオフセット
         voffset     auto      コンター間隔を表す文字列の位置のオフセット
         =========== ========= ======================================================

        **Examples**
         .. plot:: ../examples/contour.py

        """
        ax = self.ax or self._check_ax()
        kwargs.setdefault('colors', 'k')
        kwargs.setdefault('latlon', True)
        skip = kwargs.pop('skip', 1)

        cint = kwargs.pop('cint', None)        
        if cint != None:
            kwargs.setdefault('locator',matplotlib.ticker.MultipleLocator(cint))
            if kwargs.pop('cintlab',True):
                hoffset = kwargs.pop('hoffset', 0)
                voffset = -FontProperties(size=rcParams['font.size']).get_size_in_points() - rcParams['xtick.major.pad'] - rcParams['xtick.major.size']
                voffset = kwargs.pop('voffset', voffset)                
                labunit = kwargs.pop('labunit', '')
                xlabel_fontsize = ax.xaxis.get_label().get_fontsize()                
                ax.annotate('interval = %g %s' % (cint, labunit),
                            xy=(1, 0), xycoords='axes fraction',
                            xytext=(hoffset, voffset), textcoords='offset points',
                            ha='right',va='top', fontsize=xlabel_fontsize*0.9)
                
        if np.ndim(lon)==1 or np.ndim(lat)==1:
            lon, lat = np.meshgrid(lon, lat)
        else:
            if np.ndim(lon)==2: lon = lon[0,:]
            if np.ndim(lat)==2: lat = lat[:,0]                
            lon, lat = np.meshgrid(lon, lat)
            
        # 境界に応じてサイクリックにする
        if ((self.urcrnrlon%360.==self.llcrnrlon%360. or self.projection in ['npstere', 'spstere'])
                and (lon[0,-1]-lon[0,-2])>=(lon[0,0]+360.-lon[0,-1])):
            lon = np.c_[lon, lon[:,0]+360.]
            lat = np.c_[lat, lat[:,0]]
            data = np.c_[data, data[:,0]]

        if kwargs.pop('zoom', False):
            data[data.mask] = data.mean()
            lon  = scipy.ndimage.zoom(lon, 3)
            lat  = scipy.ndimage.zoom(lat, 3)
            data = scipy.ndimage.zoom(data, 3)

        return Basemap.contour(self, lon[::skip,::skip], lat[::skip,::skip],
                                data[::skip,::skip], *args, **kwargs)

    def contourf(self, lon, lat, data, *args, **kwargs):
        u"""
        塗りつぶしコンターをプロットする。

        :Arguments:
         **lon** : array_like
          緯度。
         **lat** : array_like
          経度。
         **data** : 2darray
          プロットするベクトルのx,y成分。

        :Returns:
         **CF** : QuadContourSet object
          
        **Keyword**
         独自キーワード
                 
         ========== ======= ======================================================
         Value      Default Description
         ========== ======= ======================================================
         skip       1       データの表示間隔。data[::skip,::skip]が使われる。
         cint       None    コンター間隔。指定しない場合は自動で決定される。         
         ========== ======= ======================================================

         basemapと共通のキーワード(デフォルトを独自に設定しているもの)
         
         ========== ======= ======================================================
         Value      Default Description
         ========== ======= ======================================================
         exntend    both
         latlon     True    lon,latを経度,緯度で解釈する場合はTrue
         ========== ======= ======================================================
                  
        **Examples**
         .. plot:: ../examples/contourf.py
        
        """
        skip = kwargs.pop('skip',1)
        kwargs.setdefault('latlon',True)
        kwargs.setdefault('extend', 'both')        
        
        cint = kwargs.pop('cint', None)
        if cint != None:
            kwargs.setdefault('locator',matplotlib.ticker.MultipleLocator(cint))

        if np.ndim(lon)==1 or np.ndim(lat)==1:
            lon, lat = np.meshgrid(lon, lat)
        else:
            if np.ndim(lon)==2: lon = lon[0,:]
            if np.ndim(lat)==2: lat = lat[:,0]                
            lon, lat = np.meshgrid(lon, lat)
            
        if kwargs.pop('zoom', False):
            lon  = scipy.ndimage.zoom(lon, 3)
            lat  = scipy.ndimage.zoom(lat, 3)
            data = scipy.ndimage.zoom(data, 3)
            
        # 境界に応じてサイクリックにする
        if ((self.urcrnrlon%360.==self.llcrnrlon%360. or self.projection in ['npstere', 'spstere'])
                and (lon[0,-1]-lon[0,-2])>=(lon[0,0]+360.-lon[0,-1])):
            lon = np.c_[lon, lon[:,0]+360.]
            lat = np.c_[lat, lat[:,0]]
            data = np.c_[data, data[:,0]]        

        return Basemap.contourf(self, lon[::skip,::skip], lat[::skip,::skip],
                                data[::skip,::skip], *args,**kwargs)

    def quiver(self,lon,lat,u,v,*args,**kwargs):
        u"""
        ベクトルをプロットする。

        :Arguments:
         **lon** : array_like
          緯度。
         **lat** : array_like
          経度。
         **u, v** : 2darray
          プロットするベクトルデータ。

        :Returns: 
         **QV** 

        **Keyword**
         独自キーワード
         
         =============== ======= ======================================================
         Value           Default Description
         =============== ======= ======================================================
         skip            1       データの表示間隔。data[::skip,::skip]が使われる。
         zoom            False   smoothingをかけるかどうか
         data_scale              基準長(dot_scale)あたりの実データの長さ。デフォルトは
                                 max(sqrt(u**2+v**2))                                   
         dot_scale       30      基準長さ(ドット,ピクセル)         
         =============== ======= ======================================================

         デフォルトを独自に設定しているキーワード
         
         =============== ======= ======================================================
         Value           Default Description
         =============== ======= ======================================================
         scale_units     width
         scale                   scale_unitsの1単位に対するデータ長
         headwidth
         headaxislength
         headlength
         linewidth
         width
         zorder
         =============== ======= ======================================================
         
        **Examples**
         .. plot:: ../examples/quiver.py
         
        """
        ax = self.ax or self._check_ax()
        skip = kwargs.pop('skip',1)

        if np.ndim(lon)==1 or np.ndim(lat)==1:
            lon, lat = np.meshgrid(lon, lat)
        else:
            if np.ndim(lon)==2: lon = lon[0,:]
            if np.ndim(lat)==2: lat = lat[:,0]                
            lon, lat = np.meshgrid(lon, lat)
            
        # 境界に応じてサイクリックにする
        if ((self.urcrnrlon%360.==self.llcrnrlon%360. or self.projection in ['npstere', 'spstere'])
                and (lon[0,-1]-lon[0,-2])>=(lon[0,0]+360.-lon[0,-1])):
            lon = np.c_[lon, lon[:,0]+360.]
            lat = np.c_[lat, lat[:,0]]
            u = np.ma.c_[u, u[:,0]]
            v = np.ma.c_[v, v[:,0]]

        # 緯度経度座標で与えられたu,vを地図座標での座標方向に合わせて回転
        kwargs.setdefault('latlon', True)
        if self.projection!='cyl':
            u, v = self.rotate_vector(u,v,lon,lat)
            
        # スケールを決める
        ax_width = ax.bbox.width                     # 図の横幅(dot)        
        data_scale = kwargs.pop('data_scale', None)        
        dot_scale  = kwargs.pop('dot_scale', ax_width/10.)
        
        if not kwargs.has_key('scale_units'):
            kwargs.setdefault('scale_units', 'width')            
            # データの最大値を有効数字1桁で丸めてデータの基準長とする
            if data_scale is None:
                mask = (lon <= self.urcrnrlon) & (lon >= self.llcrnrlon) & (lat <= self.urcrnrlat) & (lat >= self.llcrnrlat)
                data_scale = tools.roundoff(np.ma.hypot(u[mask], v[mask]).max(), digit=1)
            scale = data_scale * ax_width / dot_scale           # (test) dot_scale(dot)がdata_scaleになるようにする
            kwargs.setdefault('scale', scale)

        # 矢印の形状の設定
        width = kwargs.setdefault('width',0.5/ax_width)    # shaft width (ax_width*<width> dot) => default 0.5pt        
        kwargs.setdefault('headwidth',10/width/ax_width)   # headwidth   width*(<dot>/width/ax_width)
        kwargs.setdefault('headaxislength',0)
        kwargs.setdefault('headlength',5/width/ax_width)       
        kwargs.setdefault('linewidth',0.5/width/ax_width)

        # 色
        color = kwargs.get('color', None)
        if color and not 'facecolor' in kwargs and not 'edgecolor' in kwargs:
            kwargs['facecolor'] = color
            kwargs['edgecolor'] = color
                    
        QV = Basemap.quiver(self,lon[::skip,::skip],lat[::skip,::skip],
                            u[::skip,::skip], v[::skip,::skip],*args, **kwargs)

        # quiverkeyのために基準長さを残しておく
        if data_scale!=None:
            QV.data_scale = data_scale

        return QV
            
    def colorbar(self,mappable=None,location='right',size="5%",pad='2%',fig=None,ax=None,**kwargs):
        u"""
        カラーバーを描画する。        
        """
        clabel = kwargs.pop('clabel',None)
        clabelsize = kwargs.pop('clabelsize',None)
        fontsize = kwargs.pop('fontsize',None)
        CB = Basemap.colorbar(self,mappable,location,size,pad,fig,ax,**kwargs)
        if clabel!=None:
            CB.set_label(clabel,fontsize=clabelsize)
        if fontsize!=None:
            if CB.orientation=='horizontal':
                for t in CB.ax.get_xticklabels():
                    t.set_fontsize(fontsize)
            else:
                for t in CB.ax.get_yticklabels():
                    t.set_fontsize(fontsize)
        return CB                        
         
    def quiverkey(self, QV, **kwargs):
        u"""
        ベクトルプロットのreferrence arrowを描画する。

        :Arguments:
         **QV** :

        **Keyword**
         独自キーワード
         
         =============== =========== ======================================================
         Value           Default     Description
         =============== =========== ======================================================
         data_scale                  referrence arrowのデータ座標での長さ
         unit                        referrence arrowにつける単位
         =============== =========== ======================================================

         デフォルトを独自に設定しているキーワード
         
         =============== ======= ======================================================
         Value           Default Description
         =============== ======= ======================================================
         labelpos        'N'
         =============== ======= ======================================================
        """
        ax = self.ax or self._check_ax()

        # referrence arrow の実データでの長さの優先順位
        # 1. data_scaleで与える
        # 2. QVのdata_scale
        # 3. 最大値
        data_scale = kwargs.pop('data_scale', None)
        xlabel_fontsize = ax.xaxis.get_label().get_fontsize()                        
        if data_scale==None:
            data_scale = getattr(QV, 'data_scale', tools.roundoff(np.sqrt(QV.U**2 + QV.V**2).max(), digit=1))

        # referrence arrow の中心位置の決定        
        labelpos = kwargs.get('labelpos', 'N')
        coordinates = kwargs.setdefault('coordinates', 'axes')
        if coordinates == 'axes':
            width_dot  = ax.bbox.width  # (dot)
            height_dot = ax.bbox.height # (dot)
            data_scale_dot = data_scale / QV.scale * width_dot # (dot)
            if labelpos == 'N':
                loc = kwargs.pop('loc', (1-0.5*data_scale_dot/width_dot, 1+10/height_dot))
                kwargs.setdefault('labelsep', 10/height_dot)
            elif labelpos == 'W':
                loc = kwargs.pop('loc', (1-0.5*data_scale_dot/width_dot, 1+10/height_dot))
                kwargs.setdefault('labelsep', 10/width_dot)

        # fontpropeties
        fontproperties = kwargs.pop('fontproperties', {})
        for key in kwargs.keys():
            if key.find('font')>0:
                fontproperties[key[4:]] = kwargs.pop(key)
        xlabel_fontsize = ax.xaxis.get_label().get_fontsize()
        fontproperties.setdefault('size', xlabel_fontsize*0.9)
        kwargs['fontproperties'] = fontproperties

        unit = kwargs.pop('unit', '')
        ax.quiverkey(QV,loc[0],loc[1],data_scale,'%g %s'%(data_scale,unit),**kwargs)

    def drawbox(self,lon1,lat1,lon2,lat2,**kwargs):
        u"""
        地図上にボックスを描く。

        :Arguments:
         **lon1, lat1** : float
          ボックスの左下の座標(degree)。
         **lon2, lat2** : float
          ボックスの右上の座標(degree)。
        :Returns:

        **Examples**
         >>>
        """
        kwargs.setdefault('color','k')
        kwargs.setdefault('lw', 2)
        self.plot([lon1,lon2,lon2,lon1,lon1],[lat1,lat1,lat2,lat2,lat1],**kwargs)
        
def mapinterp(datain, lonin, latin, lonout, latout):
    if lonin.ndim == 2:
        lonin = lonin[0,:]
    if latin.ndim == 2:
        latin = latin[:,0]
    return mpl_toolkits.basemap.interp(datain, lonin, latin, lonout, latout)    
    
