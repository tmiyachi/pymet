# coding:utf-8
u"""
"""
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

import ticker

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
      ======== =================================================================================
      Keyword  Description
      ======== =================================================================================
      lon      地図の経度の範囲。デフォルトは(-180, 180)。
      lat      地図の緯度の範囲。デフォルトは(-90, 90)。npstere,spstereでは(lat[0],90),
               (-90,lat[-1])になる。
      lon_0    地図の中心経度。cyl,mercでは無視される。デフォルトは(lon[0]+lon[1])/2。
      lat_0    地図の中心緯度。cyl,mercでは無視される。デフォルトは(lat[0]+lat[1])/2。
      ======== =================================================================================
      
     **軸、枠**  
      ======== =================================================================================
      Keyword  Description
      ======== =================================================================================
      xlint    経度メモリの間隔(degrees)。
      ylint    緯度メモリの間隔(degrees)。
      grid     
      ======== =================================================================================

     **その他**
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
     >>>
     >>>

    .. seealso::
       see more deatail <http://matplotlib.org/basemap/>


        
    """
    def __init__(self,lon=None,lat=None,xlint=None,ylint=None,nocoast=False,
                 grid=False,**keys):
        if lon != None:
            keys['llcrnrlon']=lon[0]
            keys['urcrnrlon']=lon[1]
        if lat != None:
            keys['llcrnrlat']=lat[0]
            keys['urcrnrlat']=lat[1]

        projection = keys.setdefault('projection','cyl')

        # 地図に合わせたデフォルト値
        if projection == 'merc':
            keys['llcrnrlat']=lat[0]-0.00001
        elif projection == 'npstere':
            keys.setdefault('boundinglat',max(lat[0], 0.))
            keys.setdefault('lon_0',0.5*(lon[0]+lon[1]))
            keys.setdefault('suppress_ticks',True)
            keys.setdefault('round',True)
        elif projection == 'spstere':
            keys.setdefault('boundinglat',min(lat[1], 0.))
            keys.setdefault('lon_0',0.5*(lon[0]+lon[1]))
            keys.setdefault('suppress_ticks',True)
            keys.setdefault('round',True)
        elif projection == 'lcc':
            keys.setdefault('lon_0',0.5*(lon[0]+lon[1]))
            keys.setdefault('lat_0',0.5*(lat[0]+lat[1]))
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
            self.drawmeridians(np.arange(0,360,xlint),labels=[0,0,0,1],linewidth=0.5)
            self.drawparallels(np.arange(0,90.1,ylint),labels=[1,0,0,0],linewidth=0.5)
        else:
            ax.xaxis.set_major_formatter(ticker.BasemapXaxisFormatter(self))   
            ax.yaxis.set_major_formatter(ticker.BasemapYaxisFormatter(self))
            if xlint != None:
                self.set_xlint(xlint)
            if ylint != None:
                self.set_ylint(ylint)
       
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

    def set_xlint(self,xlint):
        u"""
        x軸目盛りの間隔(degrees)を指定する。

        :Arguments:
         **xlint**: float
          目盛り間隔(degrees)。
        """
        ax = self.ax or self._check_ax()
        xlocs = np.arange(0,360.0001,xlint)
        ax.xaxis.set_major_locator(ticker.BasemapXaxisLocator(self,xlocs))

    def set_ylint(self,ylint):
        u"""
        y軸目盛りの間隔(degrees)を指定する。

        :Arguments:
         **ylint**: float
          目盛り間隔(degrees)。
        """
        ax = self.ax or self._check_ax()
        ylocs = np.arange(-90, 90.001, ylint)
        ax.yaxis.set_major_locator(ticker.BasemapYaxisLocator(self,ylocs))

    def plot(self,lon,lat,*args,**kwargs):
        u"""
        プロット。
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


        **Keyword**
         指定可能なキーワード(一部)
         
         ========== ======= ======================================================
         Value      Default Description
         ========== ======= ======================================================
         colors     'k'     コンターの色。
         skip       1       データの表示間隔。data[::skip,::skip]が使われる。
         cint       None    コンター間隔。指定しない場合は自動で決定される。
         cinttext   True    図の右下にcontour interval= 'cint unit'の形式で
                            コンター間隔の説明を表示。cintを指定した場合に有効。
         unit       ''      cinttext=Trueの場合に表示する単位。
         textpos
         zorder
         ========== ======= ======================================================

        **Examples**
         .. plot:: ../examples/contour.py
        """
        ax = self.ax or self._check_ax()
        kwargs.setdefault('colors','k')
        kwargs.setdefault('latlon',True)
        skip = kwargs.pop('skip',1)
        cint = kwargs.pop('cint',None)
        cinttext = kwargs.pop('cinttext',True)
        textpos = kwargs.pop('textpos',(1,-0.1))
        unit = kwargs.pop('unit','')
        if cint != None:
            kwargs.setdefault('locator',matplotlib.ticker.MultipleLocator(cint))
            if cinttext:
                ax.text(textpos[0],textpos[1],'contour interval = %g %s' % (cint, unit), 
                        transform=ax.transAxes,ha='right')

        if (lon.ndim<2 and lat.ndim<2): 
            lon, lat = np.meshgrid(lon, lat)
            
        # 境界に応じてサイクリックにする
        if self.urcrnrlon%360. == self.llcrnrlon%360.:
            lon = np.c_[lon, lon[:,0]+360.]
            lat = np.c_[lat, lat[:,0]]
            data = np.c_[data, data[:,0]]        

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
         exntend    'both'
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

        if (lon.ndim<2 and lat.ndim<2): 
            lon, lat = np.meshgrid(lon, lat)

        # 境界に応じてサイクリックにする
        if self.urcrnrlon%360. == self.llcrnrlon%360.:
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
          プロットするデータ。

        :Returns: 
         **QV** 

        **Keyword**
         独自キーワード
         
         =============== ======= ======================================================
         Value           Default Description
         =============== ======= ======================================================
         skip            1       データの表示間隔。data[::skip,::skip]が使われる。
         =============== ======= ======================================================

         デフォルトを独自に設定しているキーワード
         
         =============== ======= ======================================================
         Value           Default Description
         =============== ======= ======================================================
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
        if (lon.ndim<2 and lat.ndim<2): 
            lon, lat = np.meshgrid(lon, lat)

        # 境界に応じてサイクリックにする
        if self.urcrnrlon%360. == self.llcrnrlon%360.:
            lon = np.c_[lon, lon[:,0]+360.]
            lat = np.c_[lat, lat[:,0]]
            u = np.c_[u, u[:,0]]
            v = np.c_[v, v[:,0]]

        ax = self.ax or self._check_ax()
        skip = kwargs.pop('skip',1)
        #kwargs.setdefault('latlon',True)
        kwargs.setdefault('headwidth',15)
        kwargs.setdefault('headaxislength',0)
        kwargs.setdefault('headlength',10)
        kwargs.setdefault('linewidth',0.5)
        kwargs.setdefault('width',0.001)
        return Basemap.quiver(self,lon[::skip,::skip],lat[::skip,::skip],
                              u[::skip,::skip],v[::skip,::skip],*args,**kwargs)
    
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

         **unit** : str, optional
          referrence arrowの値につける単位。
         **loc** : tuple, optional
          referrence arrowの位置。
        """
        ax = self.ax or self._check_ax()
        unit = kwargs.pop('unit', '')
        loc = kwargs.pop('loc', (0.95, 1.05))
        reflen = np.abs(QV.U+1j*QV.V).max()
        reflen = np.around(reflen,decimals=-int(np.log10(reflen)))
        ax.quiverkey(QV,loc[0],loc[1],reflen,'%g %s'%(reflen,unit),labelpos='W')

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
    
