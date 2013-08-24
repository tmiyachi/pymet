# coding:utf-8
u"""
=======================================================
格子点データ内挿モジュール (:mod:`pymet.interpolate`)
=======================================================


========

.. autosummary::
   interp_center

"""
import numpy as np
from scipy.interpolate import RectBivariateSpline, griddata

def interp_center(x, y, data, cx, cy, rmax_domain, rint=None, mode='linear'):
    u"""
    (cx,cy)を中心とした相対座標に内挿する

    :Arguments:
     **x, y** : array_like
      x,y座標
     **data** : ndarray
      データ配列
     **cx,cy** : float
      内挿の中心座標
     **:rmax_domain**
      出力する座標の範囲. x=[cx-rmax_domain,cx+rmax_domain],
      y=[cy-rmax_domain,cy+rmax_domain]が出力範囲になる．
     **rlint** : float, optional
      出力する座標の間隔．指定しない場合は入力座標に応じて決定．
     **mode** : {'linear','cubic'}, optional
      内挿方法

    :Returns:
     **xs_out, ys_out** : 2darray
      出力座標．(cx,cy)からの相対座標．
     **out** : ndarray
      内挿したデータ配列
    """
    if np.ndim(x) == 2: x = x[0,:]
    if np.ndim(y) == 2: y = y[:,0]

    ndim = np.ndim(data)
    if ndim > 2:
        oldshape = data.shape        
        data     = data.reshape(-1, oldshape[-2], oldshape[-1])
    elif ndim < 2:
        raise ValueError, "input data dimension size must be larger than 2. input ndim is {}".format(ndim)
        
    # make output grid
    if rint is None:
        rint = max(np.diff(x).max(), np.diff(y).max())
    rmax  = rmax_domain - rmax_domain%rint
    x_out = np.arange(-rmax, rmax+rint/2., rint)
    y_out = np.arange(-rmax, rmax+rint/2., rint)
    x_in  = x - cx
    y_in  = y - cy
                      
    # interpolate
    if ndim == 2:
        if mode == 'linear':
            out = RectBivariateSpline(y_in, x_in, data, kx=1, ky=1)(y_out, x_out)
        elif mode == 'cubic':
            out = RectBivariateSpline(y_in, x_in, data, kx=3, ky=3)(y_out, x_out)
    else:
        zn, _, _ = data.shape
        xn, yn   = len(x_out), len(y_out)
        out = np.ma.empty((zn,yn,xn), dtype=np.float32)
        for z in range(zn):
            if mode == 'linear':
                out[z,:,:] = RectBivariateSpline(y_in, x_in, data[z,:,:], kx=1, ky=1)(y_out, x_out)
            elif mode == 'cubic':
                out[z,:,:] = RectBivariateSpline(y_in, x_in, data[z,:,:], kx=3, ky=3)(y_out, x_out)
        out = out.reshape(oldshape[:-2] + (yn,xn))
        
    xs_out, ys_out = np.meshgrid(x_out, y_out)
    
    return xs_out, ys_out, out

def interp_cross(x_in, y_in, data, x0, x1, y0, y1, rint=None, mode='linear'):
    u"""
    (x0,y0)-(x1,y1)に沿った直線上の点に内挿する
    """
    if np.ndim(x_in) == 2: x_in = x_in[0,:]
    if np.ndim(y_in) == 2: y_in = y_in[:,0]

    ndim = np.ndim(data)
    if ndim > 2:
        oldshape = data.shape        
        data     = data.reshape(-1, oldshape[-2], oldshape[-1])
    elif ndim < 2:
        raise ValueError, "input data dimension size must be larger than 2. input ndim is {}".format(ndim)
        
    # make output grid
    if rint is None:
        rint = max(np.diff(x_in).max(), np.diff(y_in).max())
    rmax  = np.hypot(x1-x0, y1-y0)
    nr    = rmax//rint + 1
    theta = np.arctan2(y1-y0, x1-x0)
    r_out = np.linspace(0, rmax, nr)
    x_out = r_out*np.cos(theta) + x0
    y_out = r_out*np.sin(theta) + y0

    xrev = yrev = 1
    if x0 > x1:
        x_out = x_out[::-1]
        xrev  = -1
    if y0 > y1:
        y_out = y_out[::-1]
        yrev  = -1
        
    # interpolate    
    if ndim == 2:
        if mode == 'linear':
            out = RectBivariateSpline(y_in, x_in, data, kx=1, ky=1)(y_out, x_out)[::xrev,::yrev].diagonal()
        elif mode == 'cubic':
            out = RectBivariateSpline(y_in, x_in, data, kx=3, ky=3)(y_out, x_out)[::xrev,::yrev].diagonal()
    else:
        zn, yn, xn = data.shape
        out = np.ma.empty((zn,nr), dtype=np.float32)
        for z in range(zn):
            if mode == 'linear':
                out[z,:] = RectBivariateSpline(y_in, x_in, data[z,:,:], kx=1, ky=1)(y_out, x_out)[::xrev,::yrev].diagonal()
            elif mode == 'cubic':
                out[z,:] = RectBivariateSpline(y_in, x_in, data[z,:,:], kx=3, ky=3)(y_out, x_out)[::xrev,::yrev].diagonal()
        out = out.reshape(oldshape[:-2] + (-1,))
                
    return x_out, y_out, r_out, out
    
def grid2cyl(x, y, value, r, theta, method='linear'):
    u"""
    格子座標 (x,y) から(x,y) を中心とした
    円筒座標系 (r,theta) へのスカラーの内挿を行う
    """
    x = np.asarray(x)
    y = np.asarray(y)
    value = np.asarray(value)
    
    # Check input array size
    if np.ndim(x) > 2 or np.ndim(y) > 2:
        raise ValueError, "x, y is must be 2d grid array or 1d array"
    if np.ndim(x) == 1 and np.ndim(y) ==  2:        
        x, y = np.meshgrid(x, y[:,0])
    elif np.ndim(x) == 2 and np.ndim(y) ==  1:
        x, y = np.meshgrid(x[0,:], y)
    elif np.ndim(x) == 1 and np.ndim(y) == 1:
        x, y = np.meshgrid(x, y)

    rs, thetas = np.meshgrid(r, theta)

    # 円筒座標の局所直交座標表現
    x_cyl = rs*np.cos(d2r(thetas))
    y_cyl = rs*np.sin(d2r(thetas))

    value_cyl = griddata((x.ravel(),y.ravel()), value.ravel(),  
                         (x_cyl.ravel(), y_cyl.ravel()),
                         method=method).reshape(rs.shape)

    return value_cyl

def cyl2grid(r, theta, value, x, y, method='linear'):
    u"""
    (0,0) を中心とした円筒座標系 (r,theta) から、
    格子座標 (x,y) へのスカラーの内挿を行う
    """
    x = np.asarray(x)
    y = np.asarray(y)
    value = np.asarray(value)
    
    # Check input array size
    if np.ndim(x) > 2 or np.ndim(y) > 2:
        raise ValueError, "x, y is must be 2d grid array or 1d array"
    if np.ndim(x) == 1 and np.ndim(y) ==  2:        
        x, y = np.meshgrid(x, y[:,0])
    elif np.ndim(x) == 2 and np.ndim(y) ==  1:
        x, y = np.meshgrid(x[0,:], y)
    elif np.ndim(x) == 1 and np.ndim(y) == 1:
        x, y = np.meshgrid(x, y)

    rs, thetas = np.meshgrid(r, theta)
   
    # 円筒座標の局所直交座標表現
    x_cyl = rs*np.cos(d2r(thetas))
    y_cyl = rs*np.sin(d2r(thetas))

    value_grd = griddata((x_cyl.ravel(),y_cyl.ravel()), value.ravel(),  
                         (x.ravel(), y.ravel()),
                         method=method).reshape(x.shape)

    value_grd = np.ma.array(value_grd, mask=np.isnan(value_grd))
    return value_grd

def interp_subgrid(lon, lat, data, clon, clat, typ='min', dxy=None):
    u"""
    最大・最小値をサブグリッドスケールに内挿する。
    """
    dlon = np.diff(lon)[0]
    dlat = np.diff(lat)[0]
    if (not np.all(np.diff(lon)==dlon)) or (not np.all(np.diff(lat)==dlat)): 
        raise ValueError, "grid must be equally spaced"

    if not dxy is None:
        lons, lats = np.meshgrid(lon, lat)
        mask = _domainmask(lons, lats, clon, clat, dxy, dxy)                    
        if typ == 'min':
            cidx = data[mask].argmin()
        elif typ == 'max':
            cidx = data[mask].argmax()
        clon, clat = lons[mask][cidx], lats[mask][cidx]
    print clon, clat
    cx, cy = _searchidx(clon, clat, lon, lat)

    pm = data[cy,cx]
    pa = data[cy+1,cx]
    pb = data[cy,cx-1]
    pc = data[cy-1,cx]
    pd = data[cy,cx+1]

    dy = pa - pc
    dx = pd - pb
    sy = pa + pc - 2*pm
    sx = pb + pd - 2*pm

    if typ == 'min' and min([pa,pb,pc,pd])<pm:
        raise ValueError, 'minimum data grid is on the border of searched domain.'
    elif typ == 'max' and max([pa,pb,pc,pd])>pm:
        raise ValueError, 'maximum data grid is on the border of searched domain.'

    if sx == 0:
        xfactor = 0
    else:
        xfactor = dx/sx
    if sy == 0:
        yfactor = 0
    else:
        yfactor = dy/sy
        clon  = clon - 0.5*dlon*xfactor
        clat  = clat - 0.5*dlat*yfactor
        cdata = pm - 0.125*dlon*dx*xfactor - 0.125*dlat*dy*yfactor

    return clon, clat, cdata

def _searchidx(lon, lat, lon_array, lat_array):
    u"""
    lon,latに最も近いインデックス
    """
    try:
        xidx = np.argmin(np.abs(lon_array - lon))
        yidx = np.argmin(np.abs(lat_array - lat))
    except:
        raise IndexError, "(lon,lat) = ({0},{1}) is out of bounds".format(lon,lat)
    return xidx, yidx

def _domainmask(lons, lats, clon, clat, dx, dy):
    if clon-dx < 0:
        mask = (lons<=clon+dx) | (lons>=clon+360-dx)
        mask = mask & (lats>=clat-dy) & (lats<=clat+dy)
    elif clon+dx >= 360.:
        mask = (lons>=clon-dx) | (lons<=clon-360+dx)
        mask = mask & (lats>=clat-dy) & (lats<=clat+dy)        
    else:
        mask = (lons>=clon-dx) & (lons<=clon+dx) & (lats>=clat-dy) & (lats<=clat+dy)
    return mask
