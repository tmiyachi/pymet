# coding:utf-8
import numpy as np
from scipy.interpolate import RectBivariateSpline, griddata

def interp_center(x, y, data, cx, cy, rmax_domain, rint=None, mode='linear'):
    u"""
    (cx,cy)を中心とした相対座標に内挿する
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
    
