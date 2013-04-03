# coding: utf-8
"""
================================================
時空間・統計解析モジュール (:mod:`pymet.stats`)
================================================
.. autosummary::

   runave
   regression
   eof

---------------   
"""
import numpy as np
import scipy.signal as signal
import tools, constants
import scipy.linalg as linalg
import scipy.stats as stats
PI = constants.pi
NA = np.newaxis

__all__ = ['runave', 'regression', 'eof']

def runave(a, length, axis=0, bound='mask'):
    u"""
    移動平均。

    :Arguments:
     **a**      : array_like
      
     **length** : int
      移動平均の項数
     **axis**   : int, optional
      移動平均をする軸。デフォルトは 0。
     **bound**  : {'mask', 'valid'}, optional
      'mask':
        By default, mode is same. This returns output of first dimension length N.
        Boundary componets are masked.
      'valid':
       Mode valid returns output of length N - length + 1.
       running average is only given for points where the signals overlap completely.

    :Returns:
      **out** : ndarray or MaskedArray


    **Notes**

    lengthが奇数の場合は前後length/2項、偶数の場合はlength+1項を使って先頭と末尾に0.5の重みをつける。
    すなわち、length=5のときは、out[i] = (a[i-2]+a[i-1]+a[i]+a[i+1]+a[i+2])/5. 、    
    length=4のときは、out[i] = (0.5*a[i-2]+a[i-1]+a[i]+a[i+1]+0.5*a[i+2])/5. となる。    

    
    **Examples**
       
     >>> runave([1, 2, 3, 4, 5], 3, bound='mask')
     masked_array(data = [-- 2.0 3.0 4.0 --], 
                  mask = [ True False False False True ],
            fill_value = 1e+20)

     >>> runave([1, 2, 3, 4, 5], 3, bound='valid')
     array([ 2.  3.  4.])
    """
    rdata = np.array(a)
    ndim = rdata.ndim
    rdata = np.rollaxis(rdata, axis, ndim)
    ntim = rdata[0]
    if ntim < length:
        raise ValueError, "input array first dimension length must be larger than length."

    if ndim == 1:
        if length%2 != 0:
            w = np.ones(length)/float(length)
        else:
            w = np.r_[0.5, np.ones(length-1), 0.5]/float(length)                
        out = signal.convolve(rdata, w, mode='same')
        if bound=='mask':
            out = np.ma.asarray(out)
            out[:length/2] = np.ma.masked
            out[-(length/2):] = np.ma.masked
        elif 'valid':
            out = out[(length/2):-(length/2)]
        else:
            raise ValueError, "unexpected boud option '%{0}'".format(bound)
    else:
        rdata, oldshape = tools.unshape(rdata) 
        if length%2 != 0:
            w = np.vstack((np.zeros(length),np.ones(length),np.zeros(length)))/float(length)
        else:
            w = np.vstack((np.zeros(length+1),np.r_[0.5, np.ones(length-1), 0.5],\
                           np.zeros(length+1)))/float(length)
        out = signal.convolve2d(rdata, w, mode='same')
        if bound=='mask':
            out = np.ma.asarray(out)
            out[:length/2,:] = np.ma.masked
            out[-(length/2):,:] = np.ma.masked
        elif 'valid':
            out = out[(length/2):-(length/2),:]
        else:
            raise ValueError, "unexpected boud option '%{0}'".format(bound)
        out = tools.deunshape(out, oldshape)
        out = np.rollaxis(out, ndim-1, axis)
        
    return out

def lancoz(rdata,cutoff,mode='same'):
    """calculate low-pass filter (Duchon, 1979)    

    Arguments:

       'array'  -- (N,...) first dimension must be time. 

       'cutoff' -- cutoff frquency 

       'mode'   -- {'valid', 'same'}, optional
              'same':
                    By default, mode is same. This returns output of first dimension length N.
                    Boundary componets are masked.
              'valid':
                    Mode valid returns output of length N - 2*cutoff.
                    running average is only given for points where the signals overlap completely. 

    Return:
       'out'    -- filtered array

    """
    ntim = rdata.shape[0]
    if ntim < 2*cutoff+1:
        print "input array first dimension length must be larger than 2*fl+1."
        sys.exit()
    n = cutoff
    i = np.arange(-n,n+1)
    w = np.sin(2.*PI*cutoff*i)*np.sin(PI*i/n)/PI/i/PI/pi*n
    w = w/w.sum()

    if rdata.ndim == 1:
        out = np.ma.array(signal.convolve(rdata, w, mode=mode))
        if mode=='same':
            out[:n] = np.ma.masked
            out[-n:] = np.ma.masked
    elif rdata.ndim == 2:
        w = np.vstack([np.zeros(len(w)),w,np.zeros(len(w))])
        out = np.ma.array(signal.convolve2d(rdata, w, mode=mode))
        if mode=='same':
            out[:n,:] = np.ma.masked
            out[-n:,:] = np.ma.masked
    else :
        rdata, oldshape = unshape(rdata) 
        w = np.vstack([np.zeros(len(w)),np.ones(length),np.zeros(len(w))])
        out = np.ma.array(signal.convolve2d(rdata, w, mode='same'))
        if mode=='same':
            out[:n,:] = np.ma.masked
            out[-n:,:] = np.ma.masked
        else:
            out = out[n:-n,:]
        out = out.reshape((-1,)+oldshape[1:])

    return out

def regression(x, y, axis=0, dof=None):
    u"""
    線形回帰式を求める。

    :Arguments:
     **x**    : array_like

     **y**    : ndarray

     **axis** : int, optional
       
     **dof** : int, optional
       自由度。

    :Returns:
     **a, b** : floats or ndarray
       線形回帰式 y = a*x + b の係数 a, b。yと同じ形状になる。
     **prob** : floats or ndarray
       t-検定における??????

    **Notes**

    **Examples**    
     >>> x = 
    """
    x = np.asarray(x)
    y = np.asarray(y)
    ndim = y.ndim
    if x.ndim > 1:
        raise ValueError, "error x must be 1-D array"
    if dof == None: dof=len(x)
        
    if ndim == 1:
        # 1次線形回帰
        p = np.polyfit(x, y, 1)
        a = p[0]
        b = p[1]
        
        # 自由度dofのt-検定
        e = y - a*x - b
        t = a*x.var()/np.sqrt(e.var()/(dof-2.))
        prob = 1-stats.t.sf(np.abs(t),dof-2)*2        
    else:
        # 回帰する軸を先頭に
        y = np.rollaxis(y, axis, 0)
        y, oldshape = tools.unshape(y)

        # 1次線形回帰
        p = np.polyfit(x, y, 1)
        a = p[0]
        b = p[1]

        # t-検定
        e = y - x[:,NA]*a[NA,:] - b[NA,:]
        t = a*x.var()/np.sqrt(e.var(axis=0)/(dof-2.))
        prob = 1-stats.t.sf(np.abs(t),dof-2)*2
    
        # 元の形状に戻す
        a = tools.deunshape(a, oldshape)
        b = tools.deunshape(b, oldshape)
        prob = tools.deunshape(prob, oldshape)
        a = np.rollaxis(a, 0, axis+1)
        b = np.rollaxis(b, 0, axis+1)
        prob = np.rollaxis(prob, 0, axis+1)
        
    return a, b, prob

def eofold(data,lat=None,ano=True):
    """ 
    EOF解析。

    :Arguments: 
      **data** : ndarray

      **tdim** : int, optional
      **lat**  : array_like, optional

      **ydim** : int, optional

    :Returns:
     **EOFs**
      mth eof mode
     **PCs**
       mth of pc
     **explained**
       explains of mth mode (%)

    **Notes**
       
    **Referrences**
     
    **Examples**
     >>>
     
    """
    print data.shape
    if data.ndim == 2:
        ndim = 2
        tn, xn = data.shape
    elif data.ndim == 3:
        ndim = 3
        tn, yn, xn = data.shape
    else:
        raise ValueError, 'array must be 2D or 3D.'
        
    #calculate anomaly
    if not ano:
        data = data - data.mean(axis=0)
    # weight cos(lat)
    if lat != None:
        factor = np.sqrt(np.cos(PI/180.*lat))
        data = data * factor[NA,:,NA]
    if ndim == 3:
        data = data.reshape(tn,xn*yn)

    #force double precision
    dtype = data.dtype
    if dtype == np.float32:
        data = np.float64(data)
                
    #SVD
    A, Lh, E = linalg.svd(data, full_matrices=False)

    lambdas = Lh*Lh/(len(data)-1)
    EOFs = np.transpose(E)
    PCs = A*Lh

    #renormalized
    EOFs = EOFs * np.sqrt(lambdas)[NA,:]
    PCs = PCs / np.sqrt(lambdas)[NA,:]

    #reshape
    if ndim == 3:
        EOFs = EOFs.reshape(yn,xn,-1)
    #cos(lat)
    if lat != None:
        EOFs = EOFs / factor[:,NA,NA]
    #kiyoritu
    explained = lambdas / np.add.reduce(lambdas) *100.
    EOFs = np.rollaxis(EOFs,-1,0)
    PCs = np.rollaxis(PCs,-1,0)    

    if dtype == np.float32:
        explained = np.float32(explained)
        EOFs = np.float32(EOFs)
        PCs = np.float32(PCs)
    
    return EOFs, PCs, explained

def eof(data, tdim=0, lat=None, ydim=None):
    """ 
    EOF解析。

    :Arguments: 
      **data** : ndarray
       入力する2次元以上のデータ配列。
      **tdim** : int, optional
       入力データの時間次元の軸。デフォルトは0、すなわち先頭。
      **lat**  : array_like, optional
       指定すると緯度に応じた面積重みをつける。       
      **ydim** : int, optional
       latを指定した場合の緯度次元の軸
       
    :Returns:
     **EOFs** : ndarray
       m番目のEOFモードの空間構造。形状(M,Xn),Xnは空間方向のデータ数、M=min(Xn,Tn)
     **PCs** : ndarray
       m番目のEOFモードの時系列。形状(M,Tn),Tnは時間方向のデータ数。M=min(Xn,Tn)
     **lambdas** : 1darray
       m番目のEOFモードの寄与率[%]。長さM。

    .. note::
     計算の際にはSVD解析を用いる。
       
    **Referrences**
     
    **Examples**
     >>>
     
    """    
    data = np.asarray(data)
    ndim = data.ndim
    dtype = data.dtype
    if ndim<2:
        raise ValueError, "input data must have more than two dimendin."

    #緯度重みを考慮
    if lat!=None and ydim!=None:
        factor = tools.expand(np.sqrt(np.cos(PI/180.*lat)), ndim, axis=ydim)
        data = data * factor
        
    #時間軸を先頭に
    data = np.rollaxis(data, tdim, 0)
    data, oldshape = tools.unshape(data)
    tn, xn = data.shape
        
    #SVD
    A, Lh, E = linalg.svd(data, full_matrices=False)

    lambdas = Lh*Lh/tn

    EOFs = E
    PCs = np.transpose(A*Lh)
    
    #規格化
    EOFs = EOFs * np.sqrt(lambdas)[:,NA]
    PCs = PCs / np.sqrt(lambdas)[:,NA]

    #寄与率
    lambdas = lambdas / lambdas.sum() * 100.
    
    # 形状を戻す
    eofshape = list(oldshape)
    eofshape.pop(tdim)
    eofshape.insert(0, len(lambdas))
    EOFs = EOFs.reshape(eofshape)           #(M,Tn,...)
    print eofshape
    #緯度重みを考慮
    if lat!=None and ydim!=None:
        factor = tools.expand(np.sqrt(np.cos(PI/180.*lat)), ndim, axis=ydim)
        EOFs = EOFs / factor

    return EOFs, PCs, lambdas

def ceof(data,tdim=0):
    """
    複素EOF解析。

    :Arguments:
     **data** : ndarray

     **tdim** : int, optional


    :Returns:



    .. note::


    **Examples**
     >>>
     
    """
    data = np.asarray(data)
    ndim = data.ndim
    if ndim<2:
        raise ValueError, "input data have more than two dimension."

    #データ行列X(t,x)の作成
    data = np.rollaxis(data, tdim, ndim)
    X, oldshape = tools.unshape(data)
    tn, xn = X.shape
    
    #精度向上のためdouble型にする
    
    #Hilbert変換
    Z = signal.hilbert(X, axis=0)

    #共分散行列の計算
    C = signal.convolve
    
    out = tools.deunshape(d, oldshape)
    np.rollaxis(out, tdim-1, ndim)

    
