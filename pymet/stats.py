# coding: utf-8
"""
================================================
時空間・統計解析モジュール (:mod:`pymet.stats`)
================================================
.. autosummary::

   runave
   lancoz
   regression
   rmse
   acc
   eof
   corr

---------------   
"""
import numpy as np
import scipy.signal as signal
import tools, constants
import scipy.linalg as linalg
import scipy.stats
PI = constants.pi
NA = np.newaxis

__all__ = ['runave', 'regression', 'lancoz', 
           'rmse', 'acc', 'corr',
           'eof']

def runave(a, length, axis=0, bound='mask'):
    ur"""
    移動平均。

    :Arguments:
     **a**      : array_like
      
     **length** : int
      移動平均の項数
     **axis**   : int, optional
      移動平均をする軸。デフォルトは 0。
     **bound**  : {'mask', 'valid'}, optional
      'mask':
        入力データと同じ形状の配列で結果を返す。フィルタがかからない両端はマスクされる。
      'valid':
        フィルタがかからない両端の無効な値を除いた結果を返す。入力データの時系列の長さをNとすると、
        N - length + 1. (偶数の場合は N - length) が出力の時系列の長さとなる。
    :Returns:
      **out** : ndarray or MaskedArray

    .. note::
      lengthが奇数の場合は前後length/2項、偶数の場合はlength+1項を使って先頭と末尾に0.5の重みをつける。

      すなわち、lengthがn=2k+1(奇数)のとき、

      .. math:: \bar{a}_{i} = \frac{a_{i-k} + a_{i-k+1} + \ldots + a_{i} + \ldots + a_{i+k}}{n}

      n=2k(偶数)のとき、

      .. math:: \bar{a}_{i} = \frac{0.5*a_{i-k} + a_{i-k+1} + \ldots + a_{i} + \ldots + a_{i+k-1} + 0.5*a_{i+k}}{n}

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
    nt = rdata.shape[axis]
    rdata = np.rollaxis(rdata, axis, 0)
    if nt < length:
        raise ValueError, "input array first dimension length must be larger than length."

    if ndim == 1:
        if length%2 == 1:
            w = np.ones(length)/float(length)
        else:
            w = np.r_[0.5, np.ones(length-1), 0.5]/float(length)                
        out = signal.convolve(rdata, w, mode='same')
        if bound == 'mask':
            out = np.ma.asarray(out)
            out[:length/2] = np.ma.masked
            out[-(length/2):] = np.ma.masked
        elif bound == 'valid':
            out = out[(length/2):-(length/2)]
        else:
            raise ValueError, "unexpected bound option '%{0}'".format(bound)
    else:
        rdata, oldshape = tools.unshape(rdata) 
        if length%2 != 0:
            w = np.c_[(np.zeros(length),np.ones(length),np.zeros(length))]/float(length)
        else:
            w = np.c_[(np.zeros(length+1),np.r_[0.5, np.ones(length-1), 0.5],
                           np.zeros(length+1))]/float(length)
        out = signal.convolve2d(rdata, w, mode='same')
        if bound=='mask':
            out = np.ma.asarray(out)
            out[:length/2,:] = np.ma.masked
            out[-(length/2):,:] = np.ma.masked
        elif 'valid':
            out = out[(length/2):-(length/2),:]
        else:
            raise ValueError, "unexpected boud option '%{0}'".format(bound)
        out = tools.deunshape(out, (-1,) + oldshape[1:])
        out = tools.mrollaxis(out, 0, axis+1)
        
    return out

def lancoz(a, cutoff, cutoff2=None, length=None, axis=0, bound='mask', mode='lowpass'):
    ur"""
    時系列データに対してLancozフィルタをかける。

    :Arguments:
     **a** : array_like
      入力データ
     **cutoff** : int
      カットオフ(イン)周期。入力データのステップ数で指定
     **cutoff2** : int
      バンドパスフィルタをかける場合の2つ目のカットオフ(イン)周期。入力データのステップ数で指定
     **length** : int
      項数。デフォルトは 2*max(cutoff,cutoff2)+1
     **axis**   : int, optional
      移動平均をする軸。デフォルトは 0。
     **bound**  : {'mask', 'valid'}, optional
      'mask':
        入力データと同じ形状の配列で結果を返す。フィルタがかからない両端はマスクされる。
      'valid':
        フィルタがかからない両端の無効な値を除いた、(入力データの時間次元の長さ) - length + 1.
        の配列で結果を返す。
     **mode** : {'lowpass', 'highpass', 'bandpass'}, optional
      'lowpass':
       low-passフィルタとして施す。デフォルト。
      'hoghpass'
       high-passフィルタとして施す。
      'bandpass'
       band-passフィルタとして施す。
    :Return:
      **out** : array_like

    .. note::
     カットオフ周波数fc、項数nのLancoz低周波フィルタの重み関数は、

     .. math:: w_{k} = 2f_{c}\,{\rm sinc}(2f_{c}k)\,{\rm sinc}(k/n) \hspace{3em} k = -n, \ldots, n

     で表される。ここでsinc(x)は正規化されたsinc関数で、

     .. math:: {\rm sinc}(x) = \frac{\sin(\pi x)}{\pi x} \hspace{3em} {\rm sinc}(0) = 1

    **Examples**
     サンプリング間隔1日のデータに、カットオフ周期10日、項数21のLancoz低周波フィルタをかける。
      >>> data.shape
      (100, 21, 73, 144)
      >>> lowpass = lancoz(data, 10, axis=0, length=21, bound='same', mode='lowpass')
      
    **Refferences**
      Duchon, C. E., 1979: Lanczos Filtering in One and Two Dimensions. J. Appl. Meteor., 18, 1016–1022.
      doi: <http://dx.doi.org/10.1175/1520-0450(1979)018<1016:LFIOAT>2.0.CO;2>
       
    """
    rdata = np.array(a)
    ndim = rdata.ndim
    nt = rdata.shape[0]
    rdata = np.rollaxis(rdata, axis, 0)
    
    if length == None: length = 2*max(cutoff,cutoff2) + 1
    if nt < length:
        raise ValueError, "input array time dimension length must be larger than 2*cutoff+1."
    elif length%2 == 0:
        raise ValueError, "length keyword must be odd number"

    fc = 1./cutoff  # cutoff(or in) frequency
    n = length/2
    k = np.arange(-n, n+1)
    if mode == 'lowpass':
        w = 2. * fc * np.sinc(2.*fc*k) * np.sinc(k/n)
    elif mode == 'highpass':
        w = -2. * fc * np.sinc(2.*fc*k) * np.sinc(k/n)
    elif mode == 'bandpass':
        if cutoff2 == None: raise ValueError, "cutoff2 value is required in bandpass mode"
        fc1 = max(fc, 1./cutoff2) # cut off
        fc2 = min(fc, 1./cutoff2) # cut in
        w = 2. * (fc1*np.sinc(2.*fc1*k) - fc2*np.sinc(2.*fc2*k)) * np.sinc(k/n)
    else:
        raise ValueError, "mode '{0}' is invalid".format(mode)
            
    if ndim == 1:
        out = signal.convolve(rdata, w, mode='same')
        if bound == 'mask':
            out = np.ma.asarray(out)
            out[:length/2] = np.ma.masked
            out[-(length/2):] = np.ma.masked
        elif bound == 'valid':
            out = out[(length/2):-(length/2)]
        else:
            raise ValueError, "unexpected bound option '%{0}'".format(bound)
    else:
        rdata, oldshape = tools.unshape(rdata) 
        w = np.c_[(np.zeros(length), w, np.zeros(length))]
        out = signal.convolve2d(rdata, w, mode='same')
        if bound=='mask':
            out = np.ma.asarray(out)
            out[:length/2,:] = np.ma.masked
            out[-(length/2):,:] = np.ma.masked
        elif 'valid':
            out = out[(length/2):-(length/2),:]
        else:
            raise ValueError, "unexpected boud option '%{0}'".format(bound)
        out = tools.deunshape(out, (-1,) + oldshape[1:])
        out = tools.mrollaxis(out, 0, axis+1)

    return out

def runcum(a, length, axis=0, bound='mask'):
    ur"""
    移動積算値を計算する。

    :Arguments:
     **a**      : array_like
      入力配列
     **length** : int
      積算する項数
     **axis**   : int, optional
      積算する軸。デフォルトは 0。
     **bound**  : {'mask', 'valid'}, optional
      'mask':
        入力データと同じ形状の配列で結果を返す。積算する項が足りない先端はマスクされる。
      'valid':
        積算する項が足りない先端を除いた長さの結果を返す。入力データの時系列の長さをNとすると、
        N - length + 1. が出力の時系列の長さとなる。
    :Returns:
      **out** : ndarray or MaskedArray

    **Examples**
     >>> a = np.arnge(10)
     >>> a
     array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
     >>> runcum(a, 3, bound='mask')
     masked_array(data = [-- -- 3.0 6.0 9.0 12.0 15.0 18.0 21.0 24.0],
                  mask = [ True  True False False False False False False False False],
                         fill_value = 1e+20)
    >>> runcum(a, d, mode='valid')
    array([  3.,   6.,   9.,  12.,  15.,  18.,  21.,  24.])
    """
    rdata = np.array(a)
    ndim = rdata.ndim
    nt = rdata.shape[axis]
    rdata = np.rollaxis(rdata, axis, 0)
    if nt < length:
        raise ValueError, "input array first dimension length must be larger than length."

    if ndim == 1:
        w = np.ones(length)
        out = signal.convolve(rdata, w, mode='full')
        if bound == 'mask':
            out = out[:-length+1]
            out = np.ma.asarray(out)
            out[:length-1] = np.ma.masked
        elif bound == 'valid':
            out = out[length-1:-length+1]
        else:
            raise ValueError, "unexpected bound option '%{0}'".format(bound)
    else:
        rdata, oldshape = tools.unshape(rdata) 
        w = np.c_[(np.zeros(length),np.ones(length),np.zeros(length))]
        out = signal.convolve2d(rdata, w, mode='full')
        if bound=='mask':
            out = out[:-length+1,1:-1]
            out = np.ma.asarray(out)
            out[:length-1,:] = np.ma.masked            
        elif 'valid':
            out = out[length-1:-length+1,1:-1]
        else:
            raise ValueError, "unexpected boud option '%{0}'".format(bound)
        out = tools.deunshape(out, (-1,) + oldshape[1:])
        out = tools.mrollaxis(out, 0, axis+1)
        
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
        prob = 1-scipy.stats.t.sf(np.abs(t),dof-2)*2        
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
        prob = 1-scipy.stats.t.sf(np.abs(t),dof-2)*2
    
        # 元の形状に戻す
        a = tools.deunshape(a, oldshape)
        b = tools.deunshape(b, oldshape)
        prob = tools.deunshape(prob, oldshape)
        a = np.rollaxis(a, 0, axis+1)
        b = np.rollaxis(b, 0, axis+1)
        prob = np.rollaxis(prob, 0, axis+1)
        
    return a, b, prob

def corr(a, b, axis=None):
    u"""
    相関係数を計算する。

    :Arguments:
     **a, b** : array_like
      入力
     **axis** : int, optional
      a, b の一方が1次元で、もう一方が多次元の場合に相関を計算する軸。
    :Returns:
     **r** : array_like
           
    """
    a = np.asarray(a)
    b = np.asarray(b)

    if axis == None:
        if not a.shape == b.shape:
            raise ValueError, "input array must have shame shape"
    else:
        if a.ndim == 1 and b.ndim>1:
            a = tools.expand(a, b.ndim, axis=axis)            
        elif b.ndim == 1 and a.ndim>1:
            b = tools.expand(b, a.ndim, axis=axis)
    
    a_anom = a - a.mean(axis=axis)
    b_anom = b - b.mean(axis=axis)

    xy = (a_anom*b_anom).sum(axis=axis)
    xx = np.sqrt((a_anom**2).sum(axis=axis))
    yy = np.sqrt((b_anom**2).sum(axis=axis))

    r = xy/xx/yy
    
    return r

def rmse(var, basis, axes=None):
    u"""
    二乗平均誤差(Root Mean Square Error)を計算する。

    :Arguments:
     **var** : ndarray
      予報値
     **basis** : ndarray
      基準となる値
    :Returns:
     **out** : float
      二乗平均誤差
    """
    ## if not exaxes==None:
    ##     for axis in list(exaxes).sort():
    ##         basis = np.expand_dims(basis, axis=axis)
    var = np.ma.asarray(var)
    basis = np.ma.asarray(basis)
    if var.ndim!=basis.ndim:
        raise ValueError, "input array size is incorrect"
    out = (var - basis)**2
    if axes==None:
        out = np.ma.sqrt(np.ma.mean(out))
    elif isinstance(axes, int):
        out = np.ma.sqrt(np.ma.mean(out,axis=axes))
    else:
        axes = list(axes).sort()
        for i, axis in enumerate(axes):
            out = np.ma.mean(axis=axis-i)
        out = np.ma.sqrt(out)
    return out
        
def acc(fcst, anal, clim):
    u"""
    アノマリー相関を計算する。

    :Arguments:
     **fcst** : ndarray
      予報値
     **anal** : ndarray
      解析値
     **clim** : ndarray
      気候値
    :Returns:
     **out**  : float
      アノマリー相関。

    **Examples**
    """
    

    
    if not fcst.ndim==anal.ndim==clim.ndim:
        raise ValueError, "input array size is incorrect"
    return np.corrcoef((fcst-clim),(anal-clim))[0,1]


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

    
