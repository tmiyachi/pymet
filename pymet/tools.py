# coding: utf-8
u"""
==============================================
ユーティリティモジュール (:mod:`pymet.tools`)
==============================================

-------------------
配列を扱うツール
-------------------
.. autosummary::

    unshape
    deunshape
    expand

-------------------
文字列を扱うツール
-------------------
.. autosummary::
    d2s
    s2d
    lon2txt
    lat2txt
    
--------------    
"""
import numpy as np
from datetime import datetime
from dateutil.relativedelta import relativedelta

__all__ = ['unshape', 'deunshape', 'expand',
           'lon2txt', 'lat2txt', 'd2s', 's2d']

def unshape(a):
    u"""
    ndarray の第0軸以外を次元をつぶした 2d-array のコピーを返す。

    :Arguments:
     **a** : array_like
      2次元以上。

    :Returns:
     **array2d**  : 2darray
      2次元のndarray
     **oldshape** : tuple
      入力配列 a の形状
    
     
    **Examples**

    >>> a = np.arange(40).reshape(2,4,5)
    >>> a.shape
    (2, 4, 5)
    >>> b, oldshape = tools.unshape(a)
    >>> b.shape
    (2, 20)    
    >>> c = tools.deunshape(b, oldshape)
    >>> c.shape
    (2, 4, 5)
    """
    a  = np.array(a)
    if a.ndim < 2:
        raise ValueError, "a must be at least 2 dimension"
    
    oldshape = a.shape
    array2d  = a.reshape(oldshape[0],-1)
    
    return array2d, oldshape

def deunshape(a, oldshape):
    u"""
    unshape で reshape したndarrayを元に戻すための関数。

    :Arguments:
     **a**        : array_like
     **oldshape** : tuple

    :Returns:
     **arraynd**  : ndarray
      
    **Examples**

    >>> a = np.arange(40).reshape(2,4,5)
    >>> a.shape
    (2, 4, 5)
    >>> b, oldshape = tools.unshape(a)
    >>> b.shape
    (2, 20)
    >>> c = tools.deunshape(b, oldshape)
    >>> c.shape
    (2, 4, 5)
    """

    arraynd = np.array(a)
    arraynd = arraynd.reshape(oldshape)
    
    return arraynd

def expand(a, ndim, axis=0):
    u"""
    1次元配列にnp.newaxisを挿入して次元を拡張する。

    :Arguments:
     **a**    : array_like
      入力配列
     **ndim** : int
      出力配列の次元数
     **axis** : int, optional
      入力配列の要素を残す位置。デフォルトはゼロ(先頭)で後ろにnp.newaxisが挿入される。

    :Returns:
     **res** : ndarray
      出力配列。次元数はndim。

    **Examples**
     >>> x = np.array([1, 2, 3])
     >>> y = tools.expand(x, 3, axis=1)
     >>> y.shape
     (1, 3, 1)
     >>> y
     array([[[1],
             [2],
             [3]]])

    これは以下と同じである。
     >>> NA = np.newaxis
     >>> y = x[NA, :, NA]     
    """
    res = np.array(a)
    if res.ndim != 1:
        raise ValueError, "input array must be one dimensional array"
    idx = range(ndim)
    idx.remove(axis)
    for i in idx:
        res = np.expand_dims(res, axis=i)
    return res


__months__ = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']    
def d2s(d, fmt='%H:%MZ%d%b%Y'):
    u"""
    datetimeオブジェクトを文字列に変換する。

    :Arguments:
     **d** : datetime object
      変換するdatetimeオブジェクト
     **fmt** : str
      変換する文字列のフォーマット。書式はdatetime.strftimeに従う。localeに関わらず、%bは英語大文字の月名に変換される。
      デフォルトはGrADS形式の日付文字列'hh:mmZddmmmyyyy'
    """
    fmt = fmt.replace('%b', __months__[d.month-1])
    if d.year < 1900:
        fmt = fmt.replace('%Y', '{:04d}'.format(d.year))
        d = d + relativedelta(year=1900)
    return d.strftime(fmt)


def s2d(datestring):
    u"""
    GrADS形式の日付文字列'hh:mmZddmmmyyyy'をdatetimeオブジェクトに変換する。

    :Arguments:
     **datestring** : str
     
    """
    time, date = datestring.upper().split('Z')
    if time.count(':')>0:
        hh, mm = time.split(':')
    else:
        hh = time
        mm = 0
    dd  = date[:-7]
    mm = __months__.index(date[-7:-4])+1
    yyyy = date[-4:]
    return datetime(int(yyyy), int(mm), int(dd), int(hh), int(mm))

def lon2txt(lon):
    u"""
    経度の値を文字列に変換する。

    0度からの相対経度を東経、西経の文字列(30\N{DEGREE SIGN}E,
    140\N{DEGREE SIGN}Wなど)に変換する。

    :Argumets:
     **lon** : int
      経度(degrees)。0度からの相対経度で表す。
     
    :Returns:
     **lonlab** : str
      経度のラベル。
      
    **Examples**
     >>> lon2txt(135)
     '135\N{DEGREE SIGN}E'
     >>> lon2txt(-30)
     '30\N{DEGREE SIGN}W'
     >>> lon2txt(250)
     '110\N{DEGREE SIGN}W'
    """
    fmt = '%g'
    lon = (lon+360) % 360
    if lon>180:
        lonlabstr = u'%s\N{DEGREE SIGN}W'%fmt
        lonlab = lonlabstr%abs(lon-360)
    elif lon<180 and lon != 0:
        lonlabstr = u'%s\N{DEGREE SIGN}E'%fmt
        lonlab = lonlabstr%lon
    else:
        lonlabstr = u'%s\N{DEGREE SIGN}'%fmt
        lonlab = lonlabstr%lon
    return lonlab

def lat2txt(lat):
    u"""
    緯度の値を文字列に変換する。

    緯度を北緯、南緯の文字列(30\N{DEGREE SIGN}N,60\N{DEGREE SIGN}Sなど)
    に変換する。

    :Argumets:
     **lon** : int
      緯度(degrees)。南緯はマイナス、北緯はプラス。
     
    :Returns:
     **lonlab** : str
      緯度のラベル。
      
    **Examples**
     >>> lat2txt(60)
     '60\N{DEGREE SIGN}N'
     >>> lat2txt(-30)
     '30\N{DEGREE SIGN}S'
    """
    fmt = '%g'
    if lat<0:
        latlabstr = u'%s\N{DEGREE SIGN}S'%fmt
        latlab = latlabstr%abs(lat)
    elif lat>0:
        latlabstr = u'%s\N{DEGREE SIGN}N'%fmt
        latlab = latlabstr%lat
    else:
        latlabstr = u'%s\N{DEGREE SIGN}'%fmt
        latlab = latlabstr%lat
    return latlab
