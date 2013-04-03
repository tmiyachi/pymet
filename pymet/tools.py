# coding: utf-8
u"""
==============================================
ユーティリティモジュール (:mod:`pymet.tools`)
==============================================

.. autosummary::

    unshape
    deunshape
    expand

--------------    
"""
import numpy as np

__all__ = ['unshape', 'deunshape', 'expand']

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

## class Timemask():
##     #mask array for sequential time series manipuration
##     def __init__(self, start, end):
##         '''
##         time sampling is days only!!
        
##         Arguments:

##             start -- start date string  ex. 1979-6-3

##             end   -- end date string
##         '''
##         startdate = parser.parse(start)
##         enddate = parser.parse(end)
##         deltadate = timedelta(days=1)
##         tnum = (enddate-startdate).days+1

##         self.timeint = numpy.array([int((startdate + deltadate*i).strftime('%Y%m%j')) for i in range(tnum)])
##         self.year = self.timeint/100000
##         self.month = self.timeint%100000/1000
##         self.days = self.timeint%1000

## class Time():
##     #mask array for sequential time series manipuration
##     def __init__(self, start, end, delta):
##         '''
##         time sampling is days only!!
        
##         Arguments:

##             start -- start date string  ex. 1979-6-3

##             end   -- end date string
##         '''
##         years, months, days, hours = 0, 0, 0, 0
##         string = ''
##         for i in delta:
##             if i == 'y' or i == 'Y':
##                 years = int(string)
##                 string = ''
##             elif i == 'm' or i == 'M':
##                 months = int(string)
##                 string = ''
##             elif i == 'd' or i == 'D':
##                 days = int(string)
##                 string = ''
##             elif i == 'h' or i == 'H':
##                 hours = int(string)
##                 string = ''
##             else:
##                 string = string + i
            
##         startdate = parser.parse(start)
##         enddate = parser.parse(end)
##         deltadate = relativedelta(years=years,months=months,days=days,hours=hours)

##         time = []
##         timeint = []
##         date = startdate
##         while(date <= enddate):
##             time.append(date)
##             timeint.append(int(date.strftime('%Y%m%j')))
##             date += deltadate

##         self.time = numpy.array(time)
##         self.timeint = numpy.array(timeint)
##         self.year = self.timeint/100000
##         self.month = self.timeint%100000/1000
##         self.days = self.timeint%1000
