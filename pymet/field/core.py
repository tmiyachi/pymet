# coding:utf-8
import numpy as np
import copy
from datetime import datetime

__all__ = ['McGrid', 'McField', 'join']

def testgrid():
    name = 'test_grid'
    lon = np.arange(0,360,2.5)
    lat = np.arange(-90,90.1,2.5)
    lev = [1000,500,200]
    tyme = [datetime(2009,10,11), datetime(2009,10,12)]
    dims = ['tyme', 'lev', 'lat', 'lon']
    return McGrid(name=name,lon=lon,lat=lat,lev=lev,tyme=tyme,dims=dims)

def testfield():
    grid = testgrid()
    data = np.arange(73*144*3*2, dtype=np.float32).reshape(2,3,73,144)
    return McField(data, name='test_field', grid=grid)
    
    
class McGrid:
    u"""
    グリッド情報を扱うためのクラス

    :Parameters:
     **name**

     **lon**

     **lat**

     **lev**

     **tyme**

     **dims**

     **punit** : float, optional
      等圧面の気圧をPaに変換するためのパラメータ。デフォルトは100.でhPa->Paへの変換。


    **Notes**

    **Examples**
     >>> grid = pymet.McGrid(name='test_grid', lon=np.arange(0,360,2.5), lat=np.arange(-90,90.1,2.5)),
                             lev=[1000.,500.,200.], tyme=[datetime(2009,10,11),datetime(2009,10,12),
                             dims=['tyme','lev','lat','lon'])

    **Attributes**
    
    ======= ======================
    name
    lon
    lat
    lev
    tyme
    dims
    punit
    xdim
    ydim
    zdim
    tdim
    ======= ======================

    **Methods**

    .. currentmodule:: pymet.mcfield.core.McGrid
    
    ..  autosummary::
    
        copy
        gridmask
        __init__
        __getattr__
    
    """
    def __init__(self,name=None,lon=None,lat=None,lev=None,tyme=None,dims=None,punit=100.):
        self.name = name
        if hasattr(lon, '__iter__'):
            lon = np.asarray(lon)
            if len(lon)==1: lon=lon[0]
        self.lon = lon
        if hasattr(lat, '__iter__'):
            lat = np.asarray(lat)
            if len(lat)==1: lat=lat[0]
        self.lat = lat
        if hasattr(lev, '__iter__'):
            lev = np.asarray(lev)
            if len(lev)==1: lev=lev[0]
        self.lev = lev
        if hasattr(tyme, '__iter__'):
            tyme = np.asarray(tyme)
            if len(tyme)==1: tyme=tyme[0]
        self.tyme = tyme
        self.dims = dims
        self.punit = 100.
    def __getattr__(self,name):
        u"""
        xdim, ydim, zdim, tdim, edim を属性で呼び出したときに、dims属性に対応する
        次元が登録されていればその次元を返す。
        """
        try:
            if name == 'xdim':
                dname = 'lon'
                return self.dims.index('lon')
            elif name == 'ydim':
                dname = 'lat'
                return self.dims.index('lat')
            elif name == 'zdim':
                dname = 'lev'
                return self.dims.index('lev')
            elif name == 'tdim':
                dname = 'tyme'
                return self.dims.index('tyme')
            elif name == 'edim':
                dname = 'ens'
                return self.dims.index('ens')
        except ValueError:
            raise AttributeError, "McGrid instance has no dimension '{0}'".format(dname)
        raise AttributeError, "McGrid instance has no attribute '{0}'".format(name)
    def copy(self):
        u"""
        McGridオブジェクトの深いコピーを返す。
        """
        grid = McGrid(self.name)
        for a in self.__dict__:
            v = self.__dict__[a]
            if v is not None:
                grid.__dict__[a] = copy.deepcopy(v)
        return grid
    def gridmask(self, **kwargs):
        u"""
        指定した範囲の値をインデキシングするためのbool型配列を求める。

        :Arguments:
         **lon, lat, lev** : tuple or list of floats, or float, optional
           経度、緯度、鉛直次元の指定する領域。
         **tyme** : tuple or list of datetime object, or datetime object, optional
           時間次元の範囲。

        :Returns:
         **mask** : bool type ndarray

        **Examples**    
         範囲を指定する場合はタプルまたはリストで指定する。
         >>> grid = pymet.McGrid(name='test_grid', lon=np.arange(0,360,2.5), lat=np.arange(-90,90.1,2.5)),
                                  lev=[1000.,500.,200.], tyme=[datetime(2009,10,11),datetime(2009,10,12),
                                  dims=['tyme','lev','lat','lon'])
         >>> mask = grid.gridmask(lon=(0., 180.), )
        """
        for kwd in kwargs:
            if not kwd in self.dims:
                raise ValueError, "McGrid instance has no dimension {0}".format(kwd)
        mask = []
        for dim in self.dims:
            dimvalue = self.__dict__[dim]
            kwdvalue = kwargs.get(dim, None)
            if kwdvalue == None:
                mask.append(slice(None,None,None))
            elif isinstance(kwdvalue, list) or isinstance(kwdvalue, tuple):
                kwdmin, kwdmax = min(kwdvalue), max(kwdvalue)
                mask.append((kwdmin<=dimvalue) & (kwdmax>=dimvalue))
            else:
                if kwdvalue<dimvalue.min() or kwdvalue>dimvalue.max():
                    raise ValueError, "{0}={1} is out of domain".format(dim, kwdvalue)
                mask.append(np.argmin(np.abs(dimvalue-kwdvalue)))
        return mask
    
class McField(np.ma.MaskedArray):
    u"""
    格子点データを扱うためのクラス。

    :Arguments:
     **data** : array_like
      格子点データ
     **name**


    **Methods**

    .. currentmodule:: pymet.mcfield.core.McField
    
    .. autosummary::
       :toctree: generated/

       copy
       anom
       mean
        

    """
    def __new__(cls, data, **kwargs):
        print 'called  __new__'
        cls.name = None
        cls.grid = McGrid()
        return super(McField, cls).__new__(cls, data, **kwargs)
    
    def __init__(self, data, name=None, grid=None, **kwargs):
        print 'called __init__'
        super(McField,self).__init__(self, data, **kwargs)
        self.name = name
        if grid == None:
            grid = McGrid(name)
            self.grid = grid
        elif isinstance(grid,McGrid):
            self.grid = grid
        else:
            raise TypeError, "grid must be a McGrid instance"

    def copy(self):
        """
        Returns a copy of McField instance.
        """
        return McField(self.data.copy(), name=self.name,
                       grid=self.grid.copy(), mask=self.mask.copy())
                       
        
    def __getitem__(self,keys):
        print keys
        data = super(np.ma.MaskedArray,self).__getitem__(keys)
        if not isinstance(data, np.ma.MaskedArray):
            return data
        grid = self.grid.copy()
        for i, key in enumerate(keys):
            print i
            dimname = self.grid.dims[i]
            print dimname, grid.__dict__[dimname][key]
            grid.__dict__[dimname] = grid.__dict__[dimname][key]
            # スライスの結果、次元が1要素になったらその次元を削除
            if not hasattr(grid.__dict__[dimname], '__iter__'):
                grid.dims.remove(dimname)
            elif len(grid.__dict__[dimname])==1:
                grid.dims.remove(dimname)
                grid.__dict__[dimname] = grid.__dict__[dimname][0]
        return McField(data, name=self.name, grid=grid, mask=data.mask)

    # axisを引数にもつMaskedArrayのメソッドに対するラッパー
    def _axis_oper_wrapper(oper):
        def wrapper(self, *args, **kwargs):
            axis = kwargs.get('axis', None)
            grid = self.grid.copy()
            if axis!=None:
                grid.__dict__[grid.dims[axis]] = None
                grid.dims.remove(grid.dims[axis])
            result = oper(self, *args, **kwargs)
            if hasattr(result, '__iter__'):
                return McField(result, name=self.name, grid=grid, **kwargs)
            else:
                return result
        return wrapper

    # Arithmetics
    anom    = _axis_oper_wrapper(np.ma.MaskedArray.anom)
    cumprod = _axis_oper_wrapper(np.ma.MaskedArray.cumprod)
    cumsum  = _axis_oper_wrapper(np.ma.MaskedArray.cumsum)
    mean    = _axis_oper_wrapper(np.ma.MaskedArray.mean)
    prod    = _axis_oper_wrapper(np.ma.MaskedArray.prod)
    std     = _axis_oper_wrapper(np.ma.MaskedArray.std)
    sum     = _axis_oper_wrapper(np.ma.MaskedArray.sum)
    var     = _axis_oper_wrapper(np.ma.MaskedArray.var)

    #Minimum/maximum
    max = _axis_oper_wrapper(np.ma.MaskedArray.max)
    min = _axis_oper_wrapper(np.ma.MaskedArray.max)
    

    
    # ラッパーのためのクロージャー
    def _oper_wrapper(oper):
        def wrapper(self, *args, **kwargs):
            return McField(oper(self, *args, **kwargs),
                           name=self.name, grid=self.grid.copy())
        return wrapper

    __lt__        = _oper_wrapper(np.ma.MaskedArray.__lt__)
    __le__        = _oper_wrapper(np.ma.MaskedArray.__le__)
    __eq__        = _oper_wrapper(np.ma.MaskedArray.__eq__)
    __ne__        = _oper_wrapper(np.ma.MaskedArray.__ne__)
    __gt__        = _oper_wrapper(np.ma.MaskedArray.__gt__)
    __ge__        = _oper_wrapper(np.ma.MaskedArray.__ge__)
    
    __add__       = _oper_wrapper(np.ma.MaskedArray.__add__)           # +
    __sub__       = _oper_wrapper(np.ma.MaskedArray.__sub__)           # -
    __mul__       = _oper_wrapper(np.ma.MaskedArray.__sub__)           # *
    __floordiv__  = _oper_wrapper(np.ma.MaskedArray.__floordiv__)
    __mod__       = _oper_wrapper(np.ma.MaskedArray.__mod__)           #
    __divmod__    = _oper_wrapper(np.ma.MaskedArray.__divmod__)
    __pow__       = _oper_wrapper(np.ma.MaskedArray.__pow__)
    __lshift__    = _oper_wrapper(np.ma.MaskedArray.__lshift__)
    __rshift__    = _oper_wrapper(np.ma.MaskedArray.__rshift__)
    __and__       = _oper_wrapper(np.ma.MaskedArray.__and__)
    __xor__       = _oper_wrapper(np.ma.MaskedArray.__xor__)

    __div__       = _oper_wrapper(np.ma.MaskedArray.__div__)
    __truediv__   = _oper_wrapper(np.ma.MaskedArray.__truediv__)

    __radd__      = _oper_wrapper(np.ma.MaskedArray.__radd__)
    __rsub__      = _oper_wrapper(np.ma.MaskedArray.__rsub__)
    __rmul__      = _oper_wrapper(np.ma.MaskedArray.__rmul__)
    __rdiv__      = _oper_wrapper(np.ma.MaskedArray.__rdiv__)
    __rtruediv__  = _oper_wrapper(np.ma.MaskedArray.__rtruediv__)
    __rdivmod__   = _oper_wrapper(np.ma.MaskedArray.__rdivmod__)
    __rpow__      = _oper_wrapper(np.ma.MaskedArray.__rpow__)
    __rlshift__   = _oper_wrapper(np.ma.MaskedArray.__rlshift__)
    __rrshift__   = _oper_wrapper(np.ma.MaskedArray.__rrshift__)
    __rand__      = _oper_wrapper(np.ma.MaskedArray.__rand__)
    __rxor__      = _oper_wrapper(np.ma.MaskedArray.__rxor__)
    __ror__       = _oper_wrapper(np.ma.MaskedArray.__ror__)

    __iadd__      = _oper_wrapper(np.ma.MaskedArray.__iadd__)
    __isub__      = _oper_wrapper(np.ma.MaskedArray.__isub__)
    __imul__      = _oper_wrapper(np.ma.MaskedArray.__imul__)
    __idiv__      = _oper_wrapper(np.ma.MaskedArray.__idiv__)
    __itruediv__  = _oper_wrapper(np.ma.MaskedArray.__itruediv__)
    __ifloordiv__ = _oper_wrapper(np.ma.MaskedArray.__ifloordiv__)
    __imod__      = _oper_wrapper(np.ma.MaskedArray.__imod__)
    __ipow__      = _oper_wrapper(np.ma.MaskedArray.__ipow__)
    __ilshift__   = _oper_wrapper(np.ma.MaskedArray.__ilshift__)
    __irshift__   = _oper_wrapper(np.ma.MaskedArray.__irshift__)
    __iand__      = _oper_wrapper(np.ma.MaskedArray.__iand__)
    __ixor__      = _oper_wrapper(np.ma.MaskedArray.__ixor__)
    __ior__       = _oper_wrapper(np.ma.MaskedArray.__ior__)

    __neg__       = _oper_wrapper(np.ma.MaskedArray.__neg__)
    __pos__       = _oper_wrapper(np.ma.MaskedArray.__pos__)
    __abs__       = _oper_wrapper(np.ma.MaskedArray.__abs__)
    __invert__    = _oper_wrapper(np.ma.MaskedArray.__invert__)

    __int__       = _oper_wrapper(np.ma.MaskedArray.__int__)
    __long__      = _oper_wrapper(np.ma.MaskedArray.__long__)
    __float__     = _oper_wrapper(np.ma.MaskedArray.__float__)


def join(args, axis=0):
    u"""
    McFieldオブジェクトを指定した次元方向に結合する。

    :Arguments:
     **args** : tuple or list of McField objects
      結合するMcFieldオブジェクトのリストまたはタプル。結合する次元以外は同じ形状でなければならない。
     **axis** : int ot str, optional
      結合する軸を指定する。次元名('lon','lat','lev','tyme')で指定することもできる。デフォルトは0で先頭。

    :Returns:
     **out** : McField instance
      

    **Examples**
     >>> field1.grid.dims
     ['tyme', 'lev', 'lat', 'lon']
     >>> field2.grid.dims
     ['tyme', 'lev', 'lat', 'lon']
     >>> field1.shape, field2.shape
     ((2, 3, 73, 144), (4, 3, 73, 144))
     >>> field1.grid.tyme
     array([2009-10-11 00:00:00, 2009-10-12 00:00:00], dtype=object)
     >>> field2.grid.tyme
     array([2009-10-13 00:00:00, 2009-10-14 00:00:00, 2009-10-15 00:00:00,
            2009-10-16 00:00:00], dtype=object)
     >>> result = pymet.mcfield.join((field1, field2), axis='tyme')
     >>> result.shape
     (6, 3, 73, 144)
     >>> result.grid.tyme
     array([2009-10-11 00:00:00, 2009-10-12 00:00:00, 2009-10-13 00:00:00,
            2009-10-14 00:00:00, 2009-10-15 00:00:00, 2009-10-16 00:00:00], dtype=object)
            
    """
    if not isinstance(args, list) and not isinstance(args, tuple):
        raise TypeError, "input must be tuple or list of McField instance."
    if len(args)<2:
        raise ValueError, "more than two fields must be inputed."

    for arg in args:
        if not isinstance(arg, McField):
            raise TypeError, "input field must be McField instance."
        
    grid = args[0].grid.copy()
    if isinstance(axis, str):
        dim = axis
        axis = grid.dims.index(axis)
    else:
        dim = grid.dims[axis]

    data = np.ma.concatenate(args, axis=axis)
    value = np.hstack([arg.grid.__dict__[dim].copy() for arg in args])
    grid.__dict__[dim] = value
    
    return McField(data, name=grid.name, grid=grid, mask=data.mask)
