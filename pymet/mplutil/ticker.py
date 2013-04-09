# coding:utf-8
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.offsetbox import AnchoredText
from matplotlib.patheffects import withStroke
import matplotlib.colors as colors

__all__ = ['zeroclevs', 'zeroclevs_multi']

def zeroclevs(maxval,hnum=10,skip=2):
    u"""
    0を中心としたコンターレベルの配列を計算する。

    :Arguments:
     **args** : float
      求めたい配列の最大値。
     **hnum** : int, optional
      等分数。
     **skip** : int, optional

    :Returns:
     **levs** : array
      2*hnum+1要素の配列。
     **ticks** : array
    """
    if maxval<=0:
        raise ValueError, "maxval must be positive"
    levs = np.delete(np.linspace(-maxval,maxval,2*hnum+1),hnum)
    tick_half = np.linspace(0,maxval,hnum+1)[1::skip]
    tick = np.r_[-tick_half[::-1],tick_half]
    return levs, tick
        
def zeroclevs_multi(*args,**kwargs):
    u"""
    0を中心としたコンターレベルの配列を計算する。

    :Arguments:
     **data** : array_like
      コンターレベルの決定に使う全データ。複数与えることができる。
     **skip** : int, optional
      ticksラベル配列の間引き間隔。デフォルトは2。
     **factor** : float, optional
      自動決定されたコンターレベルにかけるファクター。

    :Returns:
     **levs** : array
      コンターレベルの配列。0を中心に正負各10等分した長さ21の配列。
     **ticks** : array
      levsを0を含むようにskipだけ間引いたticksラベルのための配列。

    .. notes::
     コンターラベルは入力した配列の絶対値での最大値から適切な間隔を決定する。
    """
    skip = kwargs.get('skip',2)
    factor = kwargs.get('factor',1)
    maxval = 0.
    for arg in args:
        maxval = max(np.abs(arg).max(), maxval)
    decimal = np.ceil(np.log10(maxval)-1)
    maxhead = np.ceil(maxval/10**decimal)
    levs = np.delete(np.linspace(-maxhead,maxhead,21),10)*10**decimal * factor
    tick_half = np.linspace(0,maxhead,11)[1::skip]*10**decimal * factor
    tick = np.r_[-tick_half[::-1],tick_half]
    return levs, tick
       
