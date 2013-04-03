# coding:utf-8
u"""
=================================================================
プロットのためのモジュール (:mod:`pymet.plot`)
=================================================================

地図描画
--------

.. autosummary::
   MyBasemap

   MyBasemap.fillcontinents
   MyBasemap.plot
   MyBasemap.contour
   MyBasemap.contourf
   MyBasemap.quiver
   MyBasemap.colorbar
   MyBasemap.quiverkey

   MyBasemap.set_xlint
   MyBasemap.set_ylint
   
緯度経度ラベル
---------------

.. autosummary::
   BasemapXaxisLocator
   BasemapYaxisLocator
   BasemapXaxisFormatter   
   BasemapYaxisFormatter
   lon2txt
   lat2txt
   
matplotlib
-----------
   
-----------------------
"""

import mybasemap
from mybasemap import *

__all__ = []
__all__ += mybasemap.__all__

