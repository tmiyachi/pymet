#coding:utf-8
"""
======================================================
IOモジュール (:mod:`pymet.io`)
======================================================

GrADSを用いた読み込み
----------------------

.. autosummary::

   GradsIO

netCDF形式のIO
---------------

.. autosummary::

   NetcdfIO
   NetcdfWrite

"""
import gradsio
import netcdfio
from gradsio import *
from netcdfio import *

__all__ = []
__all__ += gradsio.__all__
__all__ += netcdfio.__all__
