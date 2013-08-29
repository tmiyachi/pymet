# coding:utf-8
u"""
=================================================================
格子点データを扱うためのモジュール (:mod:`pymet.field`)
=================================================================

.. currentmodule:: pymet.field

------------------------
格子点情報を扱うクラス
------------------------

.. autosummary::

   McGrid
   McGrid.copy
   McGrid.latlon
   McGrid.dimindex
   McGrid.dimshape
   McGrid.gridmask
   McGrid.getgrid
   

------------------------
格子点データを扱うクラス
------------------------

.. autosummary::
   McField
   McField.get
   McField.runave
   McField.lowfreq   
   McField.mean
   McField.sum
   
-----------------------

.. autosummary::

   join

-----------------------   
pymet.gridへのラッパー
-----------------------

.. autosummary::

   dvardx
   dvardy
   dvardp
   d2vardx2
   d2vardy2
   div
   rot
   grad
   skgrad   

---------------------------   
pymet.dynamicsへのラッパー
---------------------------

.. autosummary::

   pottemp
   ertelpv
   stability
   tnflux2d
   tnflux3d

----------------------------

"""
