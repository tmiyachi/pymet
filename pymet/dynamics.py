# coding: utf-8
u"""
=========================================================
物理量計算モジュール (:mod:`pymet.dynamics`)
=========================================================

.. autosummary::

   pottemp
   ertelpv
   stability
   tnflux2d
   tnflux3d

-------------------   
"""
import numpy as np
import constants, tools
from grid import *

__all__ = ['pottemp', 'ertelpv', 'tnflux2d', 'tnflux3d', 'stability']

NA=np.newaxis
kappa = constants.air_kappa
rd = constants.air_rd
g = constants.earth_gravity

def pottemp(temp, lev, zdim, punit=100., p0=100000.):
    ur"""
    温位を計算する。

    :Arguments:
     **temp**  : array_like
      気温 [K]
     **lev**   : array_like
      気圧 [punit*Pa]
     **zdim**  : int
      鉛直次元の軸
     **punit** : float, optional

     **p0**    : float, optional
      基準気圧。デフォルトは1000hPa。
      
    :Returns:
      **out** : ndarray

    .. note::
       温位は次のように定義される。

       .. math:: \theta = T\left( \frac{p_{0}}{p}  \right)^{\frac{R_d}{C_p}}

    **Referrences**
    
    
    **Examples**
    
    >>>
    """    
    temp = np.asrray(temp)
    ndim = temp.ndim
    p = tools.expand(lev, ndim, axis=zdim)

    out = temp * ((p0/lev/unit)**kappa)

    return out
    
def ertelpv(uwnd, vwnd, temp, lon, lat, lev, xdim, ydim, zdim, cyclic=True, punit=100.):
    ur"""
    エルテルのポテンシャル渦度を計算する。

    :Arguments:
     **uwnd, vwnd** : ndarray
      東西風、南北風 [m/s]
     **temp**       : ndarray
      気温 [K]
     **lon, lat** : array_like
      経度、緯度 [degrees]
     **lev**   : array_like
      気圧 [punit*Pa]
     **xdim, ydim, zdim** : int
       東西、南北、鉛直次元の軸
     **cyclic** : bool, optional
      東西境界を周期境界で扱うか。デフォルトはTrue。
     **punit**  : float, optional
      Paに換算するための定数。デフォルトは100、すなわちlevはhPaで与えられると解釈する。
      
    :Returns:
      **out** : ndarray

    .. note::
    
       p座標系のポテンシャル渦度は次のように定義される。
      
       .. math:: P = -g (f\mathbf{k}+\nabla_{p}\times\mathbf{v})\cdot\nabla_{p}\theta      
       鉛直風を含む項は他の項に比べ小さいから、次の近似的値を計算している。

       .. math:: P \simeq -g\left\{ (f + \zeta)\frac{\partial\theta}{\partial p} - \frac{\partial v}{\partial p} \frac{\partial\theta}{\partial x} + \frac{\partial u}{\partial p}\frac{\partial\theta}{\partial y}\right\}
      
    **Referrences**
      Hoskins, B.J., M.E. McIntyre and A.E. Robertson, 1985: On the use and significance of isentropic potential vorticity maps, `QJRMS`, 111, 877-946, <http://onlinelibrary.wiley.com/doi/10.1002/qj.49711147002/abstract> 
    
    **Examples**   
     >>>
    """
    u, v, t = np.asarray(uwnd), np.asarray(vwnd), np.asarray(temp)
    ndim=u.ndim

    # potential temperature
    theta = pottemp(t, lev, zdim, punit=punit)
    #
    dthdp = dvardp(theta, lev, zdim, punit=punit)
    dudp  = dvardp(u, lev, zdim, punit=punit)
    dvdp  = dvardp(v, lev, zdim, punit=punit)

    dthdx = dvardx(theta, lon, lat, xdim, ydim, cyclic=cyclic)
    dthdy = dvardy(theta, lat, ydim)

    # absolute vorticity
    vor  = rot(u, v, lon, lat, xdim, ydim, cyclic=cyclic)
    f    = tools.expand(constants.const_f(lat), ndim, axis=ydim)
    avor = f + vor

    out = -g * (avor*dthdp - (dthdx*dvdp-dthdy*dudp))

    return out

def stability(temp, lev, zdim, punit=100.):
    ur"""
    p座標系での静的安定度(Brunt-Vaisala振動数)を計算する。

    :Arguments:
     **temp**  : array_like
       気温 [K]
     **lev**   : array_like
       等圧面の気圧
     **zdim**  : int
       鉛直次元の位置
     **punit** : float, optional
       Paに変換するためのファクター。デフォルトは100.。
    :Returns:
     **N** : ndarray
      静的安定度

    .. note::
       p面上の静的安定度は次式で定義される。
     
       .. math:: N^2 = g\frac{\partial (\ln\theta)}{\partial z} = -\alpha\frac{\partial(\ln\theta)}{\partial p}

    **References**    

    **Examples**
     >>>
    
    """
    temp = np.asarray(temp)
    ndim = temp.ndim
    p = tools.expand(lev, ndim, axis=zdim)*punit
    theta = pottemp(temp, lev, zdim, punit=punit)
    alpha = rd*temp/p
    N = -alpha * dvardp(np.log(theta), lev, zdim, punit=punit)

    return N

def tnflux2d(U, V, strm, lon, lat, xdim, ydim, cyclic=True, limit=100):
    ur"""
    Takaya & Nakamura (2001) の波活動度フラックスの水平成分をp面上で計算する。

    :Arguments:
     **U, V** : ndarray
      気候値東西風、南北風 [m/s]
     **strm** : ndarray
      流線関数 [m^2/s]
     **lon, lat** : array_like
      経度、緯度 [degrees]
     **xdim, ydim** : int
       東西、南北次元の軸
     **cyclic** : bool, optional
      東西境界を周期境界で扱うか。デフォルトはTrue。
     **limit** : float, optional
      デフォルトは100.
     
    :Returns:
     **tnx, tny** : MaskedArray
      フラックスの東西、南北成分
     
    .. note::
       p座標系のwave-activityフラックスのx,y成分は次で定義される。

       .. math:: \mathbf{W}=\frac{1}{2|\mathbf{U}|}
          \left( \begin{array}{c}
          U({\Psi'}_{x}^{2}-\Psi'{\Psi'}_{xx})        + V({\Psi'}_{x}{\Psi'}_{y}-\Psi'{\Psi'}_{xy})\\
          U({\Psi'}_{x}{\Psi'}_{y}-\Psi'{\Psi'}_{xy}) + V({\Psi'}_{x}^{2}-\Psi'{\Psi'}_{xx})
          \end{array} \right)

    **References**
      Takaya, K and H. Nakamura, 2001: A formulation of a phase-independent wave-activity flux for stationary and migratory quasigeostrophic eddies on a zonally varying basic flow, `JAS`, 58, 608-627.
      
      <http://journals.ametsoc.org/doi/abs/10.1175/1520-0469%282001%29058%3C0608%3AAFOAPI%3E2.0.CO%3B2>

    **Examples**
     >>>
     
    """
    U, V = np.asarray(U), np.asarray(V)
    ndim=U.ndim

    dstrmdx    = dvardx(strm, lon, lat, xdim, ydim, cyclic=cyclic)
    dstrmdy    = dvardy(strm, lat, ydim)
    d2strmdx2  = d2vardx2(strm, lon, lat, xdim, ydim, cyclic=cyclic)
    d2strmdy2  = d2vardy2(strm, lat, ydim)
    d2strmdxdy = dvardy(dvardx(strm,lon,lat,xdim,ydim,cyclic=cyclic),lat,ydim)

    tnx = U * (dstrmdx**2 - strm*d2strmdx2) + V*(dstrmdx*dstrmdy - strm*d2strmdxdy)
    tny = U * (dstrmdx*dstrmdy - strm*d2strmdxdy) + V * (dstrmdy**2 - strm*d2strmdy2)

    tnx = 0.5*tnx/np.abs(U + 1j*V)
    tny = 0.5*tny/np.abs(U + 1j*V)

    # 大きさがlimit以上のグリッドと東風領域をマスクする。
    tnxy = np.sqrt(tnx**2 + tny**2)
    tnx = np.ma.asarray(tnx)
    tny = np.ma.asarray(tny)
    tnx[tnxy>limit] = np.ma.masked
    tny[tnxy>limit] = np.ma.masked
    tnx[U<0] = np.ma.masked
    tny[U<0] = np.ma.masked
    
    return tnx, tny

def tnflux3d(U, V, T, strm, lon, lat, lev, xdim, ydim, zdim, cyclic=True, limit=100, punit=100.):
    ur"""
    Takaya & Nakamura (2001) の波活動度フラックスをp面上で計算する。

    :Arguments:
     **U, V** : ndarray
      気候値東西風、南北風 [m/s]
     **T**    : ndarray
      気候値気温 [K]
     **strm** : ndarray
      流線関数偏差 [m^2/s]
     **lon, lat** : array_like
      経度、緯度 [degrees]
     **lev** : array_like
      等圧面の気圧
     **xdim, ydim, zdim** : int
      東西、南北、鉛直次元の軸
     **cyclic** : bool, optional
      東西境界を周期境界で扱うか。デフォルトはTrue。
     **limit** : float, optional
      デフォルトは100.
     **punit** : float, optional
      等圧面の気圧levをPaに変換するファクター。デフォルトは100.、すなわちhPaからPaへ変換。
     
    :Returns:
     **tnx, tny, tnz** : MaskedArray
      フラックスの東西、南北、鉛直成分

    .. note::
       p座標系のwave-activityフラックスは次で定義される。
      
       .. math:: \mathbf{W}=\frac{1}{2|\mathbf{U}|}
        \left( \begin{array}{c}
         U({\Psi'}_{x}^{2}-\Psi'{\Psi'}_{xx})        + V({\Psi'}_{x}{\Psi'}_{y}-\Psi'{\Psi'}_{xy})\\
         U({\Psi'}_{x}{\Psi'}_{y}-\Psi'{\Psi'}_{xy}) + V({\Psi'}_{x}^{2}-\Psi'{\Psi'}_{xx}) \\
         \frac{f^2}{S^2}[U({\Psi'}_{x}{\Psi'}_{p}-\Psi'{\Psi'}_{xp}) + V({\Psi'}_{y}{\Psi'}_{p}-\Psi'{\Psi'}_{yp})]
         \end{array} \right)
       ここでSは静的安定度。
        
    **References**
      Takaya, K and H. Nakamura, 2001: A formulation of a phase-independent wave-activity flux for stationary
      and migratory quasigeostrophic eddies on a zonally varying basic flow, `JAS`, 58, 608-627.

      <http://journals.ametsoc.org/doi/abs/10.1175/1520-0469%282001%29058%3C0608%3AAFOAPI%3E2.0.CO%3B2>

    **Examples**
     >>>
     
    """
    U, V, T = np.asarray(U), np.asarray(V), np.asarray(T)
    ndim=U.ndim
    S = stability(T, lev, zdim, punit=punit)
    f = tools.expand(constants.const_f(lat), ndim, axis=ydim)

    dstrmdx    = dvardx(strm, lon, lat, xdim, ydim, cyclic=cyclic)
    dstrmdy    = dvardy(strm, lat, ydim)
    dstrmdp    = dvardp(strm, lev, zdim, punit=punit)
    d2strmdx2  = d2vardx2(strm, lon, lat, xdim, ydim, cyclic=cyclic)
    d2strmdy2  = d2vardy2(strm, lat, ydim)
    d2strmdxdy = dvardy(dstrmdx, lat, ydim)
    d2strmdxdp = dvardx(dstrmdp, lon, lat, xdim, ydim, cyclic=True)
    d2strmdydp = dvardy(dstrmdp, lat, ydim)             

    tnx = U * (dstrmdx**2 - strm*d2strmdx2) + V * (dstrmdx*dstrmdy - strm*d2strmdxdy)
    tny = U * (dstrmdx*dstrmdy - strm*d2strmdxdy) + V * (dstrmdy**2 - strm*d2strmdy2)
    tnz = f**2/S**2 * ( U*(dstrmdx*dstrmdp - strm*d2strmdxdp) - V*(dstrmdy*dstrmdp - strm*d2strmdydp) )
    
    tnx = 0.5*tnx/np.abs(U + 1j*V)
    tny = 0.5*tny/np.abs(U + 1j*V)

    # 水平成分の大きさがlimit以上のグリッドと東風領域をマスクする。
    tnxy = np.sqrt(tnx**2 + tny**2)
    tnx = np.ma.asarray(tnx)
    tny = np.ma.asarray(tny)
    tnz = np.ma.asarray(tnz)
    tnx[(U<0) | (tnxy>limit)] = np.ma.masked
    tny[(U<0) | (tnxy>limit)] = np.ma.masked
    tnz[(U<0) | (tnxy>limit)] = np.ma.masked
    
    return tnx, tny, tnz

def qgpv(strm, t, lon, lat, lev, xdim, ydim, zdim, cyclic=True, punit=100.):
    ur"""
    準地衡ポテンシャル渦度を計算する。
    """
    lapstrm = laplacian(strm, lon, lat, xdim, ydim, cyclic=True)
     
