# coding:utf-8
"""
grib1,grib2形式のデータをMcFieldオブジェクトとして読み込むためのモジュール。
gribの読み込みにはpygribを使用。
"""
import pygrib
import numpy as np
from pymet.field import McGrid, McField

class gribfile:
    def __init__(self):
        self.gribslist = []
    def open(self, fname):
        grbs = pygrib.open(fname)




def getgrib():
    pass



def gettigge(center,idate,var,lev,ftyp='cntl'):
    YYYYMM=idate.strftime('%Y%m')
    YYYYMMDDHH=idate.strftime('%Y%m%d%H')
    varname = var+str(lev)
    fname=FNAMEBASE.format(center=center.lower(),
                           ftyp=ftyp,
                           varname=varname,
                           YYYYMM=YYYYMM,
                           YYYYMMDDHH=YYYYMMDDHH)
    print fname
    #
    shortname=var
    if var=='z': shortname='gh'
    grbs = pygrib.open(fname)
    grb = grbs.select(shortName=shortname,level=lev)
    xn = grb[0]['Ni']
    yn = grb[0]['Nj']
    tn = len(grb)
    lat,lon = grb[0].latlons()
    out = np.empty((tn,yn,xn))
    tyme = np.empty(tn,dtype=np.object)
    for t, gr in enumerate(grb):
        out[t,...] = gr.values
        tyme[t] = gr.validDate
    grbs.close()
    dims = ['tyme','lat','lon']
    grid = McGrid(out, lon=lon[0,:], lat=lat[::-1,0], lev=lev, tyme=tyme, dims=dims)
    out = McField(out,name=var,grid=grid)
    return out
