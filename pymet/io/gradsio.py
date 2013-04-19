# coding:utf-8
"""
"""
from grads import GaCore
from pymet.field import McField, McGrid
import numpy as np
from datetime import datetime
import os, os.path
import locale

__all__ = ['GradsIO']

class GradsIO:
    u"""
    pygradsを用いてMcFieldデータを読み込むためのクラス。

    :Arguments:
      **Echo** : bool, optional
         gradsの出力を標準出力に表示するかどうか。デフォルトはFalse。

    .. note::


    **Examples**
     >>> gaio = pymet.io.GradsIO()

    **Attrubutes**

    ======= ====================================
    ga
    ======= ====================================

    **Methods**
    .. currentmodule:: pymet.io.gradsio.GradsIO
    
    ..  autosummary::

        open
        command
        setdim
        get
        __init__

    """
    def __init__(self, Echo=False):
        locale.setlocale(locale.LC_ALL,'en_US')
        ga = GaCore('grads -b', Echo=Echo)
        self.ga = ga
        self.fn = 0   # 開いているファイル数
        self.vars = []
        self._first = True
    def open(self, fname, Quiet=True):
        u"""
        ctlファイル、またはnetCDFファイルを開く。

        :Arguments:
         **fname** : str
          ctlファイル、またはnetCDFファイルパス。netCDF形式は拡張子(.nc or .netcdf)で判断する。
         **Quiet** : bool, optional          
        """
        if not os.path.isfile(fname):
            raise IOError, "cannot find {0}".format(fname)
        root, ext = os.path.splitext(fname)
        if ext.upper() == '.NC' or ext.upper() == '.NETCDF':
            ftype='SDF'
        else:
            ftype='default'
        self.ga.open(fname, ftype=ftype, Quiet=Quiet)
        self.fn = self.fn + 1
        self.ga('set dfile '+str(self.fn))
        qfile = self.ga.query('file')
        self.vars.append(qfile.vars)
        if self._first:
            coords = self.ga.coords()
            if coords.lon[-1] - coords.lon[0] == 360.:            
                self.ga('set lon %f %f' % (coords.lon[0], coords.lon[-2]))
            self._first = False
    def close(self, fid=-1):
        u"""
        fidで指定したファイル番号のファイルを閉じる。

        :Arguments:
         **fid** : int, optional
          閉じるファイル番号。指定しない場合は、最後に開いたファイルを閉じる。
        """
        if fid > self.fn:
            raise ValueError, "the number of opened files is less than "+str(fid)
        elif fid == -1:
            fid = self.fn
        self.ga('close '+str(fid))
        self.vars.pop(fid-1)
        self.fn = self.fn - 1
    def allclose(self):
        u"""
        開いているすべてのファイルを閉じる。
        """
        self.ga('allclose')
        self.fn = 0
        self._first = True
        self.vars = []
    def command(self,command_string):
        u"""
        コマンドをGrADSに送る。

        :Arguments:
         **command_string** : string
          コマンド。
        """
        self.ga(command_string)
    def setdim(self, time=None, lon=None, lat=None, lev=None, ens=None):
        u"""
        次元を指定する。        
        """
#        qh = self.ga.query('ctlinfo')
        locale.setlocale(locale.LC_ALL,'en_US')
        if isinstance(lon, tuple):
            self.ga('set lon %f %f' % lon)
        elif lon==None:
            coords = self.ga.coords()
            if coords.lon[-1] - coords.lon[0] == 360.:            
                self.ga('set lon %f %f' % (coords.lon[0], coords.lon[-2]))
        else:
            self.ga('set lon %f' % lon)
        if isinstance(lat, tuple):
            self.ga('set lat %f %f' % lat)
        elif lat!=None:
            self.ga('set lat %f' % lat)
        if isinstance(lev, tuple):
            self.ga('set lev %f %f' % lev)
        elif lev!=None:
            self.ga('set lev %f' % lev)
        if isinstance(time, tuple):
            self.ga('set time %s %s' % (time[0].strftime('%Hz%d%b%Y'), time[1].strftime('%Hz%d%b%Y')))
        elif time != None:
            self.ga('set time %s' % time.strftime('%Hz%d%b%Y'))
        if isinstance(ens, tuple):
            self.ga('set e %d %d' % ens)
        elif type(ens) == type(int(1)):
            self.ga('set e %d' % ens)
        elif ens=='all':
            self.ga('set e 1 %d' % self.ga.query('file',Quiet=True).ne)
        
    def get(self, var, fid=None):
        u"""
        指定した変数をMcFieldオブジェクトとして取得する。
        """
        if fid==None:
            for i, fvars in enumerate(self.vars):
                if var in fvars:
                    fid = i+1
                    break
        self.ga('set dfile '+str(fid))
        dh = self.ga.query('dims', Quiet=True)
        info = self.ga.coords()
        ne, nt, nz, ny, nx = dh.ne, dh.nt, dh.nz, dh.ny, dh.nx
        e1, e2 = dh.ei
        t1, t2 = dh.ti
        z1, z2 = dh.zi
        
        lon  = np.asarray(info.lon, dtype=np.float32)
        lat  = np.asarray(info.lat, dtype=np.float32)
        lev  = np.asarray(info.lev, dtype=np.float32)
        time = np.asarray([datetime.strptime(time,'%HZ%d%b%Y') for time in info.time])
        ens  = np.asarray(info.ens)        

        out = np.zeros((ne,nt,nz,ny,nx))
        ## 4次元以上はgradsでは同時に扱えないのでループする
        ## 少ない次元を優先的にループ
        if nz == max(ne, nt, nz):
            try:
                for i, e in enumerate(range(e1,e2+1)):
                    self.ga.cmd('set e %d' % e)
                    for j, t in enumerate(range(t1,t2+1)):
                        self.ga.cmd('set t %d' % t)
                        out[i,j,:,:,:]  = np.asarray(self.ga.eval(var), dtype=np.float32).reshape((nz,ny,nx))
                        self.ga.Writer.flush()
            except:
                self.ga.flush()
                self.ga.setdim(dh)
                raise
        elif nt == max(ne, nt, nz):
            try:
                for i, e in enumerate(range(e1,e2+1)):
                    self.ga.cmd('set e %d' % e)
                    for k, z in enumerate(range(z1,z2+1)):
                        self.ga.cmd('set z %d' % z)
                        out[i,:,k,:,:]  = np.asarray(self.ga.eval(var), dtype=np.float32).reshape((nt,ny,nx))
                        self.ga.Writer.flush()
            except:
                self.ga.flush()
                self.ga.setdim(dh)
                raise
        else:
            try:
                for j, t in enumerate(range(t1,t2+1)):
                    self.ga.cmd('set t %d' % t)
                    for k, z in enumerate(range(z1,z2+1)):
                        self.ga.cmd('set z %d' % z)
                        out[:,j,k,:,:]  = np.asarray(self.ga.eval(var), dtype=np.float32).reshape((ne,ny,nx))
                        self.ga.Writer.flush()
            except:
                self.ga.flush()
                self.ga.setdim(dh)
                raise
        out  = np.squeeze(out)
        out  = np.ma.array(out, mask=(out==info.undef))

        self.ga.flush()
        self.ga.setdim(dh)

        # create McGrid
        grid = McGrid(name=var, lon=lon, lat=lat, lev=lev, time=time, ens=ens)
        field = McField(out, name=var, grid=grid, mask=out.mask)

        return field

        
