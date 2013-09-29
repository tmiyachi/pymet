# coding:utf-8
"""
"""
from grads import GaCore, GrADSError
from pymet.field import McField, McGrid
from pymet.tools import d2s, s2d
import numpy as np
from datetime import datetime
import os, os.path

__all__ = ['GradsIO']

class GradsIO:
    u"""
    pygradsを用いてMcFieldデータを読み込むためのクラス。

    :Arguments:
      **Echo** : bool, optional
         gradsの出力を標準出力に表示するかどうか。デフォルトはFalse。

    **Examples**
     >>> gaio = pymet.io.GradsIO(Echo=False)
     >>> gaio.open('u.ctl')
     >>> gaio.open('v.ctl')
     >>> ga.setdim(lon=(0,180), lat=(-90,90), lev=500, time=datetime(2009,10,13,12))
     >>> u = ga.get('u')
     >>> v = ga.get('v')
     >>> gaio.allclose()

    **Attributes**
    
     ======= ===================================================================
     ga      grads.GaCoreクラスのインスタンス
     fn      現在開いているファイル数
     fnames  現在開いているファイル名
     vars    現在開いているファイルに含まれる変数名
     gdim
     ======= ===================================================================

    **Methods**
     ..  autosummary::

         open
         close
         allclose
         command
         setdim
         setgradsdim
         get
    """
    def __init__(self, Echo=False):
        ga = GaCore(Window=False, Echo=Echo)
        self.ga = ga
        self.fn = 0         # 開いているファイル数
        self.fn_total = 0   # インスタンスが開いたファイル数の合計
        self.vars = {}
        self.fnames = {}
        self._first = True
        
    def open(self, fname, Quiet=True):
        u"""
        ctlファイル、またはnetCDFファイルを開く。

        :Arguments:
         **fname** : str
          ctlファイル、またはnetCDFファイルパス。netCDF形式は拡張子(.nc or .netcdf)で判断する。
         **Quiet** : bool, optional          
        """
        ga = self.ga
        if fname[:7]!='http://' and not os.path.isfile(fname):
            raise IOError, "cannot find {0}".format(fname)
    
        root, ext = os.path.splitext(fname)
        if ext.upper() == '.NC' or ext.upper() == '.NETCDF':
            ftype='SDF'
        else:
            ftype='default'
        ga.open(fname, ftype=ftype, Quiet=Quiet)
        
        self.fn = self.fn + 1
        self.fn_total = self.fn_total + 1        
        ga.cmd('set dfile %d' % self.fn)
        qfile = ga.query('file')
        fid = self.fn_total
        self.vars[fid] = qfile.vars
        self.fnames[fid] =fname

        if self._first:
            coords = ga.coords()
            denv = coords.denv
            lon, lat, lev, time, ens = denv.lon, denv.lat, denv.lev, denv.tyme, denv.e
            if coords.lon[-1] - coords.lon[0] == 360.:
                lon = (coords.lon[0], coords.lon[-2])
            self.gdim = GradsDim(lon=lon,lat=lat,lev=lev,time=time,ens=ens)
            self._first = False
    
        return fid
            
    def close(self):
        u"""
        一番最後に開いたファイルを閉じる。
        """
        self.ga.cmd('close %d' % self.fn)
        lastfid = max(self.fnames.keys())
        self.vars.pop(lastfid)
        self.fnames.pop(lastfid)
        self.fn = self.fn - 1
        
    def allclose(self):
        u"""
        開いているすべてのファイルを閉じる。
        """
        for fid in range(self.fn,0,-1):
            self.ga.cmd('close %d' % fid)
        self.fn = 0
        self._first = True
        self.vars = {}
        self.fnames = {}
        
    def command(self,command_string):
        u"""
        コマンドをGrADSに送る。

        :Arguments:
         **command_string** : string
          コマンド。
        """
        self.ga.cmd(command_string)

    def setdim(self, **kwargs):
        u"""
        次元を設定する。

        :Arguments:
         **lon** : float or tuple
          経度次元
         **lat** : float or tuple
          緯度次元
         **lev** : float or tuple
          鉛直次元
         **time** : datetime object or tuple
          時間次元
         **ens** : int or tuple or 'all'
          アンサンブル次元。'all'を指定するとファイルのもつ全てのアンサンブル次元が指定される。
        .. note::
         次元情報は、GradsIOオブジェクトのgdim属性にGradsDimクラスのオブジェクトとして設定される。
         
        """
        gdim = self.gdim
        for name, value in kwargs.items():
            setattr(gdim, name, value)
        
    def setgradsdim(self):
        u"""
        gdim属性がもつGradsDimオブジェクトの情報を次元をGrADSに設定する。
        """
        ga = self.ga
        gdim = self.gdim
        if gdim.lon: ga.cmd('set lon %f %f' % gdim.lon)
        if gdim.lat: ga.cmd('set lat %f %f' % gdim.lat)
        if gdim.lev: ga.cmd('set lev %f %f' % gdim.lev)
        if gdim.time: ga.cmd('set time %s %s' % (d2s(gdim.time[0],fmt='%H:%MZ%d%b%Y'),d2s(gdim.time[1],fmt='%H:%MZ%d%b%Y')))
        if gdim.ens:
            if gdim.ens == 'all':
                ga.cmd('set e 1 %d' % ga.query('file',Quiet=True).ne)
            else:
                ga.cmd('set e %d %d' % gdim.ens)

        # test
        dh = ga.query('dims', Quiet=True)
        fdh = ga.query('ctlinfo', Quiet=True)

        x1, x2 = dh.xi
        if x1 < 1 or x2 > fdh.nx:                
            ga.cmd('set x %d %d' % (max(x1,1), min(x2, fdh.nx)))
        y1, y2 = dh.yi
        if y1 < 1 or y2 > fdh.ny:
            ga.cmd('set y %d %d' % (max(y1,1), min(y2, fdh.ny)))

                
    def get(self, var, fid=None, **kwargs):
        u"""
        指定した変数をMcFieldオブジェクトとして取得する。

        :Arguments:
         **var** : str
          変数名または、ave(変数名, dim=, dim=)。
         **fid** : int, optional
          ファイルを複数開いている場合にファイルidを指定する。指定しない場合は新しく開いたファイルから順に探して、
          指定した変数を含む最もファイルidの小さいファイルから読み込む。
        :Returns:
         **field** : McField object
        """
        ga = self.ga
        # 開かれているファイルの変数リストからvarを探す
        if 'ave(' in var:
            gavar = var[4:].split(',')[0]
        else:
            gavar = var

        if fid:
            if self.vars[fid].count(gavar) < 1:
                raise ValueError, "Cannot find variable {0} in opened files".format(gavar)
        else:
            # 現在開いている各ファイルの変数リストの辞書から、gavarを含むファイルのfidを取り出す
            var_in_fids = [opend_fid for opend_fid, opend_vars in self.vars.items() if gavar in opend_vars]
            if len(var_in_fids) == 0:
                raise ValueError, "Cannot find variable '{0}' in opened files".format(gavar)
            elif len(var_in_fids) > 1: 
                print "Warnig, multiple variables '{0}' are fond in opened files".format(gavar)
            fid = var_in_fids[-1]

        # デフォルトファイルを変更
        num_closed_file = self.fn_total - self.fn
        ga.cmd('set dfile %d' % (fid - num_closed_file))

        # 次元をGrADSに送って設定する
        self.setdim(**kwargs)
        self.setgradsdim()

        # 有効な次元かどうか確かめる(要検討)
        self._checkdim()

        # 次元を調べる
        dh = ga.query('dims', Quiet=True)
        info = ga.coords()

        ne, nt, nz, ny, nx = dh.ne, dh.nt, dh.nz, dh.ny, dh.nx
        e1, e2 = dh.ei
        t1, t2 = dh.ti
        z1, z2 = dh.zi

        # 次元の値を取得する
        lon  = np.asarray(info.lon, dtype=np.float32)
        lat  = np.asarray(info.lat, dtype=np.float32)
        lev  = np.asarray(info.lev, dtype=np.float32)
        time = np.asarray([s2d(d) for d in info.time])
        ens  = np.arange(e1, e2+1)

        # データを取得
        out = np.ma.zeros((ne,nt,nz,ny,nx))
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
                raise GrADSError, "Syntax Error"
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
                raise GrADSError, "Syntax Error"
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
                raise GrADSError, "Syntax Error"
        out  = np.squeeze(out)
        out  = np.ma.array(out, mask=(out==info.undef))
        #       self.ga.flush()
        self.ga.setdim(dh)

        # McFieldオブジェクトを作成
        grid = McGrid(name=var, lon=lon, lat=lat, lev=lev, time=time, ens=ens)
        field = McField(out, name=var, grid=grid, mask=out.mask)
        
        return field

    def _checkdim(self):
        u"""
        現在の次元の設定が正しくなければ例外を送出
        """
        ga = self.ga
        
        dh = ga.query('dims', Quiet=True)

        fdh = ga.query('ctlinfo', Quiet=True)

        x1, x2 = dh.xi
        if x1 > x2:
            raise GrADSError, "requested x-dimension is invalid. "
        if x1 < 1 or x2 > fdh.nx:                
            raise GrADSError, "requested x-dimension is out of bounds"
        y1, y2 = dh.yi
        if y1 < 1 or y2 > fdh.ny:
            raise GrADSError, "requested y-dimension is out of bounds"
        z1, z2 = dh.zi
        if (fdh.nz != 1) and (z1 < 1 or z2 > fdh.nz):
            raise GrADSError, "requested z-dimension is out of bounds"
        t1, t2 = dh.ti
        if t1 < 1 or t2 > fdh.nt:
            raise GrADSError, "requested t-dimension is out of bounds"
        e1, e2 = dh.ei
        if e1 < 1 or e2 > fdh.ne:
            raise GrADSError, "requested e-dimension is out of bounds"

class GradsDim(object):
    def __init__(self, **kwargs):
        self.lon  = _num2tuple(kwargs.get('lon', None))
        self.lat  = _num2tuple(kwargs.get('lat', None))
        self.lev  = _num2tuple(kwargs.get('lev', None))
        self.time = _num2tuple(kwargs.get('time', None))
        self.ens  = _num2tuple(kwargs.get('ens', None))                                

    def __setattr__(self, name, value):
        if name in ['lon', 'lat', 'lev', 'time', 'ens']:
            self.__dict__[name] = _num2tuple(value)
        else:
            self.__dict__[name] = value

    def __repr__(self):
        repr = ''
        repr += 'lon  : ' + str(self.lon) + '\n'
        repr += 'lat  : ' + str(self.lat) + '\n'
        repr += 'lev  : ' + str(self.lev) + '\n'
        repr += 'time : ' + str(self.time) + '\n'
        repr += 'ens  : ' + str(self.ens)

        return repr

def _num2tuple(a):
    u"""
    GradsDim用内部メソッド。
    """
    if isinstance(a, tuple) or isinstance(a, list):
        if len(a) != 2:
            raise ValueError, "dimension tuple must contain 2 value"
        return a
    elif not isinstance(a, str) and a != None:
        return (a, a)
    else:
        return a
