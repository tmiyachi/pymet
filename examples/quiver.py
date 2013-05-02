# coding:utf-8
import matplotlib.pyplot as plt
from pymet.metplt import MyBasemap
from pymet.io import GradsIO
from datetime import datetime

# read data from OpenDAP
io = GradsIO(Echo=False)
io.open('http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/pressure/uwnd.mon.mean.nc')
io.open('http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/pressure/vwnd.mon.mean.nc')
io.setdim(lon=(0,357.5),lat=(-90,90),lev=500,time=datetime(2010,8,1))

# read 200hPa u, v
u200 = io.get('uwnd')
v200 = io.get('vwnd')
lon, lat = u200.grid.latlon()

# close OpenDAP
io.close()

# plot
m = MyBasemap(lon=(0,360), lat=(-90,90), xlint=60, ylint=30, projection='cyl')
QV = m.quiver(lon, lat, u200, v200, skip=4, latlon=True)
m.quiverkey(QV, unit='m/s')
plt.title('NCEP-NCAR Reanalysis U,V200 Aug2010')

plt.show()
