# coding:utf-8
import matplotlib.pyplot as plt
from pymet.metplt import MyBasemap
from pymet.io import GradsIO
from datetime import datetime

# read data from OpenDAP
io = GradsIO(Echo=False)
io.open('http://db-dods.rish.kyoto-u.ac.jp/cgi-bin/nph-dods/ncep/ncep.reanalysis.dailyavgs/pressure/hgt.2010.nc')
io.setdim(lon=(0,357.5),lat=(-90,90),lev=500,time=datetime(2010,8,1))

# read 500hPa Geopotential
z500 = io.get('hgt')
lon, lat = z500.grid.latlon()

plt.subplot(1,1,1)
m = MyBasemap(lon=(0,360.01), lat=(0,90), xlint=60, ylint=30)#, fix_aspect=False)
m.contour(lon, lat, z500, cint=100)
m.fillcontinents()
plt.title('zorder=0 (default)')

plt.subplot(2,1,2)
m = MyBasemap(lon=(0,360.01), lat=(0,40), xlint=60, ylint=30)
m.contour(lon, lat, z500, cint=100, zorder=0)
m.fillcontinents(zorder=1)
plt.title('zorder=1')

plt.show()
