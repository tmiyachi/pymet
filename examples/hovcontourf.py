# coding:utf-8
import matplotlib.pyplot as plt
import pymet.metplt as metplt
from pymet.io import GradsIO
from datetime import datetime

# read data from OpenDAP
io = GradsIO(Echo=False)
#io.open('http://db-dods.rish.kyoto-u.ac.jp/cgi-bin/nph-dods/ncep/ncep.reanalysis.derived/pressure/uwnd.mon.ltm.nc')
io.open('http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/pressure/vwnd.2010.nc')
#io.setdim(lon=(0, 130),lat=(37.5,42.5),lev=200, time=datetime(1,8,1))
io.setdim(lon=(0,130),lat=40, lev=200, time=(datetime(2010,7,1),datetime(2010,8,1)))

# read data
data = io.get('vwnd')
lon = data.grid.lon
time = data.grid.time

# close OpenDAP
io.close()

# plot
CF = metplt.hovcontourf(lon, time, data, xylabel='lat', fmt='%d%b')
plt.colorbar(CF)
plt.title('NCEP-NCAR Reanalysis v200 AUG2010')

plt.show()
