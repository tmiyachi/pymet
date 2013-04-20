# coding:utf-8
import matplotlib.pyplot as plt
import pymet.metplt as metplt
from pymet.io import GradsIO
from datetime import datetime

# read data from OpenDAP
io = GradsIO(Echo=False)
#io.open('http://db-dods.rish.kyoto-u.ac.jp/cgi-bin/nph-dods/ncep/ncep.reanalysis.derived/pressure/uwnd.mon.ltm.nc')
io.open('http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/pressure/uwnd.mon.ltm.nc')
io.setdim(lon=60,lat=(-90,90),lev=(1000,10), time=datetime(1,8,1))

# read Temperature
uwnd = io.get('uwnd')
lat = uwnd.grid.lat
lev = uwnd.grid.lev

# close OpenDAP
io.close()

# plot
CF = metplt.crosscontourf(lat, lev, uwnd, xylabel='lat')
plt.colorbar(CF)
plt.title('NCEP-NCAR Reanalysis Zonal Wind Aug 1981-2010\n60E cross-section')

plt.show()
