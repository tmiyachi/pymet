# coding:utf-8
import matplotlib.pyplot as plt
import pymet.metplt as metplt
from pymet.io import GradsIO
from datetime import datetime
from pymet.tools import lon2txt

# read data from OpenDAP
io = GradsIO(Echo=False)
io.open('http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/pressure/uwnd.mon.ltm.nc')
io.setdim(lon=60,lat=(-90,90),lev=(1000,100), time=datetime(1,8,1))

# read Temperature
uwnd = io.get('uwnd')
lat = uwnd.grid.lat
lev = uwnd.grid.lev

# close OpenDAP
io.close()

# plot
CF = metplt.crosscontourf(lat, lev, uwnd, xylab='lat')
plt.colorbar(CF)
plt.xticks([-90,-60,-30,0,30,60,90])
plt.title('NCEP-NCAR Reanalysis Zonal Wind Aug 1981-2010\n %s cross-section'%lon2txt(60))

plt.show()
