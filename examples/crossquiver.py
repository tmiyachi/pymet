# coding:utf-8
import matplotlib.pyplot as plt
import pymet.metplt as metplt
from pymet.io import GradsIO
from datetime import datetime
from pymet.tools import lon2txt

# read data from OpenDAP
io = GradsIO(Echo=False)
io.open('http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/pressure/vwnd.mon.ltm.nc')
io.open('http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/pressure/omega.mon.ltm.nc')
io.setdim(lon=180,lat=(-60,60),lev=(1000,100), time=datetime(1,8,1))

# read Temperature
vwnd = io.get('vwnd')
omega = io.get('omega')
#vwnd = io.get('ave(vwnd,lon=0,lon=360)')
#omega = io.get('ave(omega,lon=0,lon=360)')
lat = vwnd.grid.lat
lev = vwnd.grid.lev

# close OpenDAP
io.close()

# plot
QV = metplt.crossquiver(lat, lev, vwnd, omega, xylab='lat', xyskip=2)
metplt.crossquiverkey(QV, uvunit='m/s', wunit='Pa/s')
#plt.xticks([-60,-60,-30,0,30,60,90])
plt.title('NCEP-NCAR Reanalysis V,omega\n %s cross-section'%lon2txt(180))

plt.show()

