# coding:utf-8
#import matplotlib.pyplot as plt
from pymet.metplt import MyBasemap
# requires netcdf4-python (netcdf4-python.googlecode.com)
from netCDF4 import Dataset

# read 500hPa Geopotential
nc = Dataset('../../examples/uv500_data.nc')
u500 = nc.variables['u'][0,:,:]
v500 = nc.variables['v'][0,:,:]
lon = nc.variables['longitude'][:]
lat = nc.variables['latitude'][:]

m = MyBasemap(lon=(0,360), lat=(-90,90), xlint=60, ylint=30)
QV=m.quiver(lon, lat, u500, v500, skip=4)
m.quiverkey(QV, unit='m/s')
plt.title('500hPa wind')

plt.show()
