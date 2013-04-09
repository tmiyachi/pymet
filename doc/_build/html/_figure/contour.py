# coding:utf-8
#import matplotlib.pyplot as plt
from pymet.metplt import MyBasemap
# requires netcdf4-python (netcdf4-python.googlecode.com)
from netCDF4 import Dataset

# read 500hPa Geopotential
nc = Dataset('../../examples/z500_data.nc')
z500 = nc.variables['z'][0,:,:]/9.8
lon = nc.variables['longitude'][:]
lat = nc.variables['latitude'][:]

m = MyBasemap(lon=(0,360), lat=(-90,90), xlint=60, ylint=30)
m.contour(lon, lat, z500, cint=100, unit='m')
plt.title('500hPa Geopotential height')

plt.show()
