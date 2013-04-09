# coding:utf-8
from matplotlib.ticker import MultipleLocator
from pymet.metplt import MyBasemap
# requires netcdf4-python (netcdf4-python.googlecode.com)
from netCDF4 import Dataset

# read 500hPa Geopotential
nc = Dataset('./examples/z500_data.nc')
z500 = nc.variables['z'][0,:,:]/9.8
lon = nc.variables['longitude'][:]
lat = nc.variables['latitude'][:]

plt.subplot(2,1,1)
m = MyBasemap(lon=(0,360.01), lat=(0,90), xlint=60, ylint=30)
m.contour(lon, lat, z500, locator=MultipleLocator(100))
m.fillcontinents()

plt.title('zorder=0 (default)')

plt.subplot(2,1,2)
m = MyBasemap(lon=(0,360.01), lat=(0,90), xlint=60, ylint=30)
m.contour(lon, lat, z500, locator=MultipleLocator(100), zorder=0)
m.fillcontinents(zorder=1)
plt.title('zorder=1')

plt.show()
