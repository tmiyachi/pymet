# coding:utf-8
from matplotlib.ticker import MultipleLocator
from pymet.plot import MyBasemap
# requires netcdf4-python (netcdf4-python.googlecode.com)
from netCDF4 import Dataset

# read 500hPa Geopotential
nc = Dataset('../../examples/z500_data.nc')
z500 = nc.variables['z'][0,:,:]/9.8
lon = nc.variables['longitude'][:]
lat = nc.variables['latitude'][:]

m = MyBasemap(lon=(0,360.01), lat=(-90,90), xlint=60, ylint=30)
CF = m.contourf(lon, lat, z500, locator=MultipleLocator(100))
m.colorbar(CF)
plt.title('500hPa Geopotential height')

plt.show()
