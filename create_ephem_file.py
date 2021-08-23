# create_ephem_file.py
# Create a simple file of positions at given times based on a TLE

import numpy as np
import datetime as dt
from propagate_tle import propagate_tle
import pymap3d as pm
from apexpy import Apex
from spacetrack import SpaceTrackClient
import spacetrack.operators as op

from space_track_credentials import *



starttime = dt.datetime(2017,11,21,18,42)
endtime = dt.datetime(2017,11,21,18,48)
deltime = 1/20.     # time steps in seconds
sid = 39265         # NORAD satellite ID
outfile = 'CERTO_ephem.txt'


startdate = starttime.date()
enddate = endtime.date()

# if start and end times on the same date, advance enddate by one day
if startdate==enddate:
    enddate = enddate + dt.timedelta(days=1)

# retrieve TLEs for entire period from space-track.org
st = SpaceTrackClient(identity=ST_USERNAME, password=ST_PASSWORD)
output = st.tle(norad_cat_id=sid, orderby='epoch asc', epoch=op.inclusive_range(startdate,enddate), format='tle')

# parse output into list of distinct TLEs
split = output.splitlines()
TLE_list = [[split[2*i],split[2*i+1]] for i in range(len(split)//2)]

# extract epoch from each TLE
TLE_epoch = [dt.datetime.strptime(tle[0][18:23],'%y%j')+dt.timedelta(days=float(tle[0][23:32])) for tle in TLE_list]
unix_TLE_epoch = np.array([(t-dt.datetime.utcfromtimestamp(0)).total_seconds() for t in TLE_epoch])

# form array of times to find the satellite location
unix_time_array = (starttime-dt.datetime.utcfromtimestamp(0)).total_seconds()+np.arange(0,(endtime-starttime).total_seconds(),deltime)
time_array = np.array([dt.datetime.utcfromtimestamp(ut) for ut in unix_time_array])

# Find index of epoch closest to each time in the time array
# this works, but causes HUGE memory problems when arrays large (when trying to find conjunctions for years at a time)
# closest_epoch = np.argmin(np.abs(unix_time_array[:,None] - unix_TLE_epoch[None,:]), axis=-1)
closest_epoch = np.array([np.argmin(np.abs(ut-unix_TLE_epoch)) for ut in unix_time_array])


x = np.empty((0,))
y = np.empty((0,))
z = np.empty((0,))
for i in range(len(TLE_epoch)):
    subset_times = time_array[closest_epoch==i]

    # calcualte satellite position using functions from TLE propgation script
    x0, y0, z0 = propagate_tle(subset_times,TLE_list[i])
    x = np.concatenate((x,x0))
    y = np.concatenate((y,y0))
    z = np.concatenate((z,z0))


# convert to geodetic coordinates
gdlat, gdlon, gdalt = pm.ecef2geodetic(x, y, z)

# convert to Apex coordinates
A = Apex(starttime)
alat, alon = A.geo2apex(gdlat, gdlon, height=gdalt/1000.)

# save data to txt file
outdata = np.array([unix_time_array, x, y, z, gdlat, gdlon, gdalt, alat, alon]).T
header = ''.join([f'{t:^16}' for t in ['Unix Time','X','Y','Z','GDLAT','GDLON','GDALT','ALAT','ALON']])
np.savetxt(outfile, outdata, fmt='%15.5f',header=header)

# # For double checking output:
# import matplotlib.pyplot as plt
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
#
# proj = ccrs.Orthographic(central_longitude=-94.9, central_latitude=74.7)
# fig = plt.figure()
# ax_map = fig.add_subplot(111,projection=proj)
# gl = ax_map.gridlines(color='dimgrey')
# ax_map.coastlines(linewidth=0.5, color='dimgrey')
# ax_map.add_feature(cfeature.LAND, color='lightgrey')
# ax_map.add_feature(cfeature.OCEAN, color='white')
# ax_map.set_extent([-110.,-80.,70.,80.])
# ax_map.plot(gdlon, gdlat, transform=ccrs.Geodetic())
# plt.show()
