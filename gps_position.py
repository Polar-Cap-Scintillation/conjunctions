# gps_position.py
# Functions to find the positions of GPS satellites and IPPs for a given satellite (PRN), time, and receiver.

import numpy as np
import datetime as dt
from propagate_tle import propagate_tle
from map_prn import prn2norad
import pymap3d as pm
from spacetrack import SpaceTrackClient
import spacetrack.operators as op

from space_track_credentials import *

def gps_sat_position(prn, time, coords='ECEF'):
    # calculate GPS satellite position based on PRN

    # find NORAD ID for this PRN at time 0
    norad_id = prn2norad(prn, time[0].date())

    startdate = time[0].date()
    enddate = time[-1].date()
    unix_time = np.array([(t-dt.datetime.utcfromtimestamp(0)).total_seconds() for t in time])

    # if start and end times on the same date, advance enddate by one day
    if startdate==enddate:
        enddate = enddate + dt.timedelta(days=1)

    # retrieve TLEs for entire period from space-track.org
    st = SpaceTrackClient(identity=ST_USERNAME, password=ST_PASSWORD)
    output = st.tle(norad_cat_id=norad_id, orderby='epoch asc', epoch=op.inclusive_range(startdate,enddate), format='tle')

    # parse output into list of distinct TLEs
    split = output.splitlines()
    TLE_list = [[split[2*i],split[2*i+1]] for i in range(len(split)//2)]

    # extract epoch from each TLE
    TLE_epoch = [dt.datetime.strptime(tle[0][18:23],'%y%j')+dt.timedelta(days=float(tle[0][23:32])) for tle in TLE_list]
    unix_TLE_epoch = np.array([(t-dt.datetime.utcfromtimestamp(0)).total_seconds() for t in TLE_epoch])

    # Find index of epoch closest to each time in the time array
    # this works, but causes HUGE memory problems when arrays large (when trying to find conjunctions for years at a time)
    # closest_epoch = np.argmin(np.abs(unix_time[:,None] - unix_TLE_epoch[None,:]), axis=-1)
    closest_epoch = np.array([np.argmin(np.abs(ut-unix_TLE_epoch)) for ut in unix_time])


    x = np.empty((0,))
    y = np.empty((0,))
    z = np.empty((0,))
    for i in range(len(TLE_epoch)):
        subset_times = time[closest_epoch==i]

        # calcualte satellite position using functions from TLE propgation script
        x0, y0, z0 = propagate_tle(subset_times,TLE_list[i])
        x = np.concatenate((x,x0))
        y = np.concatenate((y,y0))
        z = np.concatenate((z,z0))

    if coords == 'ECEF':
        return x, y, z
    elif coords == 'GEO':
        # convert to geodetic coordinates
        gdlat, gdlon, gdalt = pm.ecef2geodetic(x, y, z)
        return gdlat, gdlon, gdalt


def gps_azel(prn, time, site):
    # Calculate azimuth and elevation of GPS satellite based on PRN and site coordinates

    x, y, z = gps_sat_position(prn, time)
    az, el, _ = pm.ecef2aer(x, y, z, site[0], site[1], site[2])
    return az, el


def gps_ipp(prn, time, site, height=300.):
    # Calcualte IPP of GPS satellite based on PRN and site coordinates

    xsat, ysat, zsat = gps_sat_position(prn, time)
    xrec, yrec, zrec = pm.geodetic2ecef(site[0], site[1], site[2])

    vx = xsat - xrec
    vy = ysat - yrec
    vz = zsat - zrec

    earth = pm.Ellipsoid()
    a2 = (earth.semimajor_axis + height*1000.)**2
    b2 = (earth.semimajor_axis + height*1000.)**2
    c2 = (earth.semiminor_axis + height*1000.)**2

    A = vx**2/a2 + vy**2/b2 + vz**2/c2
    B = xrec*vx/a2 + yrec*vy/b2 + zrec*vz/c2
    C = xrec**2/a2 + yrec**2/b2 + zrec**2/c2 -1

    alpha = (np.sqrt(B**2-A*C)-B)/A

    lat, lon, alt = pm.ecef2geodetic(xrec + alpha*vx, yrec + alpha*vy, zrec + alpha*vz)

    return lat, lon, alt/1000.



if __name__=='__main__':

    time_array = np.array([dt.datetime(2017,11,21,18,0,0)+dt.timedelta(minutes=i) for i in range(60)])
    site = [74.4,-94.9,0.]

    print(gps_sat_position(18, time_array, coords='GEO'))
    print(gps_azel(18, time_array, site))
    print(gps_ipp(18, time_array, site))
