# propogate_tle.py

# Propigates satellite position at a particular time given a Two Line Element.
# Based primarially off of Vallado 2006
# Requires the sgp4 package (https://pypi.org/project/sgp4/)
# References:
#   Vallado, D. A., Crawford, P., Hujsak, R., and Kelso, T. S. (2006). "Revisiting Spacetrack Report #3",
#       presented at the AIAA/AAS Astrodynamics Specialist Converence, Keystone, CO, 2006 August 21-24.
#   Zhu, J. (1994). Conversion of Earth-centered Earth-fixed coordinates to geodetic coordinates.
#       IEEE Trans Aerosp Electron Syst, 30(3): 957-961. doi: 10.1109/7.303772

import numpy as np
import datetime as dt
import pymap3d as pm
from spacetrack import SpaceTrackClient
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
from sgp4.ext import jday

# import SpaceTrack username and password
from space_track_credentials import *

class TLEHandler(object):
    def __init__(self, sid):

        self.create_tle_library(sid)

    def create_tle_library(self, sid):
        # retrieve TLEs for entire period from space-track.org
        # The space-track.org API limits the number of calls you can make, so it's better to retreive the
        #   entire library for one satellite and save it as a class attribute than try to collect indivitual
        #   TLEs as needed.
        st = SpaceTrackClient(identity=ST_USERNAME, password=ST_PASSWORD)
        output = st.tle(norad_cat_id=sid, orderby='epoch asc', format='tle')

        # parse output into list of distinct TLEs
        split = output.splitlines()
        self.TLE_list = [[split[2*i],split[2*i+1]] for i in range(len(split)//2)]

        # extract epoch from each TLE
        self.TLE_epoch = [dt.datetime.strptime(tle[0][18:23],'%y%j')+dt.timedelta(days=float(tle[0][23:32])) for tle in self.TLE_list]

    def sat_position(self, time_array):

        # Find index of epoch closest to each time in the time array
        unix_time_array = np.array(time_array).astype('datetime64')
        unix_TLE_epoch = np.array(self.TLE_epoch).astype('datetime64')
        closest_epoch_idx = np.array([np.argmin(np.abs(ut-unix_TLE_epoch)) for ut in unix_time_array])

        sat_position = np.empty((3,0))

        for i in np.unique(closest_epoch_idx):
            subset_times = time_array[closest_epoch_idx==i]

            # calcualte satellite position using functions from TLE propgation script
            X, Y, Z = propagate_tle(subset_times,self.TLE_list[i])
            sat_position = np.append(sat_position, np.array([X, Y, Z]), axis=1)

        return sat_position


def propagate_tle(time0, TLE):
    # time0 is an array of datetime objects that the satellite position is to be calculated at
    # TLE is a list consisting of the first and second lines of the TLE as strings ([TLE line 1, TLE line 2])

    # Note: This function returns satellite position in Pseudo Earth-Fixed (PEF) coordinates, which are
    #   assumed to be approximately equal to Earth-Centered, Earth-Fixed (ECEF) coordinates.  This does NOT
    #   account for polar motion (precession, nutation).  For discussion of a "proper" PEF->ECEF transformation,
    #   please refer to Vallado et al., 2006 Appendix C or Panigrahi and Gaurav, 2015
    #   (https://mycoordinates.org/tracking-satellite-footprints-on-earth%E2%80%99s-surface/)

    X = []
    Y = []
    Z = []

    # initialize tle object
    tle = twoline2rv(TLE[0],TLE[1],wgs72)

    for t in time0:

        # calculate satellite position/velocity in True Equator, Mean Equinox [TEME] (units of km and km/s)
        position, velocity = tle.propagate(t.year,month=t.month,day=t.day,hour=t.hour,minute=t.minute,second=t.second)
        position_TEME = np.array(position)


        # convert to Pseudo Earth Fixed [PEF]

        # compute Julian centeries of UT1 - discussed in Vallado et al., 2006, sec. II.E
        JD = jday(t.year,t.month,t.day,t.hour,t.minute,t.second)
        T_UT1 = (JD - 2451545.0)/36525.

        # compute Greenwich Mean Sidereal Time (units of s) - Vallado et al., 2006, eqn. 2
        GMST = (67310.54841+(876600*60*60+8640184.812866)*T_UT1+0.093104*T_UT1**2-6.2e-6*T_UT1**3)
        # convert GMST to angle (units of rad)
        GMST = GMST*2*np.pi/86400. % (2*np.pi)
        # form rotational matrix
        Rot = np.array([[np.cos(GMST),np.sin(GMST),0.],[-np.sin(GMST),np.cos(GMST),0.],[0.,0.,1.]])
        # apply rotational matrix to TEME position to get PEF position (units of km) - Valladeo et al., 2006, eqn. 1
        position_PEF = np.dot(Rot,position_TEME)

        # add position to coordinate arrays
        X.append(position_PEF[0])
        Y.append(position_PEF[1])
        Z.append(position_PEF[2])

    return np.array(X)*1000., np.array(Y)*1000., np.array(Z)*1000.



def main():
    TLE = ['1     1U          18350.30892361  .00001123  00000-0  66525-4 0   109','2     1  85.0373 178.2871 0002550 225.5672 175.5175 15.21584957    13']
    times = np.array([dt.datetime(2018,12,17,0,0,0)+dt.timedelta(hours=h) for h in range(24)])

    X, Y, Z = propagate_tle(times,TLE)
    print(X, Y, Z)
    gdlat, gdlon, gdalt = pm.ecef2geodetic(X,Y,Z)

    print('{:^20}{:^10}{:^10}{:^10}'.format('Time','GLAT','GLON','GALT'))
    for t, lat, lon, alt in zip(times,gdlat,gdlon,gdalt):
        print('{}{:10.2f}{:10.2f}{:10.2f}'.format(t, lat, lon, alt))


if __name__ == '__main__':
    main()
