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
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
from sgp4.ext import jday


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
