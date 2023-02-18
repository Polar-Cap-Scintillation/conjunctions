# examples.py
# A few short examples for using TLE and conjunction code

import datetime as dt
import numpy as np
import pymap3d as pm
from propagate_tle import TLEHandler
from satellite_conjunction import SatConj


def satellite_position():
    print('SATELLITE POSITION')

    sat_id = 39452		# NORAD satellite ID for Swarm-A

    # Create array of times to find position at
    # This example does 1 day at 1 minute intervals
    time_array = np.array([dt.datetime(2023,2,6)+dt.timedelta(minutes=m) for m in range(24*60)])

    tle = TLEHandler(sat_id)
    x, y, z = tle.sat_position(time_array)

    # convert ECEF to geodetic
    glat, glon, galt = pm.ecef2geodetic(x, y, z)

    print('GLAT:', glat)
    print('GLON:', glon)
    print('GALT:', galt)



def pfrr_conjunctions():
    print('CONJUNCTIONS')

    pfrr_site =  [65.12, -147.47, 213.]

    starttime = dt.datetime(2023,2,6)
    endtime = dt.datetime(2023,2,7)
    # dates = [startdate+dt.timedelta(days=d) for d in range((enddate-startdate).days)]

    sidA = 39452		# SID for Swarm-A
    # sidB = 39451		# SID for Swarm-B
    # sidC = 39453		# SID for Swarm-C

    SwarmA = SatConj(pfrr_site[0], pfrr_site[1], pfrr_site[2], sidA, tolerance=60., deltime=60.0)
    # SwarmB = SatConj(pfrr_site[0], pfrr_site[1], pfrr_site[2], sidB, tolerance=60., deltime=60.0)
    # SwarmC = SatConj(pfrr_site[0], pfrr_site[1], pfrr_site[2], sidC, tolerance=60., deltime=60.0)

    passesA = SwarmA.conjunctions(starttime, endtime)
    # passesB = SwarmB.conjunctions(starttime, endtime)
    # passesC = SwarmC.conjunctions(starttime, endtime)

    print(passesA)


if __name__=='__main__':

    satellite_position()
    pfrr_conjunctions()
