# satellite_conjunctions.py
# Determine satellite conjunctions with an site within a given tolerance.
# You must have an account at space-track.org (free and easy to set up)
# Requires the following non-standard packages:
#   - apexpy (https://apexpy.readthedocs.io/en/latest/readme.html)
#   - pymap3d (https://geospace-code.github.io/pymap3d/)
#   - spacetrack (https://pythonhosted.org/spacetrack/)
#   - propagate_tle.py
# Note: This can take a very long time to calculate tight conjunctions over many years of data
# This is an excellent site for validation: https://swarm-aurora.com/satelliteFinder/


import datetime as dt
import numpy as np
import pymap3d as pm
try:
    from apexpy import Apex
except:
    print('Could not import apexpy - cannot calculate magnetic conjuntions.')
from propagate_tle import TLEHandler


class SatConj(object):

    def __init__(self, site_lat, site_lon, site_alt, sat_id, deltime=60., conjtype='zenith', tolerance=25., coords='geo'):
    # deltime = time step between calculations of the satellite position
    # conjtype = method to use for identifying conjunctions, either 'zenith' to find conjunctions within a certain angle of
    #   zenith or 'latlon' to find conjunctions within a certain lat/lon box of the site
    # tolerance = tolerance for conjunctions (angle in degrees for 'zenith' or 2-element list of dlat/dlon for 'latlon')
    # coords = whether to perform calculation in geodetic ('geo') or apex magnetic ('mag') coordinats
    #
    # The selections for deltime & tolerance can strongly impact the results.  The smaller deltime is, the finner
    #   the satellite position calculations will be, which reduces the change of "skipping over" a potential conjuction,
    #   but may increase the time it takes to run the code. Smaller values for tolerance forces much tigher conjuctions,
    #   which will likely reduce the total number of conjunctions.  Ideally, deltime should be selected so that several
    #   steps (footprint velocity x deltime) fit in the tolerance box at the radar's location.


        self.deltime = deltime
        self.conjtype = conjtype
        self.tolerance = tolerance
        self.conjcoords = coords
        self.site_coordinates(site_lat, site_lon, site_alt)
        self.tle = TLEHandler(sat_id)


    def site_coordinates(self, site_lat, site_lon, site_alt):
        # find some basic site paramters
        self.site_coords = np.array(pm.geodetic2ecef(site_lat, site_lon, site_alt))

        self.site_lat = site_lat
        self.site_lon = site_lon
        self.site_alt = site_alt
        self.site_zenith = np.array(pm.enu2uvw(0., 0., 1., site_lat, site_lon))


    def conjunctions(self, starttime, endtime):
        # form array of times to find the satellite location
        # deltime selection is important - large values may miss some conjunctions, small values will take a very long time
        # and appropriate value of deltime depends on the velocity (altitude) of the satellite
        ustarttime = (starttime-dt.datetime.utcfromtimestamp(0)).total_seconds()
        uendtime = (endtime-dt.datetime.utcfromtimestamp(0)).total_seconds()
        unix_time_array = np.arange(ustarttime, uendtime, self.deltime)
        time_array = np.array([dt.datetime.utcfromtimestamp(ut) for ut in unix_time_array])

        sat_position = self.tle.sat_position(time_array).T


        # define site
        if self.conjcoords == 'geo':
            site_lat = self.site_lat
            site_lon = self.site_lon
            zen_vec = self.site_zenith
        if self.conjcoords == 'mag':
            A = Apex(starttime)
            site_lat, site_lon = A.geo2apex(self.site_lat, self.site_lon, self.site_alt/1000.)
            _,_,_,_,_,_, d1, d2, d3, e1, e2, e3 = A.basevectors_apex(self.site_lat, self.site_lon, self.site_alt/1000.)
            zen_vec = np.array(pm.enu2uvw(e3[0], e3[1], e3[2], self.site_lat, self.site_lon))
            zen_vec = zen_vec/np.linalg.norm(zen_vec)

        if self.conjtype == 'zenith':
            # Find vector to satelite and zenith angle
            sat_vec = sat_position - self.site_coords
            zenith_angle = np.arccos(np.einsum('...i,i->...',sat_vec, zen_vec)/np.linalg.norm(sat_vec,axis=1))*180./np.pi

            # Find where zeith angle within tolerance
            conjunctions = zenith_angle < self.tolerance

        if self.conjtype == 'latlon':
            # Find satellite position in lat/lon coordinates
            sat_lat, sat_lon, sat_alt = pm.ecef2geodetic(sat_position[:,0], sat_position[:,1], sat_position[:,2])
            if self.conjcoords == 'mag':
                sat_lat, sat_lon = A.geo2apex(sat_lat, sat_lon, sat_alt/1000.)

            # find all times satellite position is within tolerance of the radar position
            conjunctions = np.all([np.abs(sat_lat-self.site_lat)<self.tolerance[0],np.abs(sat_lon-self.site_lon)<self.tolerance[1]], axis=0)


        # seperate individual passes
        diffs = unix_time_array[conjunctions][1:]-unix_time_array[conjunctions][:-1]
        jumps = np.argwhere(diffs>self.deltime).flatten()

        if not jumps.size:
            if not np.any(conjunctions):
                # if no passes in experiment
                passes = []
            else:
                # if only one pass in experiment
                passes = [{'time':time_array[conjunctions],'position':sat_position[conjunctions][:]}]
        else:
            # first pass
            passes = [{'time':time_array[conjunctions][:jumps[0]+1],'position':sat_position[conjunctions][:jumps[0]+1,:]}]
            # all intermediate passes
            for i in range(len(jumps)-1):
                passes.append({'time':time_array[conjunctions][jumps[i]+1:jumps[i+1]+1],'position':sat_position[conjunctions][jumps[i]+1:jumps[i+1]+1,:]})
            # last pass
            passes.append({'time':time_array[conjunctions][jumps[-1]+1:],'position':sat_position[conjunctions][jumps[-1]+1:,:]})

        return passes
