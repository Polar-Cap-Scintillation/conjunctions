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
from apexpy import Apex
from spacetrack import SpaceTrackClient
# import spacetrack.operators as op
from propagate_tle import propagate_tle
from space_track_credentials import *


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
        self.create_tle_library(sat_id)

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


        # Find index of epoch closest to each time in the time array
        unix_TLE_epoch = np.array([(t-dt.datetime.utcfromtimestamp(0)).total_seconds() for t in self.TLE_epoch])
        # this works, but causes HUGE memory problems when arrays large (when trying to find conjunctions for years at a time)
        # closest_epoch = np.argmin(np.abs(unix_time_array[:,None] - unix_TLE_epoch[None,:]), axis=-1)
        closest_epoch_idx = np.array([np.argmin(np.abs(ut-unix_TLE_epoch)) for ut in unix_time_array])


        sat_position = np.empty((0,3))

        # for i in range(len(self.TLE_epoch)):
        for i in np.unique(closest_epoch_idx):
            subset_times = time_array[closest_epoch_idx==i]

            # calcualte satellite position using functions from TLE propgation script
            X, Y, Z = propagate_tle(subset_times,self.TLE_list[i])
            sat_position = np.append(sat_position, np.array([X, Y, Z]).T, axis=0)

        # define site
        if self.conjcoords == 'geo':
            site_lat = self.site_lat
            site_lon = self.site_lon
            zen_vec = self.site_zenith
        if self.conjcoords == 'mag':
            A = Apex(starttime)
            site_lat, site_lon = A.geo2apex(self.site_lat, self.site_lon, self.site_alt/1000.)
            _,_,_,_,_,_, d1, d2, d3, e1, e2, e3 = A.basevectors_apex(self, self.site_lat, self.site_lon, self.site_alt/1000.)
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
            sat_lat, sat_lon, sat_alt = pm.ecef2geodetic(sat_position[0], sat_position[1], sat_position[2])
            if self.conjcoords == 'mag':
                sat_lat, sat_lon = A.geo2apex(sat_lat, sat_lon, sat_alt/1000.)

            # find all times satellite position is within tolerance of the radar position
            conjunctions = np.all([np.abs(sat_lat-radar_lat)<latlontol[0],np.abs(sat_lon-radar_lon)<latlontol[1]], axis=0)


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
