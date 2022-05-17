# amisr_sat_conj.py
# Calculate conjunctions when an AMISR was operating and a satellite passes close to the radar.
# This uses amisr_mode.py and satellite_conjunction.py.

import numpy as np
import datetime as dt
# from apexpy import Apex
import pymap3d as pm
# from spacetrack import SpaceTrackClient
# import spacetrack.operators as op
# import sqlalchemy as db
# from propagate_tle import propagate_tle
import amisr_lookup
import satellite_conjunction

import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# from space_track_credentials import *



# def conjunctions(starttime, endtime, radar, sid, deltime=60., conjtype='zenith', tolerance=25., conjcoords='geo'):
def conjunctions(starttime, endtime, radar, sid, **sat_conj_kwarg):
    # returns a list of conjuctions between a satellite (SID is the NORAD satellite ID) and ground-based AMISR (when the radar
    #    is running) as well as information about the radar experiment and mode
    # Input:
    #   starttime = begining of time period to be considered
    #   endtime = end of time period to be considered
    #   radar = AMISR to find conjuncitons with
    #   sid = NORAD satellite ID (space-track.org has a search tool)
    #   Also will take any optional keyword arguments for the SatConj class

    # Create AMISR_looup object
    al = amisr_lookup.AMISR_lookup(radar)

    # Initialize satellite conjunciton instance
    radar_site = al.site_coords()
    conj = satellite_conjunction.SatConj(radar_site[0], radar_site[1], 0., sid, **sat_conj_kwarg)

    # Get list of all AMISR experiments in specified time frame
    exp_list = al.find_experiments(starttime, endtime)

    # Find passes within each experiment
    pass_list = []
    for exp in exp_list:

        exp_start = dt.datetime.utcfromtimestamp(exp['start_time'])
        exp_end = dt.datetime.utcfromtimestamp(exp['end_time'])

        passes = conj.conjunctions(exp_start, exp_end)

        for p in passes:
            p.update({'mode':exp['mode'], 'experiment':exp['experiment_number']})

        pass_list.extend(passes)

    return pass_list


def validate(starttime, endtime, radar, sid, deltime=60., zatol=25., latlontol=None, conjcoords='geo'):

    passes = conjunctions(starttime, endtime, radar, sid, deltime=deltime, lattol=lattol, lontol=lontol, conjcoords=conjcoords, return_all=True)

    mode = am.AMISR_mode(radar)
    radar_site = mode.site_coords()


    mapproj = ccrs.LambertConformal(central_latitude=radar_site.latitude, central_longitude=radar_site.longitude)

    for p in passes:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection=mapproj)
        ax.coastlines()
        ax.gridlines()
        # ax.set_extent([radar_site.longitude-lontol, radar_site.longitude+lontol, radar_site.latitude-lattol, radar_site.latitude+lattol])

        lat, lon, _ = pm.ecef2geodetic(p['position'][0], p['position'][1], p['position'][2])

        ax.plot(lon, lat, transform=ccrs.Geodetic())
        ax.scatter(lon, lat, transform=ccrs.Geodetic())
        ax.scatter(radar_site.longitude, radar_site.latitude, marker='^', s=100, transform=ccrs.Geodetic())

        ax.set_title('{}, {}'.format(p['APT'],sid))

        fig.savefig('conjvalid_{:%Y%m%d_%H%M%S}.png'.format(p['APT']))
        plt.close(fig)

def output_file(starttime, endtime, radar, sid, deltime=60., zatol=25., latlontol=None, conjcoords='geo', filename='conjunctions.txt'):

    passes = conjunctions(starttime, endtime, radar, sid, deltime=deltime, lattol=lattol, lontol=lontol, conjcoords=conjcoords, return_all=True)

    with open(filename, 'w') as f:
        for p in passes:
            f.write('{:%Y-%m-%d %H:%M:%S}    {}\n'.format(p['APT'], p['mode']))

def main():
    st = dt.datetime(2017,10,1,0,0,0)
    et = dt.datetime(2017,11,1,0,0,0)

    # Spacecraft ID is the NORAD CAT ID
    # Can be found from the space-track.org satellite catalog (SATCAT)
    # sid = 39452		# SID for Swarm-A
    # sid = 25991     # SID for F15;
    sid = 39265     # SID for CASSIOPE/ePOP

    passes = conjunctions(st, et, 'RISRN', sid)
    print(passes)

    # validate(st, et, 'RISRN', sid)
    # output_file(st, et, 'RISRN', sid)


if __name__=='__main__':
    main()
