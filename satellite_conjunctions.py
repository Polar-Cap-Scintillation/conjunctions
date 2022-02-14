# satellite_conjunctions.py
# Generate a list of conjunctions between an AMISR and a particular satellite
# You must have an account at space-track.org (free and easy to set up)
# Requires the following non-standard packages:
#   - apexpy (https://apexpy.readthedocs.io/en/latest/readme.html)
#   - spacetrack (https://pythonhosted.org/spacetrack/)
#   - propagate_tle.py
# Also requires sqlalchemy and the appropriate AMISR SQL database
# Note: This takes a VERY long time to calculate tight conjunctions over many years of data

import numpy as np
import datetime as dt
from apexpy import Apex
import pymap3d as pm
from spacetrack import SpaceTrackClient
import spacetrack.operators as op
import sqlalchemy as db
from propagate_tle import propagate_tle

import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from space_track_credentials import *



def conjunctions(starttime, endtime, radar, sid, deltime=60., lattol=5., lontol=10., conjcoords='geo', return_all=False):
    # returns a list of conjuctions between a satellite (SID is the NORAD satellite ID) and ground-based AMISR when the radar is running
    # If return_all is set, a list of dictionaries containing times and satellite coordinates for each pass as well as
    #   the average pass time and radar mode is returned.  If not, the returned list only includes the average time of each pass.
    # deltime is time step between calculations of the satellite position
    # lattol is the latitudinal tolerance for a conjunction in degrees
    # lontol in the lognitudinal tolerance for a conjunction in degrees
    # The selections for deltime, lattol, & lontol can strongly impact the results.  The smaller deltime is, the finner
    #   the satellite position calculations will be, which reduces the change of "skipping over" a potential conjuction,
    #   but may increase the time it takes to run the code. Smaller values for lattol and lontol forces much tigher conjuctions,
    #   which will likely reduce the total number of conjunctions.  Ideally, deltime should be selected so that several
    #   steps (footprint velocity x deltime) fit in the tolerance box at the radar's location.
    # conjcoords specifies if conjunctions should be identified in geodetic ('geo') or magnetic ('mag') coordinates

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
    # deltime selection is important - large values may miss some conjunctions, small values will take a very long time
    # and appropriate value of deltime depends on the velocity (altitude) of the satellite
    unix_time_array = (starttime-dt.datetime.utcfromtimestamp(0)).total_seconds()+np.arange(0,(endtime-starttime).total_seconds(),deltime)
    time_array = np.array([dt.datetime.utcfromtimestamp(ut) for ut in unix_time_array])

    # Find index of epoch closest to each time in the time array
    # this works, but causes HUGE memory problems when arrays large (when trying to find conjunctions for years at a time)
    # closest_epoch = np.argmin(np.abs(unix_time_array[:,None] - unix_TLE_epoch[None,:]), axis=-1)
    closest_epoch = np.array([np.argmin(np.abs(ut-unix_TLE_epoch)) for ut in unix_time_array])


    sat_position = {'lat':[],'lon':[],'alt':[]}
    for i in range(len(TLE_epoch)):
        subset_times = time_array[closest_epoch==i]

        # calcualte satellite position using functions from TLE propgation script
        X, Y, Z = propagate_tle(subset_times,TLE_list[i])
        gdlat, gdlon, gdalt = pm.ecef2geodetic(X,Y,Z)
        # gdlon[gdlon>180.] -= 360.   # change range of longitude to -180-180
        sat_position['lat'].extend(gdlat)
        sat_position['lon'].extend(gdlon)
        sat_position['alt'].extend(gdalt)

    # convert to arrays for easier indexing
    for key, value in sat_position.items():
        sat_position[key] = np.array(value)



    # find all AMISR experiments files that fall within specified time
    engine = db.create_engine('sqlite:///'+amisr_dbfile_path+'{}_only_experiment_info.db'.format(radar.lower()))
    with engine.connect() as conn:
        exp = db.Table('experiment', db.MetaData(), autoload=True, autoload_with=engine)
        columns = [exp.columns.experiment, exp.columns.mode, exp.columns.start_year, exp.columns.start_month, exp.columns.start_time, exp.columns.end_time]
        unixstarttime = (starttime-dt.datetime.utcfromtimestamp(0)).total_seconds()
        unixendtime = (endtime-dt.datetime.utcfromtimestamp(0)).total_seconds()
        condition = db.and_(exp.columns.end_time>unixstarttime, exp.columns.start_time<unixendtime)
        query = db.select(columns).where(condition)
        exp_list = conn.execute(query).fetchall()

        site = db.Table('site', db.MetaData(), autoload=True, autoload_with=engine)
        query = db.select([site.columns.latitude, site.columns.longitude]).where(site.columns.shortname==radar.lower())
        radar_site = conn.execute(query).fetchall()[0]

    pass_list = []

    for exp in exp_list:

        # select indicies of time array durring this experiment
        tidx = np.argwhere((unix_time_array>=exp.start_time) & (unix_time_array<exp.end_time)).flatten()
        sat_time = unix_time_array[tidx]

        # find radar and satellite position in either geographic or magnetic coordinates, depending on conjcoord keyword
        if conjcoords == 'geo':
            radar_lat = radar_site.latitude
            radar_lon = radar_site.longitude
            sat_lat = sat_position['lat'][tidx]
            sat_lon = sat_position['lon'][tidx]
        elif conjcoords == 'mag':
            A = Apex(exp.start_year)
            radar_lat, radar_lon = A.geo2apex(radar_site.latitude, radar_site.longitude, 300.)
            sat_lat, sat_lon = A.geo2apex(sat_position['lat'][tidx], sat_position['lon'][tidx], sat_position['alt'][tidx])

        # find all times satellite position is within tolerance of the radar position
        conjunctions = np.all([np.abs(sat_lat-radar_lat)<lattol,np.abs(sat_lon-radar_lon)<lontol], axis=0)
        if not np.any(conjunctions):
            continue

        # seperate individual passes
        diffs = sat_time[conjunctions][1:]-sat_time[conjunctions][:-1]
        jumps = np.argwhere(diffs>deltime).flatten()

        if not jumps.size:
            # if only one pass in experiment
            passes = [{'time':sat_time[conjunctions],'lat':sat_lat[conjunctions],'lon':sat_lon[conjunctions]}]
        else:
            # first pass
            passes = [{'time':sat_time[conjunctions][:jumps[0]+1],'lat':sat_lat[conjunctions][:jumps[0]+1],'lon':sat_lon[conjunctions][:jumps[0]+1]}]
            # all intermediate passes
            for i in range(len(jumps)-1):
                passes.append({'time':sat_time[conjunctions][jumps[i]+1:jumps[i+1]+1],'lat':sat_lat[conjunctions][jumps[i]+1:jumps[i+1]+1],'lon':sat_lon[conjunctions][jumps[i]+1:jumps[i+1]+1]})
            # last pass
            passes.append({'time':sat_time[conjunctions][jumps[-1]+1:],'lat':sat_lat[conjunctions][jumps[-1]+1:],'lon':sat_lon[conjunctions][jumps[-1]+1:]})

        # find "average pass time"
        # this is the average between when the time when satellite latitude is closest to the radar
        #   latitude and the time when the satellite longitude is closest to the radar longitude
        # this is a little wonky, but is fast and avoids having to do geodecy calculations
        for p in passes:
            idx1 = np.argmin(np.abs(p['lat']-radar_lat))
            idx2 = np.argmin(np.abs(p['lon']-radar_lon))
            p['APT'] = dt.datetime.utcfromtimestamp((p['time'][idx1] + p['time'][idx2])/2)
            p['mode'] = exp.mode

        if return_all:
            pass_list.extend(passes)
        else:
            pass_list.extend([p['APT'] for p in passes])

    return pass_list


def validate(starttime, endtime, radar, sid, deltime=60., lattol=5., lontol=10., conjcoords='geo'):

    passes = conjunctions(starttime, endtime, radar, sid, deltime=deltime, lattol=lattol, lontol=lontol, conjcoords=conjcoords, return_all=True)

    # get site information from database
    engine = db.create_engine('sqlite:///'+amisr_dbfile_path+'{}_only_experiment_info.db'.format(radar.lower()))
    with engine.connect() as conn:
        site = db.Table('site', db.MetaData(), autoload=True, autoload_with=engine)
        query = db.select([site.columns.latitude, site.columns.longitude]).where(site.columns.shortname==radar.lower())
        radar_site = conn.execute(query).fetchall()[0]


    mapproj = ccrs.LambertConformal(central_latitude=radar_site.latitude, central_longitude=radar_site.longitude)

    for p in passes:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection=mapproj)
        ax.coastlines()
        ax.gridlines()
        ax.set_extent([radar_site.longitude-lontol, radar_site.longitude+lontol, radar_site.latitude-lattol, radar_site.latitude+lattol])

        ax.plot(p['lon'], p['lat'], transform=ccrs.Geodetic())
        ax.scatter(p['lon'], p['lat'], transform=ccrs.Geodetic())
        ax.scatter(radar_site.longitude, radar_site.latitude, marker='^', s=100, transform=ccrs.Geodetic())

        ax.set_title('{}, {}'.format(p['APT'],sid))

        fig.savefig('conjvalid_{:%Y%m%d_%H%M%S}.png'.format(p['APT']))
        plt.close(fig)

def output_file(starttime, endtime, radar, sid, deltime=60., lattol=5., lontol=10., conjcoords='geo', filename='conjunctions.txt'):

    passes = conjunctions(starttime, endtime, radar, sid, deltime=deltime, lattol=lattol, lontol=lontol, conjcoords=conjcoords, return_all=True)

    with open(filename, 'w') as f:
        for p in passes:
            f.write('{:%Y-%m-%d %H:%M:%S}    {}\n'.format(p['APT'], p['mode']))

def main():
    st = dt.datetime(2013,9,1,0,0,0)
    et = dt.datetime(2013,10,1,0,0,0)

    # Spacecraft ID is the NORAD CAT ID
    # Can be found from the space-track.org satellite catalog (SATCAT)
    # sid = 39452		# SID for Swarm-A
    # sid = 25991     # SID for F15;
    sid = 39265     # SID for CASSIOPE/ePOP

    passes = conjunctions(st, et, 'RISRN', sid)
    print(passes)

    validate(st, et, 'RISRN', sid)
    output_file(st, et, 'RISRN', sid)


if __name__=='__main__':
    main()
