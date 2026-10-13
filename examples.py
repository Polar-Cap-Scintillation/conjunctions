# examples.py
# Several examples of different ways to use this code

import os
import glob
import datetime as dt
import numpy as np
import pymap3d as pm

import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from satgroundconj import tle
from satgroundconj.conjunction import SatConj

from gnssutils.map_prn import prn2norad
from gnssutils import calc_ipp

def known_tle():
    """
    TLE is already known, either provided by the user or obtained in some other part of the code
    """
    print('known_tle')

    TLE = ['1     1U          18350.30892361  .00001123  00000-0  66525-4 0   109',
           '2     1  85.0373 178.2871 0002550 225.5672 175.5175 15.21584957    13']

    times = np.array([dt.datetime(2018,12,17,0,0,0)+dt.timedelta(hours=h) for h in range(24)])

    X, Y, Z = tle.propagate_tle(times,TLE)
    gdlat, gdlon, gdalt = pm.ecef2geodetic(X,Y,Z)

    print('{:^20}{:^10}{:^10}{:^10}'.format('Time','GLAT','GLON','GALT'))
    for t, lat, lon, alt in zip(times,gdlat,gdlon,gdalt):
        print('{}{:10.2f}{:10.2f}{:10.2f}'.format(t, lat, lon, alt))


def generate_database():
    """
    Generate TLE SQL database from source files
    """
    
    #tle_dir = '/Users/e30737/Desktop/Data/TLE/srctxt'
    tle_dir = '/Volumes/Janeway/TLE/srctxt'
    #tle_files = [os.path.join(tle_dir, f) for f in os.listdir(tle_dir)]
    #tle_files = [os.path.join(tle_dir, f) for f in sorted(glob.glob('tle2004_*.txt', root_dir=tle_dir))]
    #print(tle_files)
    tle_files = [os.path.join(tle_dir, 'tle2015.txt')]
    tle.create_sql_database(tle_files, dbfile='/Users/e30737/Desktop/Data/TLE/test9.db', satids=[39452,39451,39453])


def database_tle1():

    sat_id = 39452   # Swarm A
    #sat_id = 25544   # ISS
    starttime = dt.datetime(2015,2,1)
    endtime = dt.datetime(2015,2,10)

    dbfile = '/Users/e30737/Desktop/Data/TLE/test9.db'

    tlelib = tle.TLEHandler(dbfile=dbfile)
    epochs, tles = tlelib.select_tles(sat_id, starttime, endtime)
    print(epochs)


def database_tle():
    """
    A NORAD ID and times are given and the appropriate TLEs must be retrieved from the TLE database
    """
    print('database_tle')

    sat_id = 39452   # Swarm A
    #time_list = [dt.datetime(2015,2,18,5,50,0)+dt.timedelta(minutes=i) for i in range(300)]
    #time_list = [dt.datetime(2015,9,5,11,30,0, tzinfo=dt.timezone.utc)+dt.timedelta(minutes=i) for i in range(300)]
    time_list = [dt.datetime(2015,9,5,11,37,0)+dt.timedelta(minutes=i) for i in range(15)]

    #dbfile = '/Users/e30737/Desktop/Data/TLE/tle.db'

    #sat_id = 25544   # ISS
    #time_list = [dt.datetime(2000,2,10,13,25,0)+dt.timedelta(minutes=i) for i in range(6000)]

    dbfile = '/Users/e30737/Desktop/Data/TLE/test9.db'

    tlelib = tle.TLEHandler(dbfile=dbfile)
    X, Y, Z = tlelib.sat_position(sat_id, time_list)
    glat, glon, galt = pm.ecef2geodetic(X, Y, Z)

    proj = ccrs.AzimuthalEquidistant(central_longitude=-115, central_latitude=60)
    fig, ax = plt.subplots(subplot_kw=dict(projection=proj))
    ax.coastlines()
    ax.gridlines()
    ax.set_extent([-130, -100, 30, 80])

    ax.scatter(glon, glat, transform=ccrs.PlateCarree())

    plt.show()



def conjunction():
    """
    Calculate conjunction with a single satellite
    """
    print('conjunction')

    # Double check these
    pfrr_glat = 65.5
    pfrr_glon = -147.7
    pfrr_galt = 0.
    site = [pfrr_glat, pfrr_glon, pfrr_galt]
    
    sat_id = 39452   # Swarm A
    starttime = dt.datetime(2016,2,1)
    endtime = dt.datetime(2016,2,29)
    
    dbfile = '/Users/e30737/Desktop/Data/TLE/tle.db'

    conj = SatConj(site, sat_id, deltime=60., conjtype='zenith', tolerance=25., tledb='/Users/e30737/Desktop/Data/TLE/test8.db')

    conj_list = conj.conjunctions(starttime, endtime)
    
    print('number of conjunctions:', len(conj_list))
    
    print('{:^20}{:^10}{:^10}{:^10}'.format('Time','GLAT','GLON','GALT'))
    for c in conj_list:
        # Only show first point in conjunction for brevity
        t = c['time'][0]
        X, Y, Z = c['position'][0]
        lat, lon, alt = pm.ecef2geodetic(X, Y, Z)
        print('{}{:10.2f}{:10.2f}{:10.2f}'.format(t, lat, lon, alt))
    


def gnss_conjunction():
    """
    Calculate conjunctions wiht any GNSS satellite
    """
    print('gnss_conjunction')

    # Double check these
    pfrr_glat = 65.5
    pfrr_glon = -147.7
    pfrr_galt = 0.
    site = [pfrr_glat, pfrr_glon, pfrr_galt]
    
    starttime = dt.datetime(2016,2,1)
    endtime = dt.datetime(2016,2,2)
    
    # Find conjunctions between GPS satellites and PFRR
    for gnum in range(1,32):
        prn = f'G{gnum:02d}'
        sat_id = prn2norad(prn, starttime)
        if not sat_id:
            continue
        
        conj = SatConj(site, sat_id, deltime=60., conjtype='zenith', tolerance=25., tledb='/Users/e30737/Desktop/Data/TLE/test8.db')
    
        conj_list = conj.conjunctions(starttime, endtime)
        
        print(prn, len(conj_list))
        
        #print('{:^20}{:^10}{:^10}{:^10}'.format('Time','GLAT','GLON','GALT'))
        for c in conj_list:
            # Only show first point in conjunction for brevity
            t = c['time'][0]
            X, Y, Z = c['position'][0]
            ipp_glat, ipp_glon = calc_ipp(site, [X,Y,Z], satcoords='ecef', height=300.)
            print('{}{:10.2f}{:10.2f}'.format(t, ipp_glat, ipp_glon))



if __name__=='__main__':
    #generate_database()
    #known_tle()
    database_tle()
    #conjunction()
    #gnss_conjunction()

