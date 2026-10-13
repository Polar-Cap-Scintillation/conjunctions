# tle.py

# This module manages finding the correct TLEs (Two Line Elements) and using
#   them to compute the position of a satellite.
# There are a variety of services for providing the latest TLEs, but when using
#   this code to calculate historical conjunctions it is better to download
#   full files of TLEs and reference them locally. It violates the API use
#   agreement to most of these to make lots of repeted queries to download the
#   same information and you'll get your IP banned.
# Lists of TLEs are typically distributed as text files, which are extremely slow
#   to read for high volumes.  This in particular has become a problem with the
#   recent explosion of comercial LEO satellites, making single year files GB in
#   volume.  To improve this, this module expects you to first create a SQL
#   database from the TLE text file, which is faster than looking up TLEs from
#   text files.  To create the sql database, run the following command line:
#       $ create-tle-db tle1.txt tle2.txt tle3.txt
#   Source files can be downloaded from the space-track.org bulk data store:
#   https://ln5.sync.com/dl/afd354190/c5cd2q72-a5qjzp4q-nbjdiqkr-cenajuqu
# The sgp4 library uses Julian dates mostly to keep track of epochs, so that is
#   what's used for most internal functionality of this module.
# Most TLE propigation stuff is now handled internally by the sgp4 and the skyfield
#   packages.  These references may still be useful though, specifically for the
#   coordinate transformation from the output of the SGP4 calculations (TEME?)
#   into something useful for geospace applications


import os
import numpy as np
import datetime as dt
import pymap3d as pm
from tqdm import tqdm
import argparse

from sgp4.api import Satrec, WGS72, jday
from sgp4 import exporter
from sgp4.conveniences import sat_epoch_datetime
from skyfield.api import EarthSatellite, load, wgs84

import sqlalchemy
from sqlalchemy import Column, Integer, Float
from sqlalchemy import create_engine, desc
from sqlalchemy.orm import declarative_base, Session



# Define class that will be used to create TLE SQL database
Base = declarative_base()

class TLE(Base):
    __tablename__ = 'tle'
    id = Column(Integer, primary_key = True)
    satnum = Column(Integer, nullable=False, index=True)
    epoch = Column(Float, nullable=False)
    bstar = Column(Float, nullable=False)
    ndot = Column(Float, nullable=False)
    nddot = Column(Float, nullable=False)
    ecco = Column(Float, nullable=False)
    argpo = Column(Float, nullable=False)
    inclo = Column(Float, nullable=False)
    mo = Column(Float, nullable=False)
    no_kozai = Column(Float, nullable=False)
    nodeo = Column(Float, nullable=False)



def create_sql_database(source_files, dbfile='tle.db', satids=None):
    """
    Generate TLE SQL database from source text files.
    Text files can be downloaded in bulk from:
    https://ln5.sync.com/dl/afd354190/c5cd2q72-a5qjzp4q-nbjdiqkr-cenajuqu
    (This is the space-track.org bulk download site, not something sketchy.)
    """

    if os.path.exists(dbfile):
        raise FileExistsError(f'File {dbfile} already exists!')

    engine = create_engine(f"sqlite:///{dbfile}", echo=False)
    Base.metadata.create_all(engine)
    
    numfiles = len(source_files)    # Number of input files for progress tracking

    with Session(engine) as session:

        i = 0       # unique id counter
        for fi, srcfile in enumerate(source_files):
            print(f'[{fi+1}/{numfiles}] {srcfile}')

            # Count number of lines in file
            with open(srcfile, 'rb') as f:
                #num_lines = sum(1 for line in f)
                num_lines = sum(1 for _ in f)
 
            with open(srcfile, 'r') as f:
                for l1 in tqdm(f, total=num_lines/2):
                    l2 = f.readline()

                    # formulate sgp4 object
                    try:
                        sgp4obj = Satrec.twoline2rv(l1, l2)
                    except ValueError as e:
                        continue
                    
                    # If satellite list provided, skip satellites not in it
                    if satids:
                        if sgp4obj.satnum not in satids:
                            continue
       
                    # Add TLE to database
                    sqlobj = sgp2sql(sgp4obj, i)
                    session.add(sqlobj)

                    
                    # Commit every 1000 rows to not use too much memory???
                    if (i+1) % 10000 == 0:
                        session.commit()
                        session.expunge_all()

                    i += 1


            print('Committing to SQL database ...')
            session.commit()
            session.expunge_all()


def cli_create_database():
    """
    CLI fuction to generate TLE SQL database from source files
    """
    parser = argparse.ArgumentParser(prog='create-tle-db',
                                     description='Generate TLE database from source files.')
    parser.add_argument('input', nargs='+', help='input source txt files')
    parser.add_argument('-o', '--output', default='tle.db', help='output db file')
    parser.add_argument('-s', '--satidlist', help='file of specific satids to include (defalt is incude all)')
    args = parser.parse_args()

    if args.satidlist:
        with open(args.satidlist) as f:
            satids = f.readlines()
        satids = [int(s) for s in satids]
    else:
        satids = None

    create_sql_database(args.input, dbfile=args.output, satids=satids)



def sgp2sql(sgp4obj, i):
    """
    Convert a sgp4.Satrec objet to a TLE obect for SQL databae.
    """

    # Calculate Julian day epoch
    jd_epoch = sgp4obj.jdsatepoch + sgp4obj.jdsatepochF

    tleobj = TLE(id=i,
                 satnum   = sgp4obj.satnum,
                 epoch    = jd_epoch,
                 bstar    = sgp4obj.bstar,
                 ndot     = sgp4obj.ndot,
                 nddot    = sgp4obj.nddot,
                 ecco     = sgp4obj.ecco,
                 argpo    = sgp4obj.argpo,
                 inclo    = sgp4obj.inclo,
                 mo       = sgp4obj.mo,
                 no_kozai = sgp4obj.no_kozai,
                 nodeo    = sgp4obj.nodeo)

    return tleobj
        

def sql2sgp(sqlobj):
    """
    Convert a SQL database TLE obect to a sgp4.Satrec object.
    """

    #To compute the “epoch” argument, take the epoch’s Julian date and subtract 2433281.5 days.
    epoch = sqlobj.epoch - 2433281.5

    sgp4obj = Satrec()
    sgp4obj.sgp4init(
                     WGS72,                 # gravity model
                     'i',                   # 'a' = old AFSPC mode, 'i' = improved mode
                     sqlobj.satnum,         # satnum: Satellite number
                     epoch,                 # epoch: days since 1949 December 31 00:00 UT
                     sqlobj.bstar,          # bstar: drag coefficient (1/earth radii)
                     sqlobj.ndot,           # ndot: ballistic coefficient (radians/minute^2)
                     sqlobj.nddot,          # nddot: mean motion 2nd derivative (radians/minute^3)
                     sqlobj.ecco,           # ecco: eccentricity
                     sqlobj.argpo,          # argpo: argument of perigee (radians 0..2pi)
                     sqlobj.inclo,          # inclo: inclination (radians 0..pi)
                     sqlobj.mo,             # mo: mean anomaly (radians 0..2pi)
                     sqlobj.no_kozai,       # no_kozai: mean motion (radians/minute)
                     sqlobj.nodeo,          # nodeo: R.A. of ascending node (radians 0..2pi)
                    )

    return sgp4obj


class TLEHandler(object):
    def __init__(self, dbfile='tle.db'):

        self.load_db(dbfile)


    def load_db(self, dbfile):
        """
        Load SQL database of historical TLEs
        """

        engine = create_engine(f"sqlite:///{dbfile}", echo=False)
        self.session = Session(engine)


    def select_tles(self, sat_cat, time0, time1, jd_epoch=False):
        """
        Find all TLEs for a certain satellite between two times, plus the nearest
        one before and after the time interval.  This should create a comprehensive
        list of all TLEs that might be needed for propigating satelite motion between
        those two times.
        """

        sat_cat = int(sat_cat)

        # Calculate Julian day of start and end time
        jd0, fr0 = jday(time0.year, time0.month, time0.day, time0.hour, time0.minute, time0.second)
        jdtime0 = jd0 + fr0
        jd1, fr1 = jday(time1.year, time1.month, time1.day, time1.hour, time1.minute, time1.second)
        jdtime1 = jd1 + fr1


        # All TLEs between first and last time
        conditions = sqlalchemy.and_(TLE.satnum==sat_cat,
                                     TLE.epoch>=jdtime0,
                                     TLE.epoch<=jdtime1)
        tle_between = self.session.query(TLE).filter(conditions).order_by(TLE.epoch).all()

        # Last epoch before first time
        conditions = sqlalchemy.and_(TLE.satnum==sat_cat,
                                      TLE.epoch<jdtime0)
        tle_first = self.session.query(TLE).filter(conditions).order_by(desc(TLE.epoch)).first()

        # First epoch after last time
        conditions = sqlalchemy.and_(TLE.satnum==sat_cat,
                                      TLE.epoch>jdtime1)
        tle_last = self.session.query(TLE).filter(conditions).order_by(TLE.epoch).first()

        # Create full list
        tle_list = [tle_first] + tle_between + [tle_last]

        # Extract array of epoch times
        if jd_epoch:
            # Extract as Julian days
            epoch_list = [t.epoch for t in tle_list]
        else:
            # Extract as datetime objects
            epoch_list = [sat_epoch_datetime(sql2sgp(t)) for t in tle_list]

        # Return epochs as datetime objects UNLESS JD flag set
        return epoch_list, tle_list


    def sat_position(self, sat_cat, time_array):
        """
        Calculate satellite position at all points in given time array
        Note: This function only works for sequential times
        """

        # Calculate Julian times from time array
        jd_array = np.sum([jday(t.year, t.month, t.day, t.hour, t.minute, t.second) for t in time_array], axis=1)

        # Get all TLEs that fall in time range
        epoch_list, tle_list = self.select_tles(sat_cat, time_array[0], time_array[-1], jd_epoch=True)

        # Find index of epoch closest to each time in the time array
        closest_epoch_idx = np.array([np.argmin(np.abs(t-epoch_list)) for t in jd_array])

        # Somewhere in here check that epoch is within 10 days and if not, raise warning

        # Cycle through epochs and calculate position for satellite for each one when that epoch is the closest time
        sat_position = np.empty((3,0))

        for i in np.unique(closest_epoch_idx):
            # Select subset of times closest to a particular epoch
            subset_times = np.array(time_array)[closest_epoch_idx==i]
            # TLE for that epoch
            sgp4obj = sql2sgp(tle_list[i])
            line1, line2 = exporter.export_tle(sgp4obj)

            # calcualte satellite position using functions from TLE propgation script
            X, Y, Z = propagate_tle2(subset_times, line1, line2)
            sat_position = np.append(sat_position, np.array([X, Y, Z]), axis=1)

        return sat_position




# This can probably be made more efficient by doing the satelite propagation directly
#   with the sgp4 pacakge and only using skyfield for coordinate conversion.
#   This is sort of done in propagate_tle3(), but it's not quite right.
def propagate_tle2(time0, tleline1, tleline2):

    ts = load.timescale()
    satellite = EarthSatellite(tleline1, tleline2)
   
    tmp = np.array([[t.year, t.month, t.day, t.hour, t.minute, t.second] for t in time0])
    tstmp = ts.utc(tmp[:,0], tmp[:,1], tmp[:,2], tmp[:,3], tmp[:,4], tmp[:,5])

    geocentric = satellite.at(tstmp)

    posobj = wgs84.geographic_position_of(geocentric)
   
    return posobj.itrs_xyz.m



#def propagate_tle3(time0, tleline1, tleline2):
#
#    # initialize tle object
#    tle = Satrec.twoline2rv(tleline1, tleline2)
#
#    out = [jday(t.year, t.month, t.day, t.hour, t.minute, t.second) for t in time0]
#    jd, fr = np.array(out, order='F').T     # order='F' needed for some kind of error not contiguous in C error when passed into sgp4_array() ???
#    e, r, v = tle.sgp4_array(jd, fr)
#    #print(r.shape, v.shape)
#
#
#    #posobj = TEME(r)
#
#    r,v = sgp4lib.TEME_to_ITRF(jd+fr, r, v)
#    print(r.shape)
#
#    fig, ax = plt.subplots(subplot_kw={'projection':'3d'})
#    #ax.plot(position[0], position[1], position[2])
#    ax.plot(posobj.itrs_xyz.km[0], posobj.itrs_xyz.km[1], posobj.itrs_xyz.km[2])
#    plt.show()




