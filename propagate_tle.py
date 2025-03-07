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
from pathlib import Path

# import SpaceTrack username and password
from space_track_credentials import *

class TLEHandler(object):
    def __init__(self, sid):

        self.load_tle(sid)
        #self.create_tle_library(sid)

    def load_tle(self, sid):
        # load TLE library for a single sattelite from local
        # if satellite does not exist, download from spacetrack

        tle_file = Path(local_tle_library, f'{sid}.txt')
        #print(tle_file)

        # parse output into list of distinct TLEs
        #try:
        with open(tle_file, 'r') as f:
            output = f.read()
        #except FileNotFoundError:
        #    self.download_tle(sid)
        #    with open(tle_file, 'r') as f:
        #        output = f.read()

        split = output.splitlines()
        self.TLE_list = [[split[2*i],split[2*i+1]] for i in range(len(split)//2)]

        # extract epoch from each TLE
        # Do this as np.datetime64 instead of datetime?  May be faster?
        self.TLE_epoch = [dt.datetime.strptime(tle[0][18:23],'%y%j')+dt.timedelta(days=float(tle[0][23:32])) for tle in self.TLE_list]

#    def download_tle(self, sid):
#        # retrieve all TLEs for satellite from space-track.org
#        # The space-track.org API limits the number of calls you can make, so it's better to retreive the
#        #   entire library for one satellite and save it locally than try to collect individual
#        #   TLEs as needed.
#
#        print(f'Downloading TLE for {sid} from spacetrack.org')
#        st = SpaceTrackClient(identity=ST_USERNAME, password=ST_PASSWORD, base_url='https://for-testing-only.space-track.org')
#        output = st.gp_history(norad_cat_id=sid, orderby='epoch asc', format='tle')
#
#        tle_file = Path(local_tle_library, f'{sid}.txt')
#        with open(tle_file, 'w') as f:
#            f.write(output)
#        print('...Done')
#        #print(output)
#
#
#
#    def create_tle_library(self, sid):
#        # retrieve TLEs for entire period from space-track.org
#        # The space-track.org API limits the number of calls you can make, so it's better to retreive the
#        #   entire library for one satellite and save it as a class attribute than try to collect indivitual
#        #   TLEs as needed.
#        st = SpaceTrackClient(identity=ST_USERNAME, password=ST_PASSWORD, base_url='https://for-testing-only.space-track.org')
#        output = st.gp_history(norad_cat_id=sid, orderby='epoch asc', format='tle')
#
#        print(output)
#
#        # parse output into list of distinct TLEs
#        split = output.splitlines()
#        self.TLE_list = [[split[2*i],split[2*i+1]] for i in range(len(split)//2)]
#
#        # extract epoch from each TLE
#        self.TLE_epoch = [dt.datetime.strptime(tle[0][18:23],'%y%j')+dt.timedelta(days=float(tle[0][23:32])) for tle in self.TLE_list]

    def sat_position(self, time_array):

        # Find index of epoch closest to each time in the time array
        unix_time_array = np.array(time_array).astype('datetime64')
        unix_TLE_epoch = np.array(self.TLE_epoch).astype('datetime64')
        closest_epoch_idx = np.array([np.argmin(np.abs(ut-unix_TLE_epoch)) for ut in unix_time_array])

        # Somewhere in here check that epoch is within 10 days and if not, update TLE library from spacetrack.org?

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


def example():

    # Create a seperate file called space_track_credentials.py and add only your spacetrack credentials as shown:
    # ST_USERNAME=''
    # ST_PASSWORD=''

    # NORAD satellite ID
    #sat_id = 44628   # TLE for ICON
    sat_id = 39452

    # Create TLEHandler object
    tle = TLEHandler(sat_id)

    # Create array of desired times
    time_array = np.array([dt.datetime(2020,1,1,0,0,0)+dt.timedelta(seconds=60.*m) for m in range(60)])

    # Call tle.sat_position to calculate the satelite position at each time
    sat_position = tle.sat_position(time_array).T

    # position returned in ECEF coordinates
    print(sat_position)


def create_tle_library():

    # Swarm A, B, C; LLITED A, B
    #sid = [39451, 39452, 39453, 56219, 56220]

   # sat_id = [10684,
   #     10893,
   #     11054,
   #     11141,
   #     11690,
   #     11783,
   #     14189,
   #     15039,
   #     15271,
   #     16129,
   #     20061,
   #     19802,
   #     20830,
   #     20185,
   #     20361,
   #     20452,
   #     20302,
   #     20533,
   #     20724,
   #     22446,
   #     20959,
   #     21552,
   #     21890,
   #     22014,
   #     22108,
   #     21930,
   #     22275,
   #     24320,
   #     22581,
   #     22231,
   #     23833,
   #     22877,
   #     22779,
   #     23027,
   #     22657,
   #     25030,
   #     22700,
   #     23953,
   #     26605,
   #     24876,
   # sat_id = [26407,
   #     27704,
   #     25933,
   #     28129,
   #     32711,
   # sat_id = [34661,
   #     35752,
   #     26360,
   #     29486,
   #     28874,
   # sat_id = [26690,
   #     32260,
   #     27663,
   #     32384,
   #     29601,
   # sat_id = [28190,
   #     28361,
   #     28474,
   #     36585,
   #     37753,
   #     39533,
   #     38833,
   #     39166,
   #     39741,
   #     40105,
   #     40294,
   #     41328,
   #     40534,
   #     40730,
   #     41019,
   #     43873,
   #     44506,
   #     45854,
   #     46826,
   #     48859,
    sid = [55268,
        62339]
#        28112,
#        26987,
#        28509,
#        28916,
#        28915,
#        29672,
#        29670,
#        29671,
#        32277,
#        32276,
#        32275,
#        32393,
#        32394,
#        32395,
#        33378,
#        33379,
#        33380,
#        33466,
#        33468,
#        33467,
#        36111,
#        36400,
#        36402,
#        36112,
#        36113,
#        36401,
#        37139,
#        37138,
#        37137,
#        37829,
#        37869,
#        37867,
#        37868,
#        37938,
#        39155,
#        21216,
#        21217,
#        21218,
#        22057,
#        22513,
#        23045,
#        22514,
#        23043,
#        23044,
#        23398,
#        23396,
#        23397,
#        23511,
#        23512,
#        23203,
#        21853,
#        21854,
#        23205,
#        21855,
#        22058,
#        22512,
#        22056,
#        23204,
#        23736,
#        23513,
#        23735,
#        25593,
#        23620,
#        23621,
#        23734,
#        26566,
#        25594,
#        23622,
#        25595,
#        26564,
#        26565,
#        26988,
#        26989,
#        27617,
#        27619,
#        27618,
#        28113,
#        28114,
#        28508,
#        28510,
#        28917,
#        37372,
#        40315,
#        57517,
#        46805,
#        52984,
#        54031,
#        41330,
#        42939,
#        41554,
#        39620,
#        40001,
#        43508,
#        43687,
#        44299,
#        44850,
#        45358,
#        54377,
#        28922,
#        32781,
#        37846,
#        37847,
#        38857,
#        38858,
#        40128,
#        40129,
#        40544,
#        40545,
#        40889,
#        40890,
#        41859,
#        41175,
#        41174,
#        41550,
#        41549,
#        41860,
#        41861,
#        41862,
#        43055,
#        43056,
#        43057,
#        43058,
#        43566,
#        43567,
#        43564,
#        43565,
#        49809,
#        49810,
#        59598,
#        61183,
#        59600,
#        61182,
#        31115,
#        34779,
#        36287,
#        36590,
#        36828,
#        37210,
#        37256,
#        37384,
#        37763,
#        37948,
#        38091,
#        38250,
#        38251,
#        38774,
#        38775,
#        38953,
#        41434,
#        41586,
#        43539,
#        44231,
#        40549,
#        40749,
#        40748,
#        40938,
#        41315,
#        43001,
#        43002,
#        43107,
#        43108,
#        43207,
#        43208,
#        43245,
#        43246,
#        43581,
#        43582,
#        43602,
#        43603,
#        43622,
#        43623,
#        43647,
#        43648,
#        43683,
#        43706,
#        43707,
#        44204,
#        44337,
#        44542,
#        44543,
#        44709,
#        44793,
#        44794,
#        44864,
#        44865,
#        45344,
#        45807,
#        56564,
#        58654,
#        58655,
#        61186,
#        61187,
#        37158,
#        42738,
#        42917,
#        42965,
#        49336,
#        39199,
#        39635,
#        40269,
#        40547,
#        41241,
#        41384,
#        41469,
#        43286,
#        56759]


#    for i in range(int(len(sat_id)/5)):
#        sid = sat_id[i*5:(i+1)*5].copy()
#        print(sid)
        
    sid.sort()

    with SpaceTrackClient(identity=ST_USERNAME, password=ST_PASSWORD, base_url='https://for-testing-only.space-track.org') as client:
        output = client.gp_history(norad_cat_id=sid, orderby=['norad_cat_id','epoch asc'], format='tle')

    for i in range(len(sid)-1):
        part = output.partition(f'1 {sid[i+1]}U')
    #for p in part:
    #    print(p)

        tle_file = Path(local_tle_library, f'{sid[i]}.txt')
        with open(tle_file, 'w') as f:
            f.write(part[0])
        output = part[1]+part[2]

    tle_file = Path(local_tle_library, f'{sid[-1]}.txt')
    with open(tle_file, 'w') as f:
        f.write(output)
        #print('...Done')

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
#    main()
    example()
#    create_tle_library()
