# amisr_lookup.py
# Look up the experiment an AMISR radar was running in durring a particular period,
#   including mode/experiment information.
# Requires sqlalchemy and the appropriate AMISR SQL database

import datetime as dt
import sqlalchemy as db
import pathlib
import re

import sys
sys.path.append(pathlib.Path(__file__).parent.resolve())
from space_track_credentials import amisr_dbfile_path


class AMISR_lookup(object):

    def __init__(self, radar):
        self.radar = radar
        # find all AMISR experiments files that fall within specified time
        self.engine = db.create_engine('sqlite:///'+amisr_dbfile_path+'{}_only_experiment_info.db'.format(radar))
        # with engine.connect() as conn:
        self.conn = self.engine.connect()
        self.exp = db.Table('experiment', db.MetaData(), autoload=True, autoload_with=self.engine)
        self.params = [self.exp.columns.experiment, self.exp.columns.mode, self.exp.columns.start_time, self.exp.columns.end_time]

    def find_experiments(self, starttime, endtime):
        unixstarttime = (starttime-dt.datetime.utcfromtimestamp(0)).total_seconds()
        unixendtime = (endtime-dt.datetime.utcfromtimestamp(0)).total_seconds()

        # queary AMISR database for experiments in this time frame
        condition = db.and_(self.exp.columns.end_time>unixstarttime, self.exp.columns.start_time<unixendtime)
        query = db.select(self.params).where(condition)
        exp_list = self.conn.execute(query).fetchall()

        return [{'experiment_number':exp.experiment, 'mode':exp.mode, 'start_time':exp.start_time, 'end_time':exp.end_time} for exp in exp_list]

    def check_master_exp(self, exp, procdir='/Volumes/AMISR_PROCESSED/processed_data'):
        exp_path = self.experiment_path(exp, procdir)

        expdir = pathlib.Path(procdir, self.radar.upper(), exp['experiment_number'][0:4], exp['experiment_number'][4:6], exp['mode'], exp['experiment_number'])


    def experiment_path(self, exp, procdir='/Volumes/AMISR_PROCESSED/processed_data'):

        expdir = pathlib.Path(procdir, self.radar.upper(), exp['experiment_number'][0:4], exp['experiment_number'][4:6], exp['mode'], exp['experiment_number'])
        print(expdir)

        # check if path exists before returning
        if not expdir.exists():
            expdir = expdir.with_name(exp['experiment_number']+'.done')

        if expdir.exists():
            return expdir
        else:
            return None

        # try:
        #     path = '/Volumes/AMISR_PROCESSED/processed_data/RISR-N/{}/{}/{}/{}'.format(exp.experiment[0:4],exp.experiment[4:6],exp.mode,exp.experiment)
        #     datafiles = [f for f in os.listdir(path) if (f.endswith('h5') and '_lp_' in f)]
        # except OSError:
        #     try:
        #         path = '/Volumes/AMISR_PROCESSED/processed_data/RISR-N/{}/{}/{}/{}.done'.format(exp.experiment[0:4],exp.experiment[4:6],exp.mode,exp.experiment)
        #         datafiles = [f for f in os.listdir(path) if (f.endswith('h5') and '_lp_' in f)]
        #     except OSError:
        #         continue



    def select_experiment_file(self, experiment_path, type='lp'):

        if not experiment_path:
            return None

        exp_num = experiment_path.name
        datafiles = [f.name for f in experiment_path.glob('*.h5')]
        # print(datafiles)



        # create list of "normal" fitted files, fitcal if available, then cal, then all others
        calfiles = [f for f in datafiles if re.match(r'{}_{}_\d+min-fitcal.h5'.format(exp_num, type), f)]
        if not calfiles:
            calfiles = [f for f in datafiles if re.match(r'{}_{}_\d+min-cal.h5'.format(exp_num, type), f)]
        if not calfiles:
            calfiles = [f for f in datafiles if re.match(r'{}_{}_\d+min.h5'.format(exp_num, type), f)]
        if not calfiles:
            return None


        # Get file with minimum time resolution
        timeres = [int(re.search(r'\d+min', f).group()[:-3]) for f in calfiles]
        datafile = calfiles[timeres.index(min(timeres))]


        return experiment_path.joinpath(datafile)



        # # choose file with shortest time resolution
        # try:
        #     timeres = [int(re.search(r'\d+min', f).group()[:-3]) for f in calfiles]
        #     datafile = calfiles[timeres.index(inttime)]
        #     # fileres = inttime
        #     postint = False
        # except ValueError:
        #     try:
        #         datafile = calfiles[timeres.index(min(timeres))]
        #         # fileres = min(timeres)
        #         postint = True
        #     except ValueError:
        #         continue
        #
        # print(os.path.join(path,datafile))
        # num_files += 1



    def site_coords(self):
        site = db.Table('site', db.MetaData(), autoload=True, autoload_with=self.engine)
        query = db.select([site.columns.latitude, site.columns.longitude]).where(site.columns.shortname==self.radar.lower())
        radar_site = self.conn.execute(query).fetchall()[0]
        return radar_site


    def __del__(self):
        self.conn.close()
