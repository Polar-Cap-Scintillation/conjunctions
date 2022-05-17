# amisr_lookup.py
# Look up the experiment an AMISR radar was running in durring a particular period,
#   including mode/experiment information.
# Requires sqlalchemy and the appropriate AMISR SQL database

import datetime as dt
import sqlalchemy as db
import pathlib

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

    def site_coords(self):
        site = db.Table('site', db.MetaData(), autoload=True, autoload_with=self.engine)
        query = db.select([site.columns.latitude, site.columns.longitude]).where(site.columns.shortname==self.radar.lower())
        radar_site = self.conn.execute(query).fetchall()[0]
        return radar_site


    def __del__(self):
        self.conn.close()
