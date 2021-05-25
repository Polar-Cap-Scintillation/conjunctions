import os
import sys
from sqlalchemy import Column, ForeignKey, Integer, String, Boolean, Float
from sqlalchemy.ext.declarative import declarative_base
 
Base = declarative_base()

class Experiment(Base):
    __tablename__ = 'experiment'
    experiment = Column(String(128), primary_key=True)   # typically, yyyymmdd.xxx
    mode = Column(String(128), ForeignKey('radarmode.name'), nullable=False)  # mode the radar was running e.g. WorldDay40.v01
    start_year  = Column(Integer, nullable=False)
    end_year    = Column(Integer, nullable=False)
    start_month = Column(Integer, nullable=False)
    end_month   = Column(Integer, nullable=False)
    start_day   = Column(Integer, nullable=False)
    end_day     = Column(Integer, nullable=False)
    start_time  = Column(Integer, nullable=False)  # unix time
    end_time    = Column(Integer, nullable=False)  # unix time
    total_seconds    = Column(Integer, nullable=False) # total time experiment ran for
    start_filenumber = Column(Integer, nullable=False)  # filenumbers as output by radac e.g. d00011234
    end_filenumber   = Column(Integer, nullable=False)  # filenumbers as output by radac e.g. d00011234
    path = Column(String(1024), nullable=False) # path to experiment as mounted on Appleton
    site = Column(String(128), ForeignKey('site.name'),nullable=False)
    filesize = Column(Float, nullable=False) # total size in MB of experiments
    aeurx = Column(Float, nullable=True) # median aeurx in the experiment
    aeutx = Column(Float, nullable=True) # median aeutx in the experiment
    txpower = Column(Float, nullable=True) # median txpower in the experiment


class Site(Base):
    __tablename__ = 'site'
    name = Column(String(128), nullable=False)
    shortname = Column(String(128), primary_key=True)
    code = Column(Integer, nullable=False)
    altitude = Column(Float, nullable=False)
    latitude = Column(Float, nullable=False)
    longitude = Column(Float, nullable=False)


class RadarMode(Base):
    __tablename__ = 'radarmode'
    name = Column(String(128), primary_key=True) # something like WorldDay40.v01
    beams = Column(String(2048), nullable=False) # comma separated sorted string of beam codes e.g. '64047,64517'
    num_beams = Column(Integer, nullable=False)
    dtc0 = Column(String(128), ForeignKey('dtcconfig.name'), nullable=True)
    dtc1 = Column(String(128), ForeignKey('dtcconfig.name'), nullable=True)
    dtc2 = Column(String(128), ForeignKey('dtcconfig.name'), nullable=True)
    dtc3 = Column(String(128), ForeignKey('dtcconfig.name'), nullable=True)


class DTCConfig(Base):
    # name format is "radarmode.name_dtc#"" e.g. WorldDay40.v01_dtc0
    __tablename__ = 'dtcconfig'
    name = Column(String(128), primary_key=True)
    codetype   = Column(String(128), nullable=False)  # either coherent, incoherent, uncoded
    # pulsename  = Column(String(128), nullable=False)  # either barker, m-code, alternating code, long pulse
    # code = Column(String(128), nullable=False)  # string of the hex of the code
    processing = Column(String(128), nullable=False)  # CohCode, IncohCodeFl, S, PLFFTS, etc.
    tx_frequency  = Column(Integer, nullable=False)
    rx_frequency  = Column(Integer, nullable=False)
    has_raw       = Column(Boolean, nullable=False)
    # has_cal       = Column(Boolean, nullable=False)
    # has_noise     = Column(Boolean, nullable=False)
    # baud_length   = Column(Integer, nullable=False) # length of each baud in microseconds
    # pulse_length  = Column(Integer, nullable=False) # length of the pulse in microseconds
    sample_length = Column(Integer, nullable=False) # length of the voltage sample in microseconds



