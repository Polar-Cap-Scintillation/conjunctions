# map_prn.py
# functions for converting GPS PRN numbers to NORAD satellite ID

# IMPORTANT: PRNs are reassigned as GPS satellites are comissioned/decomissioned
# This means a particular PRN may refer to different physical satellites
#   at different points in time.  NORAD IDs are associated with physical
#   satellites and do NOT change.  Time is required to properly map PRNs
#   to NORAD IDs.

from ftplib import FTP
import datetime as dt
import numpy as np

def retrieve_prn_mapping_info():
    # Retrieve PRN table from AGI FPI site and convert it to a dictionary
    # ftp://ftp.agi.com/pub/Catalog/Almanacs/SEM/GPSData.txt

    prn_table = []

    ftp = FTP('ftp.agi.com')
    ftp.login()
    ftp.retrlines('RETR pub/Catalog/Almanacs/SEM/GPSData.txt', callback=prn_table.append)
    ftp.quit()

    # parse lines of table
    header_length = 22      # number of header lines
    prn = np.empty((0,), dtype=int)
    svn = np.empty((0,), dtype=int)
    satnum = np.empty((0,), dtype=int)
    starttime = np.empty((0,), dtype=dt.date)
    endtime = np.empty((0,), dtype=dt.date)

    for line in prn_table[header_length:]:
        ls = line.split()
        # skip empty lines
        try:
            prn = np.append(prn, int(ls[0]))
        except IndexError:
            continue
        svn = np.append(svn, int(ls[1]))
        satnum = np.append(satnum, int(ls[2]))
        starttime = np.append(starttime, dt.date.fromisoformat(ls[7]))
        # handle current satelite (no endtime given)
        # this treatment assumes input will always be a time in the past
        try:
            endtime = np.append(endtime, dt.date.fromisoformat(ls[11]))
        except IndexError:
            endtime = np.append(endtime, dt.date.today())

    # organize into dictionary based on PRN
    prn_mapping_dict = {}
    for i in range(1,32):
        idx = np.argwhere(prn==i).flatten()
        prn_mapping_dict[i] = {'SVN':svn[idx], 'NORADID':np.array(satnum[idx]), 'STARTTIME':np.array(starttime[idx]), 'ENDTIME':np.array(endtime[idx])}

    return prn_mapping_dict


def find_date_index(startdates, enddates, targdate):
    # find the index where the target date is between the start and end dates
    # raises an error if the date does not correspond to a range in the PRN table

    try:
        tidx = np.argwhere((targdate>=startdates) & (targdate<enddates)).flatten()[0]
    except IndexError:
        raise ValueError('The specified PRN was not used on {}!'.format(targdate.isoformat()))

    return tidx


def prn2norad(prn, date):
    # map PRN to NORAD SAT ID

    mapping_dict = retrieve_prn_mapping_info()
    idx = find_date_index(mapping_dict[prn]['STARTTIME'], mapping_dict[prn]['ENDTIME'], date)

    return mapping_dict[prn]['NORADID'][idx]


def prn2svn(prn, date):
    # map PRN to SVN

    mapping_dict = retrieve_prn_mapping_info()
    idx = find_date_index(mapping_dict[prn]['STARTTIME'], mapping_dict[prn]['ENDTIME'], date)

    return mapping_dict[prn]['SVN'][idx]





if __name__=='__main__':

    print(prn2norad(18, dt.date(2019,3,21)))
    print(prn2svn(18, dt.date(2019,3,21)))
