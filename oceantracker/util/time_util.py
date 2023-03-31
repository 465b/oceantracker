from datetime import  datetime, timedelta
import dateutil.parser
import math

import numpy as np
# deal with date time operations,



def ot_time_zero(): return datetime(1970,1,1)

def seconds_to_date(s): return ot_time_zero() + timedelta(seconds=s) # better than total seconds which allows for computers time zone

def seconds_to_iso8601str(s): return seconds_to_date(s).isoformat() # better than total seconds which allows for computers time zone
def seconds_to_pretty_str(s, seconds= True):
    fmt="%Y-%m-%d %H:%M"
    if seconds: fmt +=":%S"
    s_str= seconds_to_date(s).strftime(fmt)
    return s_str
def seconds_to_short_date(s):
    fmt="%Y_%m_%d"
    s_str= seconds_to_date(s).strftime(fmt)
    return s_str
def iso8601str_to_seconds(s):  return date_to_seconds(date_from_iso8601str(s))

def date_to_seconds(date): return  (date-ot_time_zero()).total_seconds()

def diff(date1,date2) :  return  (date1 - date2).total_seconds()

def float_sec_to_time_delta(s): return np.asarray(s, dtype=np.float64).astype('timedelta64[s]')


def iso8601_str(d)   :  return d.isoformat()

def  day_hms(x):
    # convert decimal seconds to days:  hms string
    days = math.floor(x/24/3600)
    x -= days*24*3600
    hours   =  math.floor(x/3600)
    x -= hours*3600
    minutes =  math.floor(x/60)

    return "%02.0f %02.0f:%02.0f" % (days, hours, minutes)