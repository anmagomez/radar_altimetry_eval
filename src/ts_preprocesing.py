import os
import time as tm
import numpy as np
from pylab import plt
from datetime import datetime as dt
import pandas as pd
import scipy.interpolate as sc

  
'''Functions developed by Sylvain Biancamaria'''
            
def decimalyear2yearmonthdayhour(self, decimaltime_curr):
    '''
    Function to convert time vector in decimal year to
    correponding year, month, day and hour for each date
    '''
    # Compute year for the input date in decimal year
    year_curr = int(np.floor(decimaltime_curr))
    # Compute number of days in current year for the input date
    #  (01 January = day 1)
    start_yearcurr = dt(year=year_curr, month=1, day=1)
    start_nextyear = dt(year=year_curr+1, month=1, day=1)
    totnumdays_yearcurr = (tm.mktime(start_nextyear.timetuple()) -
                           tm.mktime(start_yearcurr.timetuple()))/(24.0*3600.0)
    numdays_decimaltimecurr = int(np.floor((decimaltime_curr-year_curr) *
                                           totnumdays_yearcurr)) + 1
    # Compute month for the input date in decimal year
    ydm_curr = dt.strptime(str(year_curr)+" "+str(numdays_decimaltimecurr),
                           "%Y %j")
    month_curr = ydm_curr.month
    day_curr = ydm_curr.day
    # Compute hour for the input date in decimal year
    #  (numdays_decimaltimecurr-1 because 01 January
    #  19 hour is actualy day 0 from starting year + 19 hours)
    hour_curr = int((decimaltime_curr-year_curr -
                     float(numdays_decimaltimecurr-1)/totnumdays_yearcurr) *
                    totnumdays_yearcurr*24)
    return (year_curr, month_curr, day_curr, hour_curr)


def yearday2decimalyear(year_curr, day_curr):
    '''
    Function to convert date defined by its year (year_curr)
    and the number of days elapsed since the begining of the
    year (day_curr), with the covention that the first day is
    day 1 (and not 0) and the last day of the year is day 366
    for leap year and 365 for others.
    For example, February 3rd, 2019 corresponds to:
    year_curr = 2019
    day_curr = 34
    '''
    start_yearcurr = dt(year=int(year_curr), month=1, day=1)
    start_nextyear = dt(year=int(year_curr)+1, month=1, day=1)
    totnumdays_yearcurr = (tm.mktime(start_nextyear.timetuple()) -
                           tm.mktime(start_yearcurr.timetuple()))/(24.0*3600.0)
    decimalyear = year_curr + float(day_curr)/totnumdays_yearcurr
    return (decimalyear)


def yearmonthdayhourminutesec2decimalyear(year_curr, month_curr, day_curr,
                                          hour_curr, minute_curr, sec_curr):
    '''
    Function to convert time vector in year, month, day, hour,
    minute and second to correponding decimal year for each date
    '''
    date_curr = dt(int(year_curr), int(month_curr), int(day_curr),
                   int(hour_curr), int(minute_curr), int(sec_curr))
    dayyear_curr = int(date_curr.strftime("%j"))
    decimaltime_curr = yearday2decimalyear(year_curr, dayyear_curr - 1 +
                                           hour_curr/24.0 +
                                           minute_curr/(24.0*60.0) +
                                           sec_curr/(24.0*60.0*60.0))
    return decimaltime_curr

#####Other functions
def get_date_time_cols(df, date_fd):
    
    df['year'] =df[date_fd].dt.year
    df['month']=df[date_fd].dt.month
    df['day']  =df[date_fd].dt.day
    df['hour'] =12
    df['decimal_y'] = np.array(list(map(yearmonthdayhourminutesec2decimalyear,
                           df['year'].array,df['month'].array, df['day'].array, df['hour'].array, 
                             np.zeros(df['year'].shape), np.zeros(df['year'].shape))))
    return df

def convert_units(df,height_fd, origin='FEET', to='METER', check_col=True, unit_fd='unit',
                  gauge_fd='gauge_id', stations=None):
    '''Convert from ft to m and from m to f
    currently works for a dataframe that has a combination of ft, cm and m
    '''
    if origin=='FEET' and to=='METER':
        conversion_factor=0.3048
    if origin=='FEET' and to=='CM':
        conversion_factor=30.48
    if origin=='CM' and to=='METER':
        conversion_factor=0.01
    if origin=='METER' and to=='FEET':
        conversion_factor=3.28084
    if origin=='CM' and to=='FEET':
        conversion_factor=0.03281
    if origin=='METER' and to=='CM':
        conversion_factor=100
    look_for=origin
    if check_col:
        df_final=pd.DataFrame()
        if stations is None:
            stations=df[gauge_fd].unique()
            for lc in stations:
                df_local=df.loc[df['gauge_id']==lc].copy()
                units=df_local[unit_fd].iloc[0]
                if look_for==units:
                    df.loc[:, height_fd+'_'+to]=df[height_fd]*conversion_factor
                else:
                    df.loc[:, height_fd+'_'+to]=df[height_fd]
                df_final=pd.concat((df_final,df), axis=0)
        df_final=df_final.rename(columns={'height':'height_rw', 'height_'+to:'height'})
    
    return df