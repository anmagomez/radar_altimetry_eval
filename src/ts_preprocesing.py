import os
import time as tm
import numpy as np
from pylab import plt
from datetime import datetime, timezone
import datetime as dt
import pandas as pd
import pytz
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
    start_yearcurr = datetime(year=year_curr, month=1, day=1)
    start_nextyear = datetime(year=year_curr+1, month=1, day=1)
    totnumdays_yearcurr = (tm.mktime(start_nextyear.timetuple()) -
                           tm.mktime(start_yearcurr.timetuple()))/(24.0*3600.0)
    numdays_decimaltimecurr = int(np.floor((decimaltime_curr-year_curr) *
                                           totnumdays_yearcurr)) + 1
    # Compute month for the input date in decimal year
    ydm_curr = datetime.strptime(str(year_curr)+" "+str(numdays_decimaltimecurr),
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
    start_yearcurr = datetime(year=int(year_curr), month=1, day=1)
    start_nextyear = datetime(year=int(year_curr)+1, month=1, day=1)
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
    date_curr = datetime(int(year_curr), int(month_curr), int(day_curr),
                   int(hour_curr), int(minute_curr), int(sec_curr))
    dayyear_curr = int(date_curr.strftime("%j"))
    decimaltime_curr = yearday2decimalyear(year_curr, dayyear_curr - 1 +
                                           hour_curr/24.0 +
                                           minute_curr/(24.0*60.0) +
                                           sec_curr/(24.0*60.0*60.0))
    return decimaltime_curr

def interp_ts12ts2_stat(dy_ts1, data_ts1, dy_ts2, data_ts2):
    '''
    Linearly interpolate data from a first time series (ts1) to date of a
    second time series (ts2). Then compute Pearson correlation coefficient
    between interpolated ts1 and ts2, Nash-Sutcliffe coefficient of ts2 anom
    wrt ts1 anomalies, and root mean square difference from ts2 anomalies wrt
    to ts1 anomalies.
    Inputs:
    - dy_ts1: decimal time for first time series (ts1)
    - data_ts1: data for first time series (ts1)
    - dy_ts2: decimal time for second time series (ts2)
    - data_ts2: data for second time series (ts2)
    Outputs:
    - datast1_2_ts2dy: ts1 data at ts2 dates
    - corr_ts1ts2: correlation coefficient between ts1 and ts2 data
    - ns_ts2: ts2 minus time mean Nash-Sutcliffe coefficient wrt ts1 minus
              time mean
    - rmsd_ts2: root mean square difference between ts2 minus time mean and
                ts1 minus time mean
    - ampl_ts1: difference between max and min for interpolated ts1 time
                series over common time period with ts2
    '''
    # Interpolate linearly in-situ measurements to altimetry measurement times
    datast1_2_ts2dy = sc.griddata(dy_ts1, data_ts1, dy_ts2, method='linear')
    # Compute the correlation between altimetry and discharge
    icommon = ((np.isnan(datast1_2_ts2dy) == 0) &
               (np.isnan(data_ts2) == 0)).nonzero()
    if len(icommon[0]) > 1:
        datats2_commonts1 = data_ts2[icommon]
        datats1_commonts2 = datast1_2_ts2dy[icommon]
        # Correlation coefficient
        vec2corrcoef = np.zeros((2, datats2_commonts1.size))
        vec2corrcoef[0, :] = datats2_commonts1
        vec2corrcoef[1, :] = datats1_commonts2
        matcorr_ts1ts2 = np.corrcoef(vec2corrcoef)
        corr_ts1ts2 = matcorr_ts1ts2[0, 1]
        # Nash-Sutcliffe coefficient
        diffts = (datats2_commonts1 - np.nanmean(datats2_commonts1)) -\
            (datats1_commonts2 - np.nanmean(datats1_commonts2))
        ns_ts2 = 1 - np.sum(np.square(diffts)) / \
            np.sum(np.square(datats1_commonts2 -
                             np.nanmean(datats1_commonts2)))
        # RMS from ts2 wrt ts1
        rmsd_ts2 = np.linalg.norm(diffts)/np.sqrt(diffts.size)
        # Amplitude of ts1 time series over common date with ts2
        ampl_ts1 = np.max(datats1_commonts2) - np.min(datats1_commonts2)
    else:
        corr_ts1ts2 = np.nan
        ns_ts2 = np.nan
        rmsd_ts2 = np.nan
        ampl_ts1 = np.nan
    return (datast1_2_ts2dy, corr_ts1ts2, ns_ts2, rmsd_ts2, ampl_ts1)


def comp_insi_alti(finsi, insiname, falti, altivsname, ncoldate, ncolh,
                   ncolgeoid, outdir, fstat, fvolo=None, volovsname=None,
                   fstatvol=None):
    '''
    Compare in situ discharge and altimetry time series
    '''
    # Load in situ water level time series
    (indy, inh) = load_insi_volformat(finsi)
    # Load Altis altimetry time series
    (altiyear, altimonth, altiday, altihour, altiminute,
     altiwelev) = load_altis(falti, ncoldate, ncolh, ncolgeoid)
    # Interpolate in situ time series to altis dates
    altidy = np.array(list(map(yearmonthdayhourminutesec2decimalyear,
                               altiyear, altimonth, altiday, altihour,
                               altiminute, np.zeros(altiday.shape))))
    
    ##AMG Here it starts
    
    (inh_2_altidy, corr_alti_insi, ns_alti_insi, rmse_alti_insi,
     ampl_insi) = interp_ts12ts2_stat(indy, inh, altidy, altiwelev)
    # Save statistics between alti and insi time series
    if (fstat is not None) and (os.path.exists(fstat)):
        line2save = (altivsname+';'+insiname+';' +
                     str(np.round(corr_alti_insi, decimals=2))+';' +
                     str(np.round(ns_alti_insi, decimals=2))+';' +
                     str(np.round(rmse_alti_insi, decimals=2))+';' +
                     str(np.round(ampl_insi, decimals=2))+'\n')
        ffstat = open(fstat, 'a')
        ffstat.write(line2save)
        ffstat.close()
    # Plots
    plt.close('all')
    plt.ion()
    # Plot in situ discharge and altimetry water elevation time series
    plt.figure(1)
    ax = plt.gca()
    plt.plot(indy, inh-np.nanmean(inh_2_altidy), 'bo', markersize=1)
    plt.plot(altidy, altiwelev-np.nanmean(altiwelev), 'ro', markersize=2.5)
    plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    plt.grid(b='on', axis='x')
    plt.legend((insiname, altivsname), loc='best', numpoints=1,
               framealpha=0.6)
    # plt.xlim(np.floor(np.min(altidy)), np.ceil(np.max(altidy)))
    plt.xlim(2016.0, np.ceil(np.max(altidy)))
    plt.xlabel('Time (decimal year)', size=15, weight='bold')
    plt.ylabel('Water elev. anomalies (m)', size=15, weight='bold')
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(15)
        tick.label.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(15)
        tick.label.set_fontweight('bold')
    plt.title('alti/insi corr. = '+str(np.round(corr_alti_insi, 2)),
              size=15, weight='bold')
    plt.savefig(os.path.join(outdir, 'comp_'+altivsname+'_'+insiname+'.png'),
                bbox_inches='tight')
    plt.figure(2)
    ax = plt.gca()
    plt.plot(indy, inh, 'bo', markersize=1)
    plt.plot(altidy, altiwelev, 'ro', markersize=2.5)
    plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    plt.grid(b='on', axis='x')
    plt.legend((insiname, altivsname), loc='best', numpoints=1,
               framealpha=0.6)
    # plt.xlim(np.floor(np.min(altidy)), np.ceil(np.max(altidy)))
    plt.xlim(2016.0, np.ceil(np.max(altidy)))
    plt.xlabel('Time (decimal year)', size=15, weight='bold')
    plt.ylabel('Water elev. wrt geoid (m)', size=15, weight='bold')
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(15)
        tick.label.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(15)
        tick.label.set_fontweight('bold')
    plt.title('alti/insi corr. = '+str(np.round(corr_alti_insi, 2)), size=15,
              weight='bold')
    plt.savefig(os.path.join(outdir, 'abs_comp_'+altivsname+'_'+insiname +
                             '.png'), bbox_inches='tight')


def comp_multiple_altis_insi():
    '''
    Compare altimetry water elevation time series from GDR file (generated
    with Altis) to in situ water elevation time series
    '''
    # Input directory
    altidir = '/home/biancama/partagewindows/ao_projets/LOCSS/alti/ts/'
    insidir = ('/home/biancama/partagewindows/ao_projets/LOCSS/insitu/ts_'
               'insitu/timeseries/')
    # Output directory
    outdir = ('/home/biancama/partagewindows/ao_projets/LOCSS/comp_alti_insi'
              '/altis_only')
    # Altis time series information
    ncoldate = 'date'
    ncolh = 'ice1_ku_SurfHeight_alti_median'
    ncolgeoids3 = 'geoid_01_median'
    ncolgeoidj3 = 'geoid_eigen6c4d_median'
    ncolgeoid_default = ncolgeoids3
    # Load file with in situ gage name associated with corresponding VS name
    fngages = ('/home/biancama/partagewindows/ao_projets/LOCSS/altis_'
               'locsslakes.csv')
    (altilist, insilist) = np.loadtxt(fngages, delimiter=';', unpack=True,
                                      dtype='U40,U75')
    # Create the text file to store alti/insi statistics
    fstat = os.path.join(outdir, 'stat_alti_insi.csv')
    ffstat = open(fstat, 'w')
    ffstat.write(('altivs_name;insigage_name;Corr;NSE_anom;RMSE_anom_m;'
                  'amplitudeinsi_m\n'))
    ffstat.close()
    # Loop over all VS/gage couples
    for icoupl in range(len(altilist)):
        altivsname = altilist[icoupl]
        falti = os.path.join(altidir, altivsname+'.csv')
        if altivsname.find('Sentinel-3A') > -1:
            ncolgeoid = ncolgeoids3
        elif altivsname.find('Jason-3') > -1:
            ncolgeoid = ncolgeoidj3
        else:
            ncolgeoid = ncolgeoid_default
        insiname = insilist[icoupl]
        finsi = os.path.join(insidir, insiname+'.txt')
        print('Comparison '+altivsname+' vs '+insiname)
        comp_insi_alti(finsi, insiname, falti, altivsname, ncoldate, ncolh,
                       ncolgeoid, outdir, fstat)

def load_altis(falti, ncoldate, ncolh, ncolgeoid=None, nodataalti=-9999, wse_ref='g'):
    ##TODO: Pending to raise exceptions in errors
    '''
    Function to load altimetry water elevation in AlTiS csv format
    Inputs:
    - falti: AlTiS csv file
    - ncoldate: name of the column with date of each sample time series
    - ncolh: name of the column with water elevation time series
    - nodataalti: no data values in the water elevation time series
    - ncolgeoid: name of the column with geoid elevation time series. Default is None
    - wse_ref: type of water surface elevation to extract. 
               'e' means reference to the ellipsoid.
               'g' (Default) means reference to the geoid
               If wse_ref='e', ncolgeoid has to contain a valid column
    Outputs:
    - alti_year: year of all samples in the time series
    - alti_month: month of all samples in the time series
    - alti_day: day of all samples in the time series
    - alti_hour: hour of all samples in the time series
    - alti_minute: minute of all samples in the time series
    - alti_height: water elevation of all samples in the time series
                   referenced to the mission reference geoid (and not the
                   ellipsoid)
    '''
    # Retrieve in-situ water level file header
    fin = open(falti)
    isheader = fin.readline()
    fin.close()
    isheader_split = (isheader.replace('\n', '')).split(',')
    # Extract water elevation
    patterncol = ncolh
    icolh = [i for i, s in enumerate(isheader_split)
             if patterncol.lower() in s.lower()]
    if len(icolh) == 0:
        print(('Error: no column '+ncolh+' in '+falti))
        alti_height = None
    elif len(icolh) > 1:
        print(('Error: more than one column '+ncolh+' in '+falti))
        alti_height = None
    elif len(icolh) == 1:
        alti_height = np.loadtxt(falti, skiprows=1, delimiter=',',
                                 usecols=[icolh[0]])
    # Get lines with no data values
    ivaliddata = (alti_height > nodataalti).nonzero()
    if wse_ref=='e':
        if ncolgeoid is None:
            print(('Geoid column not specified. WSE will be referenced to the geoide. ncolgeoid= '+str(ncolgeoid)))
            alti_geoid=None
        else:
            # Extract geoid
            patterncol = ncolgeoid
            icolg = [i for i, s in enumerate(isheader_split)
                     if patterncol.lower() in s.lower()]
            if len(icolg) == 0:
                print(('Error: no column '+ncolgeoid+' in '+falti))
                alti_geoid = None
            elif len(icolg) > 1:
                print(('Error: more than one column '+ncolgeoid+' in '+falti))
                alti_geoid = None
            elif len(icolg) == 1:
                alti_geoid = np.loadtxt(falti, skiprows=1, delimiter=',',
                                        usecols=[icolg[0]])
    else: 
        alti_geoid=None
    
    # Extract date
    patterncol = ncoldate
    icold = [i for i, s in enumerate(isheader_split)
             if patterncol.lower() in s.lower()]
    if len(icold) == 0:
        print(('Error: no column '+ncoldate+' in '+falti))
        alti_date = None
    elif len(icold) > 1:
        print(('Error: more than one column '+ncoldate+' in '+falti))
        alti_date = None
    elif len(icold) == 1:
        alti_date = np.loadtxt(falti, skiprows=1, delimiter=',',
                               usecols=[icold[0]], dtype='U')
    if alti_date is None:
        alti_year = None
        alti_month = None
        alti_day = None
        alti_hour = None
        alti_minute = None
    else:
        split_vecdate = ' '.join(alti_date.tolist()).replace(':', ' ')\
                        .replace('-', ' ').split(' ')
        alti_year = np.array(list(map(int, split_vecdate[::6])))
        alti_month = np.array(list(map(int, split_vecdate[1::6])))
        alti_day = np.array(list(map(int, split_vecdate[2::6])))
        alti_hour = np.array(list(map(int, split_vecdate[3::6])))
        alti_minute = np.array(list(map(int, split_vecdate[4::6])))
        if alti_geoid is not None:
            alti_wse = alti_height+alti_geoid
        else:
            alti_wse = alti_height
    return (alti_year[ivaliddata], alti_month[ivaliddata],
            alti_day[ivaliddata], alti_hour[ivaliddata],
            alti_minute[ivaliddata],
            alti_wse[ivaliddata])



#####Other functions

def open_match_station_altis(g_path,altis_name, st_fd, df_gts, st_id,altis_date_fd, altis_height_fd, nodataalti=-9999, ncolgeoid=None,wse_ref='g'):
    '''Match data from altimeter with data from the station
    Inputs:
        g_path: Altis observation directory
        altis_name: Name of the file without extension
        df_gts: Dataframe with ground observations
        st_fd: station column name in df_gts
        st_id; station id
        altis_date_fd: altis date column name
        altis_height_fd: altis height column name
        nodataalti: No data value in altis file. Default:-9999
        ncolgeoid: Column with the geoid model values to get WSE relative to the ellipsoide. 
                   Default None: WSE is calculated relative to the geoid
        wse_ref: Reference to get water surface elevation (WSE). Default: 'g', relative to the geoid
                 'e', relative to the ellipsoide
    Outputs:
        1) altis dataframe
        2) ground observations for the station st_id
    '''
    #Read altis
    (altiyear, altimonth, altiday, altihour, altiminute, altiwelev)=load_altis(g_path+altis_name+'.csv', 
                                                                               altis_date_fd, altis_height_fd, 
                                                                               nodataalti=-9999, ncolgeoid=ncolgeoid, wse_ref=wse_ref)
    # Compute decimal year from year/month/day from altimetry time series
    altidy = np.array(list(map(yearmonthdayhourminutesec2decimalyear,
                           altiyear, altimonth, altiday, altihour,
                           altiminute, np.zeros(altiday.shape))))
    
    dic_altis={'decimal_y':altidy,'height':altiwelev,'year':altiyear, 'month':altimonth, 'day':altiday, 'hour':altihour}
    
    df_altis=pd.DataFrame(dic_altis)
    
    # df_altis['date']=pd.datetime(df_altis['year'].str.cat(df_altis['month'], sep='-').
    #                              cat(df_altis['day'], sep='-').cat(df_altis['hour'], sep=' ').cat(':00',sep=''), utc=True)
    df_altis[altis_date_fd]= df_altis['year'].apply(str)+'-'+df_altis['month'].apply(str)+'-'+df_altis['day'].apply(str)+' '+df_altis['hour'].apply(str)+':00'
    df_altis[altis_date_fd]=pd.to_datetime(df_altis['date'], utc=True)

    df_altis['name_altis']=altis_name
    df_gts_st=df_gts.loc[df_gts[st_fd]==st_id].copy()
    
    return (df_altis, df_gts_st)

def get_common_period(df_ts1, df_ts2, date_ts1_fd, date_ts2_fd, delta_days=False, ndays=0):
    '''Get the common period of time between two time series, df_ts1 and df_ts2
        WARNING: If not timezone information in any of the dataframe, utc is assumed 
        It can return the common period + or - a number of days based on df_ts1
        Assume both dataframe datetime is utc
        If delta_days=False, ndays is assume =0, returns a common dataframe with the interception between the two dataframes
        If delta_days=True, ndays has to be different from 0
            Rules:
                if df_ts1 starts earlier than df_ts2, common period starts at initial time of df_ts2
                if df_ts1 starts later than df_ts2, common period starts ndays earlier than initial time of df_ts1
                if df_ts1 ends earlier than df_ts2, common period ends ndays later than final time of df_ts1
                if df_ts1 ends later than df_ts2, common period end at the final time of df_ts2
    '''
    
    if delta_days==True and ndays<=0:
        #TODO:Convert this to raise exceptions
        print('Error ndays cannot be 0 or lower if delta_days=True')
        return None
    if delta_days==False:
        ndays=0
    
    utc=pytz.utc
    
    #Min and max dates in df_ts1 and df_ts2
    if df_ts1[date_ts1_fd].dt.tz is None:
        df_ts1[date_ts1_fd]=[utc.localize(date) for date in df_ts1[date_ts1_fd]]
        
    if df_ts2[date_ts2_fd].dt.tz is None:
        df_ts2[date_ts2_fd]=[utc.localize(date) for date in df_ts2[date_ts2_fd]]
        
    earlier_date_ts1=min(df_ts1[date_ts1_fd])
    final_date_ts1=max(df_ts1[date_ts1_fd])
        
    earlier_date_ts2=min(df_ts2[date_ts2_fd])
    final_date_ts2=max(df_ts2[date_ts2_fd])

    if earlier_date_ts1 >= earlier_date_ts2:
        initial_date=earlier_date_ts1 - dt.timedelta(days=ndays)
    else:
        initial_date=earlier_date_ts2

    if final_date_ts1 >= final_date_ts2:
        end_date=final_date_ts2
    else:
        end_date=final_date_ts1 + dt.timedelta(days=ndays)
    
    #filter dataframes
    df_ts1=df_ts1.loc[(df_ts1[date_ts1_fd]>=initial_date)&(df_ts1[date_ts1_fd]<=end_date)].copy()
    df_ts2=df_ts2.loc[(df_ts2[date_ts2_fd]>=initial_date)&(df_ts2[date_ts2_fd]<=end_date)].copy()
    
    return (df_ts1, df_ts2)

def moving_window_around_date(df, date, delta, v_fd, d_fd):
    '''Moving window of a value around a date +- n days defined by delta
        Inputs:
            df: Dataframe containing the dates and values
            date: date around the one the moving window will be done
            v_fd: name of the value column in df
            d_fd: name of the date column in df
        Output:
           In order in the window
           median
           mean 
           std: standard deviation 
           number of not null values used for the mean and the median
           
    '''
    df_t=df.loc[(df[d_fd]>=(date-dt.timedelta(days=delta)))&(df[d_fd]<=(date+dt.timedelta(days=delta)))]
    return df_t[v_fd].median(skipna=True),df_t[v_fd].mean(skipna=True),df_t[v_fd].std(skipna=True), df_t[v_fd].count()

def get_comp_metrics(ts_obs,ts_est):
    
    icommon = ((np.isnan(ts_est) == 0) &
               (np.isnan(ts_obs) == 0)).nonzero()
    
    if len(icommon[0]) > 1:
        datats2_commonts1 = ts_obs[icommon]
        datats1_commonts2 = ts_est[icommon]
        # Correlation coefficient
        vec2corrcoef = np.zeros((2, datats2_commonts1.size))
        vec2corrcoef[0, :] = datats2_commonts1
        vec2corrcoef[1, :] = datats1_commonts2
        matcorr_ts1ts2 = np.corrcoef(vec2corrcoef)
        corr_ts1ts2 = matcorr_ts1ts2[0, 1]
        # Nash-Sutcliffe coefficient
        diffts = (datats2_commonts1 - np.nanmean(datats2_commonts1)) -\
            (datats1_commonts2 - np.nanmean(datats1_commonts2))
        ns_ts2 = 1 - np.sum(np.square(diffts)) / \
            np.sum(np.square(datats1_commonts2 -
                             np.nanmean(datats1_commonts2)))
        # RMS from ts2 wrt ts1
        rmsd_ts2 = np.linalg.norm(diffts)/np.sqrt(diffts.size)
        # Amplitude of ts1 time series over common date with ts2
        ampl_ts1 = np.max(datats1_commonts2) - np.min(datats1_commonts2)
    else:
        corr_ts1ts2 = np.nan
        ns_ts2 = np.nan
        rmsd_ts2 = np.nan
        ampl_ts1 = np.nan
    return (corr_ts1ts2, ns_ts2, rmsd_ts2, ampl_ts1)

def get_date_time_cols(df, date_fd, has_hour=False, has_min=False):
    ##TODO: Calculate minutes
    df['year'] =df[date_fd].dt.year
    df['month']=df[date_fd].dt.month
    df['day']  =df[date_fd].dt.day
    if has_hour==False:
        df['hour'] =12
    else: 
        df['hour']=df[date_fd].dt.hour
    
    df['decimal_y'] = np.array(list(map(yearmonthdayhourminutesec2decimalyear,
                           df['year'].array,df['month'].array, df['day'].array, df['hour'].array, 
                             np.zeros(df['year'].shape), np.zeros(df['year'].shape))))
    return df

def convert_units(df,height_fd, origin='FEET', to='METER', check_col=True, unit_fd='unit',
                  gauge_fd='gauge_id', stations=None):
    '''Convert from ft to m and from m to f
    currently works for a dataframe that has a combination of ft, cm and m and just 2 units. 
    inputs:
        df: Dataframe with the observations to convert. If is has several units, unit_fd has to indicate the column with unit information
        origin: Unit origin of the conversion
        to: unit to convert
        stations: group of stations to convert. If None (default) if will convert all values in df with the origin units
        height_fd: Column with the values to convert
        gauge_fd: name of the column of the station. If None the conversion is apply to all the height_fd column
        check_col: Verify which rows need to be converted. Default True. False, will convert the whole dataframe df
    output:
        dataframe with the units converted. Old values are saved in height_rw column, new values are storage in column height 
        
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
                # print('\n in function\n')
                # print(lc, df.loc[df[gauge_fd]==lc].shape)
                df_local=df.loc[df[gauge_fd]==lc].copy()
                units=df_local[unit_fd].iloc[0]
                if look_for==units:
                    df_local.loc[:, height_fd+'_'+to]=df_local[height_fd]*conversion_factor
                else:
                    df_local.loc[:, height_fd+'_'+to]=df_local[height_fd]
                df_final=pd.concat((df_final,df_local), axis=0)
                # print(df_final.shape, lc+' '+units)
    else:
        if stations is None:
            df_final=df.copy()
            df_final.loc[:, height_fd+'_'+to]=df_final[height_fd]*conversion_factor 
        else:
            df_final=pd.DataFrame()
            for lc in stations: #TODO: Validate not all the stations in stations can be in df
                df_local=df.loc[df_local[gauge_fd]==lc].copy()
                if not df_local.empty:
                    df_local.loc[:, height_fd+'_'+to]=df_local[height_fd]*conversion_factor 
                    df_final=pd.concat((df_final,df_local), axis=0)
    df_final=df_final.rename(columns={height_fd:'height_rw', height_fd+'_'+to:'height'})
    
    return df_final