import os
import time as tm
import numpy as np
from pylab import plt
from datetime import datetime, timezone
import datetime as dt
import pandas as pd
import pytz
import scipy.interpolate as sc
import scipy.stats
# from icecream import ic

  
'''Functions developed by Sylvain Biancamaria'''
            
def decimalyear2yearmonthdayhour(decimaltime_curr):
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

def load_altis(falti, ncoldate, ncolh, ncolobs,ncolgeoid=None, nodataalti=-9999, wse_ref='g'):
    ##TODO: Pending to raise exceptions in errors
    '''
    Function to load altimetry water elevation in AlTiS csv format
    Inputs:
    - falti: AlTiS csv file
    - ncoldate: name of the column with date of each sample time series
    - ncolh: name of the column with water elevation time series
    - ncolobs: name of the column with the number of samples 
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
    - alti_obsample: samples for each median data observation from altimeter
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
        
    # Extract number of samples
    patterncol = ncolobs
    icolo = [i for i, s in enumerate(isheader_split)
             if patterncol.lower() in s.lower()]
    if len(icolo) == 0:
        print(('Error: no column '+ncolobs+' in '+falti))
        alti_obsample = None
    elif len(icolo) > 1:
        print(('Error: more than one column '+ncolobs+' in '+falti))
        alti_obsample = None
    elif len(icolo) == 1:
        alti_obsample = np.loadtxt(falti, skiprows=1, delimiter=',',
                                 usecols=[icolo[0]])
    
    
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
            alti_wse[ivaliddata], 
            alti_obsample[ivaliddata])



#####Other functions made by Angelica

def convert_tuple_2_date(date_tuple):
    ''' Recieves a tuple with year, month, day and hour and convert it to a datetime '''
    #TODO Validate function
    y=date_tuple[0]
    m=date_tuple[1]
    d=date_tuple[2]
    h=date_tuple[3]
    return dt.datetime(y,m,d,h)

def open_match_station_altis(g_path,altis_name, st_fd, df_gts, st_id,altis_date_fd, altis_height_fd, n_obs_altis,nodataalti=-9999, ncolgeoid=None,wse_ref='g'):
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
    (altiyear, altimonth, altiday, altihour, altiminute, altiwelev, altiobs)=load_altis(g_path+altis_name+'.csv', 
                                                                               altis_date_fd, altis_height_fd, n_obs_altis,
                                                                               nodataalti=-9999, ncolgeoid=ncolgeoid, wse_ref=wse_ref)
    # Compute decimal year from year/month/day from altimetry time series
    altidy = np.array(list(map(yearmonthdayhourminutesec2decimalyear,
                           altiyear, altimonth, altiday, altihour,
                           altiminute, np.zeros(altiday.shape))))
    
    dic_altis={'decimal_y':altidy,'height':altiwelev,'nobs_altis':altiobs,'year':altiyear, 'month':altimonth, 'day':altiday, 'hour':altihour}
    
    df_altis=pd.DataFrame(dic_altis)
    
    # df_altis['date']=pd.datetime(df_altis['year'].str.cat(df_altis['month'], sep='-').
    #                              cat(df_altis['day'], sep='-').cat(df_altis['hour'], sep=' ').cat(':00',sep=''), utc=True)
    df_altis[altis_date_fd]= df_altis['year'].apply(str)+'-'+df_altis['month'].apply(str)+'-'+df_altis['day'].apply(str)+' '+df_altis['hour'].apply(str)+':00'
    df_altis[altis_date_fd]=pd.to_datetime(df_altis['date'], utc=True)

    df_altis['name_altis']=altis_name
    df_gts_st=df_gts.loc[df_gts[st_fd]==st_id].copy()
    
    return (df_altis, df_gts_st)

def get_common_period(df_ts1, df_ts2, date_ts1_fd, date_ts2_fd, delta_days=False, ndays=0, decimal_date=False):
    '''Get the common period of time between two time series, df_ts1 and df_ts2
        WARNING: If working with dates and not timezone information in any of the dataframe, utc is assumed 
        It can return the common period + or - a number of days based on df_ts1
        Assume both dataframe datetime is utc or a decimal time with the same reference initial time
        If decimal_date=False (default) datetime is assumed delta_days and ndays are valid to work with. 
        If decimal_date is true, common period is +- .1 decimal units (#TODO: make this more elegant, check if decimal columns are decimal)
        If delta_days=False, ndays is assume =0, returns a common dataframe with the interception between the two dataframes
        If delta_days=True, ndays has to be different from 0
            Rules:
                if df_ts1 starts earlier than df_ts2, common period starts at initial time of df_ts2
                if df_ts1 starts later than df_ts2, common period starts ndays earlier than initial time of df_ts1
                if df_ts1 ends earlier than df_ts2, common period ends ndays later than final time of df_ts1
                if df_ts1 ends later than df_ts2, common period end at the final time of df_ts2
    '''
    if decimal_date:
        #TODO:Convert this to raise exceptions
        if delta_days:
            print('Error ndays or delta_days will be ignored if decimal_date=True')
            return None
        
        earlier_date_ts1=min(df_ts1[date_ts1_fd])
        final_date_ts1=max(df_ts1[date_ts1_fd])
        
        earlier_date_ts2=min(df_ts2[date_ts2_fd])
        final_date_ts2=max(df_ts2[date_ts2_fd])
        
        if earlier_date_ts1 >= earlier_date_ts2:
            initial_date=earlier_date_ts1 - 0.1
        else:
            initial_date=earlier_date_ts2

        if final_date_ts1 >= final_date_ts2:
            end_date=final_date_ts2
        else:
            end_date=final_date_ts1 + 0.1
        
        #filter dataframes
        df_ts1=df_ts1.loc[(df_ts1[date_ts1_fd]>=initial_date)&(df_ts1[date_ts1_fd]<=end_date)].copy()
        df_ts2=df_ts2.loc[(df_ts2[date_ts2_fd]>=initial_date)&(df_ts2[date_ts2_fd]<=end_date)].copy()
        
    else:
    
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

def moving_window_around_date(df, date, delta, v_fd, d_fd, method='median'):
    '''Moving window of a value around a date +- n days defined by delta
        Inputs:
            df: Dataframe containing the dates and values
            date: date around the one the moving window will be done
            v_fd: name of the value column in df
            d_fd: name of the date column in df
            method: Method used to get the value within the window. Default median, Other option 'linear' for linear interpolation. 
                    if option chosen is 'linear' and only one value is on the window
        Output:
           In order in the window
           median or interpolate value depending on method
           mad in the window
           mean in the window
           std: standard deviation 
           number of not null values used for the mean and the median
           criteria: method used 'median' or 'linear'
           
    '''
    df_t=df.loc[(df[d_fd]>=(date-dt.timedelta(days=delta)))&(df[d_fd]<=(date+dt.timedelta(days=delta)))].copy()
    n_cd=df_t[v_fd].count()
    if n_cd==1 and method=='linear':
        value=np.nan
        mad=np.nan
        mean_w=np.nan
        criteria=np.nan
    elif n_cd>=1:
        mean_w=df_t[v_fd].mean(skipna=True)
        mad=np.median(abs(df_t[v_fd]-df_t[v_fd].median(skipna=True)))
        if method=='median':
            value=df_t[v_fd].median(skipna=True)
            criteria='median'
        elif method=='linear':
            df_t['decimal_y']= np.array(list(map(yearmonthdayhourminutesec2decimalyear,
                                                      df_t[d_fd].array.year,
                                                      df_t[d_fd].array.month, 
                                                      df_t[d_fd].array.day, 
                                                      df_t[d_fd].array.hour,
                                                      df_t[d_fd].array.minute, 
                                                      df_t[d_fd].array.second)))
            date_decimal=yearmonthdayhourminutesec2decimalyear(date.year,
                                                      date.month, 
                                                      date.day, 
                                                      date.hour,
                                                      date.minute, date.second)
            # print(df_closer[d_fd].to_numpy(), date)
            # print(df_closer['decimal_y'].to_numpy(), date_decimal)
            # print(df_t[v_fd], 'Values\n')
            value = sc.griddata(df_t['decimal_y'].to_numpy(), df_t[v_fd].to_numpy(), date_decimal, method=method)
            # print(value, 'This is the value')
            criteria='linear'
    else:
        value=np.nan
        mad=np.nan
        mean_w=np.nan
        criteria=np.nan
    return value,mad, mean_w,df_t[v_fd].std(skipna=True), n_cd, criteria

def closer_value_around_date(df, date, delta, v_fd, d_fd, method='median'):
    '''Moving window of a value around a date (targeting date) +- n days defined by delta
        Inputs:
            df: Dataframe containing the dates and values
            date: date around the one the closer value will be estimated
            v_fd: name of the value column in df
            d_fd: name of the date column in df
            method: method implemented when more than one value is closer to the date. Default median. Options: 'linear' for linear interpolation
                    scipy.interpolate.griddata is used to do the linear interpolation
        Output:
            val_cd: value obtained from the analysis 
            type_cd: 'closed', if only one date was the closest, 'median': if more than two dates where equally closed
            ndays_cd: amount of days between the targe date and the closest date
            n_cd: number of data used to calculate value
            mea_cd: mean of the values in case ndays >1. nan if ndays==1
            std_cd: std of the values in case ndays>1. nan if ndays==1
            no_days: distance between the target date and the closest day
    '''
    df_t=df.loc[(df[d_fd]>=(date-dt.timedelta(days=delta)))&(df[d_fd]<=(date+dt.timedelta(days=delta)))].copy()
    if not df_t.empty:
        df_t=df_t.sort_values(by=[d_fd])
        df_t['diff_days']=df_t[d_fd].apply(lambda x: (date-x)/np.timedelta64(1,'D')).abs()
        closer_date=df_t['diff_days'].min()
        df_closer=df_t.loc[df_t['diff_days']==closer_date].copy()
        no_days=closer_date
        n_cd = df_closer[v_fd].count()
    else:
        
        n_cd=0
        
    if n_cd==1:
        # print('Solo hay uno que esta cerca')
        value = df_closer[v_fd].iloc[0]
        criteria='closer'
        mean=np.nan
        std=np.nan
        mad=np.nan
        
    elif n_cd==0:
        value=np.nan
        criteria=np.nan
        mean=np.nan
        std=np.nan
        mad=np.nan
        closer_date=-1
        no_days=-1
        
    else:
        # df_closer=df_closer.sort_values(d_fd).copy()
        if method=='linear':
            # Interpolate linearly in-situ measurements to altimetry measurement times
            # datast1_2_ts2dy = sc.griddata(dy_ts1, data_ts1, dy_ts2, method='linear')
            df_closer['decimal_y']= np.array(list(map(yearmonthdayhourminutesec2decimalyear,
                                                      df_closer[d_fd].array.year,
                                                      df_closer[d_fd].array.month, 
                                                      df_closer[d_fd].array.day, 
                                                      df_closer[d_fd].array.hour,
                                                      df_closer[d_fd].array.minute, df_closer[d_fd].array.second)))
            date_decimal=yearmonthdayhourminutesec2decimalyear(date.year,
                                                      date.month, 
                                                      date.day, 
                                                      date.hour,
                                                      date.minute, date.second)
            # print(df_closer[d_fd].to_numpy(), date)
            # print(df_closer['decimal_y'].to_numpy(), date_decimal)
            # print(df_closer[v_fd], 'Values')
            value = sc.griddata(df_closer['decimal_y'].to_numpy(), df_closer[v_fd].to_numpy(), date_decimal, method=method)
            # print(value, 'This is the value')
            criteria='linear'
            mean=np.nan
            
        if method=='median':    
            value =df_closer[v_fd].median(skipna=True)
            criteria='median'
        mean=df_closer[v_fd].mean(skipna=True)
        std=df_closer[v_fd].std(skipna=True)
        mad=np.median(abs(df_closer[v_fd]-value))
            
        
    # dict_output={'val_cd_d_'+str(delta)+:value,
    #              'type_cd_d_'+str(delta):criteria,
    #              'ndays_cd_d_'+str(delta): closer_date,
    #              'n_cd_d_'+str(delta): df_closer[v_fd].count(), 
    #              'mea_cd_d_'+str(delta):mean, 
    #              'std_cd_d_'+str(delta): std
    #             }
    
    val_cd = value
    type_cd = criteria
    ndays_cd= closer_date
    mea_cd=mean
    std_cd=std
    mad_cd=mad
    # val_cd,type_cd,ndays_cd,mea_cd, std_cd,n_cd
    return val_cd,type_cd,ndays_cd,mea_cd, mad_cd, std_cd,n_cd, no_days

# Pending add interpolation 

#https://stackoverflow.com/questions/32237862/find-the-closest-date-to-a-given-date
def get_comp_metrics(ts_obs,ts_est):
    ''' Compare altis and insitu data
        
        Outputs (In order):
        results: dictionary with the keys
            PR: pearson correlation
            PR_p_val: pearson correlation p-value
            RHO: spearman correlation
            RHO_p_val: spearman correlation p-value
            ns_ts2: NSC
            rmsd_ts2: RMSE
            ampl_ts1: Amplitude of ts1 time series over common date with ts2
            me: Mean Error (#TODO Pending check)
            ve: Variance of the error
            datats2_commonts1.size: Size of the observations
            datats1_commonts2.size: Size of Altis GDR
            datats2_commonts1: observations
            datats1_commonts2: GDR
    '''
    icommon = ((np.isnan(ts_est) == 0) &
               (np.isnan(ts_obs) == 0)).nonzero()
    
    if len(icommon[0]) > 1:
        datats2_commonts1 = ts_obs[icommon]
        datats1_commonts2 = ts_est[icommon]
        # Correlation coefficient
        vec2corrcoef = np.zeros((2, datats2_commonts1.size))
        vec2corrcoef[0, :] = datats2_commonts1 - np.nanmean(datats2_commonts1)
        vec2corrcoef[1, :] = datats1_commonts2 - np.nanmean(datats1_commonts2)
        matcorr_ts1ts2 = np.corrcoef(vec2corrcoef)
        corr_ts1ts2 = matcorr_ts1ts2[0, 1]
        #Adding other values
        p_corr_scipy,p_value_scipy = scipy.stats.pearsonr(vec2corrcoef[0, :], vec2corrcoef[1, :])
        s_corr_scipy,s_value_scipy = scipy.stats.spearmanr(vec2corrcoef[0, :], vec2corrcoef[1, :])
        #error at each time
        ei=vec2corrcoef[0, :]-vec2corrcoef[1, :]
        ###
        # Nash-Sutcliffe coefficient
        # diffts = (datats2_commonts1 - np.nanmean(datats2_commonts1)) - (datats1_commonts2 - np.nanmean(datats1_commonts2))
        diffts = (vec2corrcoef[0, :] - np.nanmean(vec2corrcoef[0, :])) - (vec2corrcoef[1, :] - np.nanmean(vec2corrcoef[1, :]))
        ns_ts2 = 1 - np.sum(np.square(diffts)) / \
            np.sum(np.square(vec2corrcoef[0, :] -
                             np.nanmean(vec2corrcoef[1, :])))
        # RMS from ts2 wrt ts1
        # rmsd_ts2 = np.linalg.norm(diffts)/np.sqrt(diffts.size)
        rmsd_ts2 = np.sqrt(np.nanmean(np.square(ei)))
        # Amplitude of ts1 time series over common date with ts2
        ampl_ts1 = np.max(datats1_commonts2) - np.min(datats1_commonts2)
        # Mean Error non bias removed
        me=np.nanmean((datats2_commonts1-datats1_commonts2))
        # Mean error bias removed
        me_br=np.nanmean((vec2corrcoef[0, :]-vec2corrcoef[1, :]))
        #variance of error
        ve=np.nanmean(np.square(ei - me))
        size_obs=datats2_commonts1.size
        size_est=datats1_commonts2.size
        
    else:
        datats2_commonts1=np.nan
        datats1_commonts2=np.nan
        corr_ts1ts2 = np.nan
        ns_ts2 = np.nan
        rmsd_ts2 = np.nan
        ampl_ts1 = np.nan
        me=np.nan
        ve=np.nan
        p_corr_scipy=np.nan
        p_value_scipy=np.nan
        s_corr_scipy=np.nan
        s_value_scipy=np.nan
        me_br=np.nan
        size_obs=np.nan
        size_est=np.nan
    # print ('non bias corrected', me) 
    # print ('non bias removed', me_br) 
    results={'PR':p_corr_scipy,'PR_p_val':p_value_scipy,'RHO':s_corr_scipy,'RHO_p_val':s_value_scipy,
         'NSF':ns_ts2,'RMSE_ts2':rmsd_ts2,'ampl_ts1':ampl_ts1,'me':me,'me_bias_corr':me_br,'ve':ve,
         'size_obs':size_obs, 'size_est':size_est, 
         'obs':datats2_commonts1, 'GDR':datats1_commonts2}
    return results

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
                             df[date_fd].dt.minute, df[date_fd].dt.second)))
    return df

def convert_units(df,height_fd, origin='FEET', to='METER', unit_fd='unit',
                  gauge_fd='gauge_id', stations=None):

    '''Convert from an unit origin that can be ft, cm, or m, to a unit destination that can be
        ft, cm, m
       If the original unit is different from origin convert the unit to 'to' units
       #TODO: Add an exception for the columns, logical control for now
    '''
    conv_dic={'FEET':{'METER':0.3048,'CM':30.48, 'FEET':1},
              'CM': {'METER':0.01,'FEET':0.03281,'CM':1},
              'METER':{'FEET':3.28084,'CM':100,'METER':1}}
    
    #Verify the unit field in df has the valid unit label
    flag = 0
    flag_nc=0
    
    if unit_fd in df.columns:
        if(set(df[unit_fd].unique()).issubset(set(conv_dic.keys()))):
            flag = 1
    else:
        units=origin
        flag_nc=1
        flag = 1
        
    if (height_fd in df.columns) and (gauge_fd in df.columns):
        
        if (origin in conv_dic.keys())and (to in conv_dic.keys())and(flag==1):
            conversion_factor=conv_dic[origin][to]
            look_for=origin
            
            df_final=pd.DataFrame()
            if stations is None:
                stations=df[gauge_fd].unique()

                for lc in stations:

                    df_local=df.loc[df[gauge_fd]==lc].copy()
                    if flag_nc==0:
                        units=df_local[unit_fd].iloc[0]

                    if look_for==units:
                        df_local.loc[:, height_fd+'_'+to]=df_local[height_fd]*conversion_factor
                    elif to==units:
                        df_local.loc[:, height_fd+'_'+to]=df_local[height_fd]
                    else:
                        df_local.loc[:, height_fd+'_'+to]=df_local[height_fd]*conv_dic[units][to]
                    df_final=pd.concat((df_final,df_local), axis=0)
                        # print(df_final.shape, lc+' '+units)
                df_final=df_final.rename(columns={height_fd:'height_rw', height_fd+'_'+to:'height'})
    else: 
        df_final=pd.DataFrame() #This is dangerous
    
    return df_final

#Get statistics for each value 
#Code from https://stackoverflow.com/questions/47637774/pandas-groupby-quantile-values/58535752
def rename(newname):
    def decorator(f):
        f.__name__ = newname
        return f
    return decorator
#Define quatiles at a percentage define by y
def q_at(y):
    @rename(f'q{y:0.3f}')
    def q(x):
        return x.quantile(y)
    return q
#To change names of level col name specially after stats asre calculated with describe 
def change_statfd_names(df,stats_fd, var_fd, additional_fd):
    name_vs=[]
    for vfd in var_fd:
        name_vs=name_vs+[vfd+'_'+sfd for sfd in stats_fd]
    new_colnames=additional_fd+name_vs
    df.columns=new_colnames
    return df

def filter_extreme_duplicates(df_start, st_id, date_fd,height_fd, cols, cut_off, gauge_list=None, exact_date=False):
    '''Function for LOCSS
        It summarizes the median of values that are duplicates at the same date if the standard deviation of those values is lesss than a cut_off in m.
        By default the same date is considered a time difference of 1 min
        Inputs: 
            st_id: id column name that contains the id for gauge_id
            date_fd: date column name
            height_fd: height column name. Values must be in the same units
            cols: columns to preserve in the dataframe 
            cut_off: cut off used to extract the data. 
            gauge_list: (optional) list of gauges/stations to apply the filtering. Default: None. This does the filtering at each station in the st_id column
            exact_date: (optional) If False, the time difference is considered 1 min, if True, the exact date is used
        Outputs:
            datafame without the duplicates that fill out the description. The count column allows to know if the median was calculated
            dataframe with the data removed because if did not comply the cut_off
            The amount of data discarded
    '''
    if gauge_list is None:
        gauge_list=df_start[st_id].unique()
    else:
        gauge_ids=df_start[st_id].unique()
        gauge_list = set.intersection(set(gauge_ids), set(gauge_list))
    
    df_final_np=pd.DataFrame()
    df_final_removed=pd.DataFrame()
    total_discarded=0
    total_duplicated=0
    
    for st in gauge_list:
        
        df=df_start.loc[df_start[st_id]==st].copy()
         
        df=df.sort_values(by=date_fd)
        df['diff_val']=df[height_fd].diff()
        
        # Default 1 min difference 
        
        df['diff_date']=df[date_fd].diff().apply(lambda x: x/np.timedelta64(1, 'm')).fillna(0).astype('int64') #Difference in minutes
        df_dp=pd.concat([df.loc[(df['diff_date']<1).shift(-1).fillna(True)],df.loc[(df['diff_date']<1)&(~df['diff_val'].isnull())]]).sort_values(by=date_fd).copy()
        
        if exact_date:
            date_duplicate_mask=df[date_fd].duplicated(keep=False) #inquire weather the dates are duplicated mark all duplicates at true
            df_dp=df[date_duplicate_mask].copy()

        if not df_dp.empty:
            # print(st)
            df_stats_dp=df_dp[cols+[height_fd]].groupby(cols).describe().reset_index()
            stats_fd=['count','mean','std','min','q_25', 'q_50', 'q_75','max']
            var_fd='height'
            additional_fd=cols
            df_stats_dp=change_statfd_names(df_stats_dp,stats_fd, [var_fd], additional_fd)

            
            discard_mask=~(df_stats_dp[var_fd+'_std'].isnull())&(df_stats_dp[var_fd+'_std']>=cut_off)
            total_discarded=total_discarded+df_stats_dp.loc[discard_mask].shape[0]
            
            # print('Total discarded'+str(total_discarded), df_stats_dp.loc[discard_mask])
            #ic(num_discard)
            #Those that do not satisfy the cutoff 
            df_remove_ex_dp=df_stats_dp.loc[discard_mask].copy()
            # ic(df_remove_ex_dp)
            
            
            
            #Remove the extreme duplicates

            if not df_remove_ex_dp.empty:
                #Remove those with extreme duplicates
                #ic(df.shape)
                # print(st)
                df=df.loc[~df[date_fd].isin(df_remove_ex_dp[date_fd])]
                df_final_removed=pd.concat((df_final_removed, df_remove_ex_dp), axis=0)
                # ic(df.groupby(by=cols, as_index=False).size())

            #for all the duplicates extract the median values. Those that are not duplicates the median is the same value original value
            df=df.groupby(by=cols, as_index=False).agg(height_rc=(height_fd,'median'),
                                                       height_std=(height_fd,'std'),
                                                       height_count=(height_fd,'count'))
                #ic(df.shape)
            
            df_duplicates=df.loc[df['height_count']>=2]
            n_duplicates=df_duplicates['height_count'].sum()-df_duplicates.shape[0]
            total_duplicated=total_duplicated+n_duplicates
            df_final_np=pd.concat((df_final_np, df), axis=0)
                
                
    return df_final_np, df_final_removed, total_discarded, total_duplicated

def extract_data_gauge(df, st_id_fd, st_id, date_fd, height_fd, low_lim, high_lim):
    ''' Extract the rows for a particular station in the dataframe which water elevation change exceed certain threshold 
        between low_lim and high_lim range.
        Inputs:
          df: dataframe with all the values, the station and the height
          st_id_fd: name of the column containing the station id
          st_id: station id in which values would de extracted
          date_fd: name of the column containing the date
          height_fd: name of the column containing the height
          low_lim: lower limit value 
          high_lim: higher limit value
        Output: 
          Three dataframes, two of them with two additional columns.
          diff column, which contains the difference between two consecutive heights, and 
          diff_time, difference among dates in fraction of a day
          first dataframe: dataframe with all the data from the station st_id
          second dataframe: dataframe with all the stations except st_id
          third dataframe: rows to remove which diff values are between the low_lim and high_lim
          
          
    '''
    df_ground_ex=df.loc[df[st_id_fd]==st_id].copy()
    df_ground_re=df.loc[~(df[st_id_fd]==st_id)].copy()
    df_ground_ex=df_ground_ex.sort_values(by=[date_fd])
    df_ground_ex['diff']=df_ground_ex[height_fd].diff()
    df_ground_ex['diff_time']=df_ground_ex[date_fd].diff()#.apply(lambda x: x/np.timedelta64(1, 'D')).fillna(0).astype('int64')
    # df_ground_ex['diff'].max()
    df_diff_to_remove=df_ground_ex.loc[((df_ground_ex['diff']>=high_lim)|(df_ground_ex['diff']<=low_lim))]
    return (df_ground_ex, df_ground_re, df_diff_to_remove)
    #df with only the station st_id, df without station st_id, df with df to remove



def remove_extreme_data(df,df_diff_to_remove, date_fd, height_fd, diff_fd):
    ''' #TO-DO: Do the iteration of this function backwards
        Move a segment of a time series to the start point of another segment 
        #TO-DO: Pending to explain this better
    '''
    df=df.copy()
    row_iterator = df_diff_to_remove.iterrows()
    _, this = next(row_iterator)  # take first item from row_iterator
    for index,next_to in row_iterator:
        current_date=this[date_fd]
        next_date=next_to[date_fd]
        diff_replace=df.loc[df[date_fd]==current_date,diff_fd]
        mask=(df[date_fd]>=current_date)&(df[date_fd]<next_date)
        hg=df.loc[mask,height_fd]
        df.loc[mask,height_fd]=hg-np.round(diff_replace.iloc[0],decimals=2)
        df[diff_fd]=df[height_fd].diff()
        new_diff=np.round(df.loc[df[date_fd]==next_date,diff_fd],decimals=2)
        next_to[diff_fd]=new_diff
        this=next_to

    current_date=next_to[date_fd]
    next_date=df[date_fd].max()
    diff_replace=df.loc[df[date_fd]==current_date,diff_fd]
    # print(diff_replace.values)
    mask=(df[date_fd]>=current_date)&(df[date_fd]<=next_date)
    hg=df.loc[mask,height_fd]
    df.loc[mask,height_fd]=hg-np.round(diff_replace.iloc[0],decimals=2)
    df[diff_fd]=df[height_fd].diff()
    return df

def aggregage_values_by_group(df, group_by_fds, sum_dict):
    ''' Agregate values in a dataframe, leaving one index instead of multiindex as columns
        Inputs:
        df: datatrame to aggregate
        group_by_fds: columns to group by the dataframe and do the aggregation
        sum_dict: desired aggregations for the data 
        
        TODO: Include example
    '''
    df = df.groupby(group_by_fds).aggregate(sum_dict)
    df.columns = ['_'.join([str(c) for c in col]).strip() for col in df.columns.values]
    if isinstance(df, (pd.DatetimeIndex, pd.MultiIndex)):
        df = df.to_frame(index=False)
    return df.reset_index()

def change_df_col_names(df, new_colnames):
    '''Change all the names of a dataframe by the names in new_colnames. The columns are changed in the same order of new_colnames'''
    for nn in new_colnames:
        df.rename(columns={df.columns[new_colnames.index(nn)]:nn}, inplace=True) 
    return df