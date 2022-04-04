import pandas as pd
import numpy as np
import os
from LOCSS_data_structure import GaugeCollection

class GroundObservations:
    
    
    def read_ground_data(self, source, file_type='.csv', date_fd=None, time_fd=None,
                         id_fd=None, height_fd=None, station_id=None, path=None, skip_rows=0):
        
        #if source=='USGS':
            #Create function to connect to R and get the station
        
        if source=='LOCSS':
            df=self.get_locss(source,date_fd, time_fd,id_fd, height_fd, station_id)
            return df
        if source=='ARHN':
            df=self.get_arhn(source, date_fd, id_fd, height_fd, station_id, path, skip_rows)
            return df
        if source=='USGS':
            df=self.get_usgs_file(date_fd, id_fd, height_fd, station_id, path, skip_rows)
            df['source']=source
            return df
    
    def get_usgs_file(self, date_fd=None, id_fd=None, height_fd=None, station_id=None, path=None, skip_rows=0):
        if path is None:
            full_path='../data/USGS_data_gages.csv'
        if date_fd is None:
            date_fd='Date'
        if id_fd is None:
            id_fd='site_no'
        if height_fd is None:
            height_fd='X_00065_00003'
        df=pd.read_csv(full_path, parse_dates=[date_fd])
        df[id_fd]=df[id_fd].astype(str).str.strip()
        if station_id is not None:
            if type(station_id)==list:
                df=df.loc[df[id_fd].isin(station_id)]
            else:
                df=df.loc[df[id_fd]==station_id]
        df=self._unify_cols(df, id_fd, date_fd, height_fd)
        return df
        
    def get_arhn(self, source, date_fd=None, id_fd=None, height_fd=None, station_id=None, path=None,skip_rows=0):
        if path is None:
            path='../data/sel_argentina/'
        if date_fd is None:
            date_fd='Fecha y Hora'
        if id_fd is None:
            id_fd='gauge_id'
        if height_fd is None:
            height_fd='Altura [m]'
        
        postfix='.xlsx'
        files = [f for f in os.listdir(path) if (os.path.isfile(os.path.join(path, f))&
                                                  (f.endswith(postfix)))]
        length_id=4 #To extract the id from the name of the station
        stations=[f[:-len(postfix)][-length_id:] for f in files] #Get the ids of the stations
        df_final=pd.DataFrame()
        if station_id is None:
            for f in files:
                df_final=pd.concat((df_final, self.get_one_arhn(path, f, id_fd, station_id, source, date_fd, skip_rows)), axis=0)
        
        elif type(station_id)==list:
            for st in station_id:
                sel_file=[f for f in files if f[:-len(postfix)][-length_id:]==st]
                if len(sel_file)==1:
                    df=self.get_one_arhn(path, sel_file[0], id_fd, st, source, date_fd, skip_rows)
                    df_final=pd.concat((df_final, df), axis=0)
        else:
            sel_file=[f for f in files if f[:-len(postfix)][-length_id:]==station_id]
            df_final=self.get_one_arhn(path, sel_file[0], id_fd, station_id, source, date_fd, skip_rows)
        
        df_final=self._unify_cols(df_final, id_fd, date_fd, height_fd)
        return df_final
              
        
    def get_one_arhn(self, path, f, id_fd, station_id, source, date_fd, skip_rows=0):
        df=pd.read_excel(path+f, skiprows=skip_rows)
        
        df['source']=source
        df[id_fd]=station_id
        df[date_fd]=pd.to_datetime(df[date_fd], dayfirst=True)
        return df
        
    def get_locss(self,source, date_fd=None, time_fd=None,id_fd=None, height_fd=None, station_id=None, all_fd=False):
        #Tested with source only rest values by default
        #For now the data is saved in a fix rout, this has to change later when I get the API
        dir_ts='../data/readings_up_to_20220325.csv'
        dir_loc='../data/gauges_up_to_20220330.csv'
        if date_fd is None: #Tested
            date_fd='date'
            time_fd='time'
        if id_fd is None: #Tested
            id_fd='gauge_id' 
        if height_fd is None: #Tested
            height_fd='height'
        #Load the time series
        df_locss=pd.read_csv(dir_ts, parse_dates=[date_fd])

        #Read the coordinates
        df_coord_locss=pd.read_csv(dir_loc, sep=",")
        df_coord_locss['location']=df_coord_locss[id_fd].str.slice(2,4)

        gc=GaugeCollection()
        df_locss_filtered=gc.filter_test_gages(df_locss,id_fd).copy()
        #Convert field to datetime
        df_locss_filtered[date_fd]=pd.to_datetime((df_locss_filtered[date_fd].astype(str)+' '+df_locss_filtered[time_fd].astype(str)), 
                                                  format='%Y-%m-%d %H:%M:%S')
        if all==False:
            df_locss_filtered=pd.merge(df_locss_filtered, df_coord_locss[[id_fd,'min_height','max_height', 'unit']], on=id_fd)
        else:
            df_locss_filtered=pd.merge(df_locss_filtered, df_coord_locss, on=id_fd)
        df_locss_filtered=df_locss_filtered.loc[(df_locss_filtered[height_fd]>=df_locss_filtered['min_height'])&
                                                (df_locss_filtered[height_fd]<=df_locss_filtered['max_height'])].copy()
        df_locss_filtered['source']=source
        if station_id is not None:
            if type(station_id)==list:
                df_locss_filtered=df_locss_filtered.loc[df_locss_filtered[id_fd].isin(station_id)]
            else:
                df_locss_filtered=df_locss_filtered.loc[df_locss_filtered[id_fd]==station_id]
        
        df_locss_filtered=self._unify_cols(df_locss_filtered, id_fd, date_fd, height_fd)
        return df_locss_filtered
            
        
    def _unify_cols(self,df, id_fd,date_fd,height_fd, filter_others_out=False):
        '''Pending 
        Validate the columns are not already in the dataframe
        Implement filter_out
        '''
        df=df.rename(columns={id_fd:'gauge_id',
                             date_fd:'date',
                             height_fd:'height'})
        return df