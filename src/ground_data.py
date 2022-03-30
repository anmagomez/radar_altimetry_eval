import pandas as pd
import numpy as np
import os
from LOCSS_data_structure import GaugeCollection

class GroundObservations:
    
    
    def read_ground_data(self, source, file_type='.csv', date_fd=None, id_fd=None, 
                         height_fd=None, station_id=None, path=None, skip_rows=0):
        
        #if source=='USGS':
            #Create function to connect to R and get the station
        
        if source=='LOCSS':
            df=self.get_locss(source,date_fd, id_fd, height_fd, station_id)
            return df
        if source=='ARHN':
            df=self.get_arhn(source, date_fd, id_fd, height_fd, station_id, path)
        
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
                df_final=pd.concat(df_final, self.get_one_arhn(path, f, id_fd, station_id, source, date_fd, skip_rows), axis=1)
        
        elif type(station_id)==list:
            for st in station_id:
                sel_file=[f for f in files if f[:-len(postfix)][-length_id:]==st]
                if len(sel_file)==1:
                    df_final=pd.concat(df_final, self.get_one_arhn(path, sel_file[0], id_fd, station_id, source, date_fd, skip_rows), axis=1)
        else:
            sel_file=[f for f in files if f[:-len(postfix)][-length_id:]==station_id]
            df_final=get_one_arhn(path, sel_file[0], id_fd, station_id, source, date_fd, skip_rows)
        return df_final
              
        
    def get_one_arhn(self, path, f, id_fd, station_id, source, date_fd, skip_rows=0):
        df=pd.read_excel(path+f, skiprows=skip_rows, parse_dates=[date_fd])
        df['source']=source
        df[id_fd]=station_id
        df[date_fd]=pd.datetime(df[date_fd], dayfirst=True)
        return df
        
    def get_locss(self,source, date_fd=None, id_fd=None, height_fd=None, station_id=None):
        #Tested with source only rest values by default
        #For now the data is saved in a fix rout, this has to change later when I get the API
        dir_ts='../data/readings_up_to_20220325.csv'
        dir_loc='../data/gauges_up_to_20220330.csv'
        if date_fd is None: #Tested
            date_fd='date'
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
        df_locss_filtered=gc.filter_test_gages(df_locss,id_fd)
        df_locss_filtered=pd.merge(df_locss_filtered, df_coord_locss[[id_fd,'min_height','max_height', 'unit']], on=id_fd)
        df_locss_filtered=df_locss_filtered.loc[(df_locss_filtered[height_fd]>=df_locss_filtered['min_height'])&
                                                (df_locss_filtered[height_fd]<=df_locss_filtered['max_height'])].copy()
        df_locss_filtered['source']=source
        if station_id is not None:
            if type(station_id)==list:
                df_locss_filtered=df_locss_filtered.loc[df_locss_filtered[id_fd].isin(station_id)]
            else:
                df_locss_filtered=df_locss_filtered.loc[df_locss_filtered[id_fd]==station_id]
        
        return df_locss_filtered
            
        
        