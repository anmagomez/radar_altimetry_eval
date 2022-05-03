import pandas as pd
import numpy as np
import os
from LOCSS_data_structure import GaugeCollection

class GroundObservations:
    
    def __init__(self):
        #Initialize all the path variables. Change when services can be stablished directly 
        self.usgs_path='../data/USGS_data_gages.csv'
        self.arhn_path='../data/sel_argentina/'
        self.locss_r_path='../data/readings_up_to_20220502.csv'
        self.locss_g_path='../data/gauges_up_to_20220404.csv'
        self.rvbr_path='../data/brasil/'
    
    def read_ground_data(self, source, file_type='.csv', date_fd=None, time_fd=None,
                         id_fd=None, height_fd=None, station_id=None, path=None, skip_rows=0):
        
        #if source=='USGS':
            #Create function to connect to R and get the station
        df=None
        if source=='LOCSS':
            df=self.get_locss(source,date_fd, time_fd,id_fd, height_fd, station_id)
            
        if source=='ARHN':
            df=self.get_arhn(source, date_fd, id_fd, height_fd, station_id, path, skip_rows)
            
        if source=='USGS':
            df=self.get_usgs_file(date_fd, id_fd, height_fd, station_id, path, skip_rows)
            df['source']=source
            
        if source=='RVBR':
            df=self.get_rvbr_file(source, date_fd, id_fd, height_fd, station_id, path, skip_rows)
        return df
    
    def get_rvbr_file(self, source, date_fd=None, id_fd=None, height_fd=None, station_id=None, path=None, skip_rows=0):
        #Pending deal with code names and encoding
        if path is None:
            path=self.rvbr_path
        if date_fd is None:
            date_fd='Data da Medição'
        if id_fd is None:
            id_fd='Código'
        if height_fd is None:
            height_fd='Cota (m)'
        
        postfix='.xlsx'
        files = [f for f in os.listdir(path) if (os.path.isfile(os.path.join(path, f))&
                                                  (f.endswith(postfix)))]
        #length_id=4 #To extract the id from the name of the station
        #stations=[f[:-len(postfix)][-length_id:] for f in files] #Get the ids of the stations
        
        df_final=pd.DataFrame()
        #df_files=self._get_all_lakes(files)
        start_st=3
        end_st=8
        stations=[f[start_st:end_st] for f in files]
        if station_id is None:
            for st in stations:
                df_final=pd.concat((df_final, self.get_one_arhn(path, files[stations.index(st)], id_fd,
                                                                station_id=st, source=source, 
                                                                date_fd=date_fd, skip_rows=skip_rows)), axis=0)
        elif type(station_id)==list:
            for st in station_id:
                # df_files=self._get_all_lakes(files)
                sel_file=[f for f in files if f[start_st:end_st]==st]
                print(st, sel_file)
                if len(sel_file)==1:
                    df=self.get_one_arhn(path, sel_file[0], id_fd, st, source, date_fd, skip_rows)
                    df_final=pd.concat((df_final, df), axis=0)
        else:
            #sel_file=[f for f in files if f[:-len(postfix)][-length_id:]==station_id] This is faster, but I made it with df, fix later
            sel_file=[f for f in files if f[start_st:end_st]==station_id]
            df_final=self.get_one_arhn(path, sel_file[0], id_fd, station_id, source, date_fd, skip_rows)
        
        df_final[height_fd]=pd.to_numeric(df_final[height_fd].apply(lambda x: str(x).replace(',','.')), errors='coerce')
        df_final=self._unify_cols(df_final, id_fd, date_fd, height_fd)
        
        return df_final
    
    def _get_all_lakes(self, file_list, by=None):
        '''Not using this function for now, faster to deal with the list directly'''
        df=pd.DataFrame(file_list, columns=['file_name'])
        df_ex=df['file_name'].str.split('_', expand=True)

        df_ex=df_ex.rename(columns={0:'location',1:'lake_id',2:'name'})
        df_ex['name']=df_ex['name'].apply(lambda x:x[:-4]) #Take out the xls
        df_all=pd.concat((df,df_ex), axis=1)
        if by=='lake_id':
            return df_all['lake_id']
        elif by=='name':
            return df_all['name']
        else:
            return df_all
    
    def get_usgs_file(self, date_fd=None, id_fd=None, height_fd=None, station_id=None, path=None, skip_rows=0):
        if path is None:
            full_path=self.usgs_path
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
            path=self.arhn_path
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
        #TODO: Solve the engine parsing, it can hide other engine problems
        
        df=pd.read_excel(path+f, skiprows=skip_rows)
        df['source']=source
        df[id_fd]=station_id
        df[date_fd]=pd.to_datetime(df[date_fd], dayfirst=True)
        return df
        
    def get_locss(self,source, date_fd=None, time_fd=None,id_fd=None, height_fd=None, station_id=None, all_fds=False):
        #Tested with source only rest values by default
        #For now the data is saved in a fix rout, this has to change later when I get the API
        dir_ts=self.locss_r_path
        dir_loc=self.locss_g_path
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
        if all_fds==False:
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