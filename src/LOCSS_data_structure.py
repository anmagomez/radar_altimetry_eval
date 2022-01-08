import pandas as pd


class Coordinate: 
    
    """Coordinate constructor
       Get the lattitude, lat, longitude, long, and Coordinate Reference System, CRS. CRS must be a number with the EPSG identifcator. 
       EPSG identificator can be retrieve from https://spatialreference.org/
       By default CRS=4326"""
    
    def __init__ (self, lat, long, CRS=4326): 
        self.lat=lat
        self.lon=long
        self.CRS=CRS
        
class Measurement: 
    
    def __init__ (self, date, value, type_value): 
        self.date=date
        self.value=lat
        self.type_value=type_value
        
class Elevation:
    
    def __init__(self, elev, local_GPS, source):
        self.elevation=elev
        self.local_GPS=local_GPS
        self.source=source
class Note:
    
    def __init__(self, comment, date):
        self.comment=comment
        self.date=date
    
class Gauge:
    
    def __init__(self, name, guage_id, install_date, photo, unit, lat, long, CRS, elevation, local_GPS, elev_source,notes):
        self.name=name
        self.visible_id=guage_id
        self.install_date=install_date
        self.photo=photo
        self.unit=unit
        self.coord=Coordinate(lat, long, CRS)
        self.elev=Elevation(elevation, local_GPS, elev_source)
        
        if isinstance(notes, str):
            #assume the notes corresponds to the installation date
            self.note=Note(notes, install_date)
            
        if isinstance(notes, Note) or isinstance(notes, List) :
            self.note=notes
        
        self.country=self.get_country(gauge_id)
        self.subregion=self.get_country(gauge_id)
   
    def set_id(self, gauge_id):
        #cod format LakeInitials (2 letters), _, country (2 letters), 
        self.country=self.get_country(gauge_id)
        self.subregion=self.get_country(gauge_id)
        return g_id

class GaugeCollection:
    test_gages=["MFN2","BPN2","119"]
    
    def __init__(self, gages_collection, fd_id):
        self.df=gages_collection
        self.fd_id=fd_id
        
    def __init__(self):
        self.df=pd.DataFrame()
    
    def filter_test_gages(self, df = None, fd_id=None):
        """Filter LOCSS testing gages first. This method doesn't modify the internal gage collection. 
            df is the gauge list. If it is not provided the internal df attribute is used """
        if df is None:
            df=self.df
            fd_id=self.fd_id
            
        df=df.loc[~df[fd_id].isin(self.test_gages)]
        return df
        
    #def add_gauge(self, name, guage_id, install_date, photo, unit, lat, long, CRS, elevation, local_GPS, elev_source,notes):
    
    #def add_gages_from_loccs_df(self, df_locss):
        
    def add_gages(self, gages): #gages must be a dataframe
        self.df=df.append(gages)
        