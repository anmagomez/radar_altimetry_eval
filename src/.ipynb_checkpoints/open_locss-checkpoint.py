#To configure virtual
import pandas as pd
#Test gauges 119, MFN2 and BPN2
#Total  gauges without enough data New york: 9, total from Bangladesh: 26
filter_lakes=["MFN2",
"BPN2",
"BRK2",
"CQK2",
"CHK2",
"CSK2",
"DPK2",
"LNK2",
"LPK2",
"NAK2",
"RLK2",
"119",
"ABB2",
"ADB2",
"BBB2",
"BHB2",
"BNB2",
"BPB2",
"BTB2",
"CCB2",
"CGB2",
"CRB2",
"CTB2",
"DAB2",
"DCB2",
"DHB2",
"DLB2",
"EDB2",
"GBB2",
"JSB2",
"KJB2",
"KSB2",
"KTB2",
"MGB2",
"MKB2",
"PTB2",
"RMB2",
"RUB2",
"SBB2"]

countries={'R2':'France', 'B2':'Bangladesh', 'N2': 'North Carolina'}
df_locss=pd.read_csv('data/readings.csv')
df_locss['loc_id']=df_locss['gauge_id'].str.slice(start=2)


print(df_locss['loc_id'].unique())


