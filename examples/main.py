import sys
sys.path.append('../')
import CS_Indicators as CS
from brix import Indicator, Handler
import PreCompOsmNet
import Simulation
import geopandas as gpd
import pandas as pd
import requests
import os
import numpy as np


table_name='epaversion'
#print(sys.path)

r = requests.get('https://cityscope-api.smartaraucania.org/api/table/epaversion/GEOGRID')
print(r.status_code)
# JSON_geogrid =r.json()
# print(JSON_geogrid['features'][0]['geometry']['coordinates'])


geogrid=gpd.read_file('../tables/{}/geogrid.geojson'.format(table_name))
#geogrid=gpd.read_file(r.text, driver='GeoJSON')
zones=gpd.read_file('../tables/{}/zones.geojson'.format(table_name))
zones['GEOID']=zones['GEOID'].astype('int64')
zones=zones.set_index('GEOID')
simpop_df=pd.read_csv('../tables/{}/simpop_df.csv'.format(table_name))

import random
naics_cols=[col for col in zones.columns if 'naics' in col]
simpop_df['naics']=0
for ind, row in simpop_df.iterrows():
    if zones.loc[row['work_geoid'], 'emp_total']>0:
        simpop_df.loc[ind, 'naics']=random.choices(naics_cols, weights=zones.loc[row['work_geoid'], naics_cols].values, k=1)[0]
    else:
        simpop_df.loc[ind, 'naics']=random.choices(naics_cols, k=1)[0]

simpop_df=simpop_df.drop('industry', axis=1)


networks, mode_dicts, route_lengths=PreCompOsmNet.create_2_scale_osmnx_network(
    zones.loc[zones['sim_area']], zones.loc[zones['model_area']],
    add_modes=[{'name': 'walk', 'speed': 4800/3600},
               {'name': 'cycle', 'speed': 14000/3600},
               {'name': 'pt', 'speed': 25000/3600}])
modes={mode: Simulation.Mode(mode_dicts[mode]) for mode in mode_dicts}
mob_sys=Simulation.MobilitySystem(modes=modes, networks=networks)


class Default_Mode_Choice_model():
    def __init__(self):
        pass
    
    def predict_modes(self, all_trips_df):
        options=['drive', 'cycle', 'walk', 'pt']
        all_trips_df['mode']='drive'
        ind_u1500=all_trips_df['route_distance']<1500
        ind_1500_5k=((all_trips_df['route_distance']>1500)&(all_trips_df['route_distance']<5000))
        ind_5k_10k=((all_trips_df['route_distance']>5000)&(all_trips_df['route_distance']<10000))
        ind_10kplus=all_trips_df['route_distance']>10000

        all_trips_df.loc[ind_u1500==True, 'mode']=np.random.choice(
            options, size=sum(ind_u1500), 
            replace=True, p=[0.1, 0.2, 0.6, 0.1])

        all_trips_df.loc[ind_1500_5k==True, 'mode']=np.random.choice(
            options, size=sum(ind_1500_5k), 
            replace=True, p=[0.25, 0.15, 0.3, 0.3])

        all_trips_df.loc[ind_5k_10k==True, 'mode']=np.random.choice(
            options, size=sum(ind_5k_10k), 
            replace=True, p=[0.6, 0.04, 0.01, 0.35])

        all_trips_df.loc[ind_10kplus==True, 'mode']=np.random.choice(
            options, size=sum(ind_10kplus), 
            replace=True, p=[0.9, 0.01, 0.01, 0.08])
        return all_trips_df

mode_descriptions = [{"name": 'drive',
                    'color': "#e41a1d"},
                     {"name": 'cycle',
                    'color': "#377eb8"},
                     {"name": 'walk',
                    'color': "#4daf4a"},
                     {"name": 'pt',
                    'color': "#ffff33"},
                    ]

H=Handler(table_name=table_name)

d=CS.Density_Indicator(zones=zones)
p=CS.Proximity_Indicator(zones=zones, geogrid=geogrid)
m=CS.Mobility_indicator(zones, geogrid, table_name, simpop_df, mob_sys,
                        mode_choice_model=Default_Mode_Choice_model(),
                        mode_descriptions=mode_descriptions,
                        route_lengths=route_lengths)

H.add_indicators([
    d,
    p, 
    m
])



geogrid_data=H.get_geogrid_data()

H.listen()
