from brix import Indicator, Handler
import Simulation
import PreCompOsmNet
import statsmodels
from statsmodels.distributions.empirical_distribution import ECDF
import geopandas as gpd
import urllib
import json
import requests
import math
import osmnx
from shapely.geometry import Point, shape
import datetime
from scipy import spatial
import numpy as np
import sys
import pandas as pd
import copy
import random

"""
"""

def aggregate_types_over_grid(geogrid_data, side_length, type_def, interactive_only=True):
    cell_area=side_length*side_length
    aggregated={}
    for cell in geogrid_data:
        if ((cell['interactive']=='Web') or (interactive_only==False)):
            name=cell['name']
            type_info=type_def[name]
            height=cell['height']
            if isinstance(height, list):
                height=height[-1]
                if 'sqm_pperson' in type_info:
                    sqm_pperson=type_info['sqm_pperson']
                else:sqm_pperson=50
                total_capacity=height*cell_area/sqm_pperson
                if name in aggregated:
                    aggregated[name]+=total_capacity
                else:
                    aggregated[name]=total_capacity
    return aggregated

def aggregate_attributes_over_grid(agg_types, attribute, side_length, type_def, digits=None):
    # TODO: eliminate repetition with previous function
    cell_area=side_length*side_length
    aggregated={}
    for cell_type in agg_types:
        type_info=type_def[cell_type]
        if attribute in type_info:
            if type_info[attribute] is not None:
                total_capacity=agg_types[cell_type]
                for code in type_info[attribute]:
                    if digits==None:
                        code_digits=code
                    else:
                        code_digits=code[0:digits]
                    attr_capacity=total_capacity*type_info[attribute][code]
                    if code_digits in aggregated:
                        aggregated[code_digits]+= attr_capacity
                    else:
                        aggregated[code_digits]= attr_capacity
    return aggregated

def get_central_nodes(geodf, G):
    """ takes a geodf and returns a list of nodes closest to their centroids
    returns both the nearest nodes and the distance
    """
    geodf_projected=osmnx.projection.project_gdf(geodf)
    projected_centroids=geodf_projected['geometry'].centroid.values
    projected_centroids_lst=[[c.x, c.y] for c in projected_centroids]
    G_proj=osmnx.projection.project_graph(G, geodf_projected.crs)
    G_proj_coordinates={node: [G_proj.nodes[node]['x'], G_proj.nodes[node]['y']] for node in G_proj.nodes}
    node_order=[node for node in G_proj_coordinates]
    nodes_kdtree=spatial.KDTree([G_proj_coordinates[node] for node in node_order])
    dist, ind_nearest=nodes_kdtree.query(projected_centroids_lst)
    nearest_nodes=[node_order[ind] for ind in ind_nearest]
    return nearest_nodes, dist

def prob_floor(num):
    result=int(num)
    remainder=num-result
    if random.uniform(0, 1)<remainder:
        result+=1
    return result

class Proximity_Indicator(Indicator):
    def setup(self, zones, geogrid, buffer=1200, pois=None, poi_names=[], score_dict=None, 
              range_factor=2, ref_geogrid_data=None):
        print('Setting up Proximity Indicator')
        self.name='Proximity'
        self.indicator_type = 'hybrid'
        self.zone_to_node_tolerance=500
        self.grid_to_node_tolerance=100
        self.naics_codes=[col for col in zones.columns if 'naics' in col]
        self.buffer=buffer
        self.zones=zones
        self.geogrid=geogrid
        self.sim_area_geoids=list(zones.loc[zones['sim_area']].index)
        self.all_site_ids=self.sim_area_geoids+list(range(len(geogrid))) 
        self.get_graph_reference_area()
       
        print('\t Getting central nodes')
        zones_nearest_nodes, zones_nearest_dist= get_central_nodes(self.zones, self.ref_G)
        self.zones['central_node']=zones_nearest_nodes
        self.zones['nearest_dist']=zones_nearest_dist
        
        grid_nearest_nodes, grid_nearest_dist= get_central_nodes(self.geogrid, self.ref_G)
        self.geogrid['central_node']=grid_nearest_nodes
        self.geogrid['nearest_dist']=grid_nearest_dist

        self.pois=pois
        self.poi_names=poi_names
        if pois is not None:
            poi_nearest_nodes, poi_nearest_dist= get_central_nodes(self.pois, self.ref_G)
            self.pois['central_node']=poi_nearest_nodes
            self.pois['nearest_dist']=poi_nearest_dist

        self.score_dict=score_dict
        if self.score_dict==None:
            self.score_dict={'walkable_housing': {'col':'res_total', 'from': 'source_emp'},
                            'walkable_employment': {'col':'emp_total', 'from': 'source_res'},
                            'walkable_healthcare': {'col':'emp_naics_62', 'from': 'source_res'},
                            'walkable_hospitality': {'col':'emp_naics_72', 'from': 'source_res'},
                            'walkable_shopping': {'col':'emp_naics_44-45', 'from': 'source_res'}}
            for name in poi_names:
                self.score_dict[name]={'walkable_{}'.format(name): {'col': name, 'from': 'source_res'}}

        for score in self.score_dict:
            if self.score_dict[score]['col'] not in self.zones.columns:
                self.zones[self.score_dict[score]['col']]=0

        self.get_reachable_geoms_from_all()
        self.calculate_baseline_scores(range_factor)
        self.calculate_baseline_ind(ref_geogrid_data)
            
            
    def make_ego_graph_around_geometry(self, zone, tolerance):
        """
        For geometries within the buffered geogrid only.
        Returns the graph within a walkable distance of the centre of the zone
        """
        if zone['nearest_dist']<tolerance:
            sub_graph = osmnx.graph.nx.ego_graph(self.ref_G, zone['central_node'], radius=1200, distance='length')
        else:
            sub_graph = osmnx.graph.nx.Graph()
        return sub_graph
    
    def get_graph_reference_area(self):
        reference_zones=self.zones.loc[self.zones['reference_area']]
        print('\t Downloading graph for reference area')
        reference_zone_graph=self.get_network_around_geom_buffered(reference_zones)
        self.ref_G=reference_zone_graph
        
    def calculate_baseline_scores(self, range_factor):
        # TODO: should use the get_reachable_geoms function?
        print('\t Calculating baseline scores')
#         self.base_scores={'walkable_{}'.format(x): [] for x in [
#             'employment', 'housing', 'healthcare', 'hospitality', 'shopping']}
        base_scores={}
        self.base_reachable_attributes={}
        self.score_ranges={}
        stats_to_aggregate=[col for col in self.zones.columns if (('res_' in col) or ('emp_' in col))]
        # get the baseline reachable attributes and scores for every zone
        for ind, row in self.zones.loc[self.zones['reference_area']].iterrows():
            reachable_zones=self.zone_to_reachable[ind]['zones']
            self.base_reachable_attributes[ind]=self.zones.loc[reachable_zones][stats_to_aggregate].sum().to_dict()
            self.base_reachable_attributes[ind]['source_res']=row['res_total']
            self.base_reachable_attributes[ind]['source_emp']=row['emp_total']
            # Add the aggregated point POIs if there are any
            if self.pois is not None:
                reachable_pois=self.zone_to_reachable[ind]['pois']
                agg_reachable_poi=self.pois.loc[reachable_pois].sum()
                for name in self.poi_names:
                    self.base_reachable_attributes[ind][name]=agg_reachable_poi[name]
            # get scores for individual zones- weighting cancels out
            base_scores[ind]=self.attributes_to_scores([self.base_reachable_attributes[ind]])
            
        # Create the ranges for each score using only the zones (not the grid cells
        # TODO: avoid repetition with thr Density indicator
        self.base_zones_scores=pd.DataFrame.from_dict(base_scores, orient='index')
        for score in self.base_zones_scores.columns:
            base_scores_no_nan=[x for x in self.base_zones_scores[score] if x==x]
            min_range, max_range=min(base_scores_no_nan), max(base_scores_no_nan)
            max_range=min_range+range_factor*(max_range-min_range)
            diff=max((1, max_range-min_range))
            self.score_ranges[score]={'min': min_range,
                                     'max': max_range,
                                     'range': diff}
            
        # get weighted scores across the simulation area zones 
        # (ignore the grid which is empty in reference and therefore would be weighted zero)
#         ref_scores=self.attributes_to_scores([self.base_reachable_attributes[ind] for ind in self.sim_area_geoids])
#         self.ref_ind=self.normalise_ind(ref_scores)
            
        # get the base reachable attributes for every grid cell location
        for i_c in range(len(self.geogrid)):
            reachable_zones=self.grid_to_reachable[i_c]['zones']
            self.base_reachable_attributes[i_c]=self.zones.loc[reachable_zones][stats_to_aggregate].sum().to_dict()
            if self.pois is not None:
                reachable_pois=self.grid_to_reachable[i_c]['pois']
                agg_reachable_poi=self.pois.loc[reachable_pois].sum()
                for name in self.poi_names:
                    print(self.pois.names)
                    print(i_c.keys())
                    
                    self.base_reachable_attributes[i_c][name]=agg_reachable_poi[name]
                    print(self.base_reachable_attributes.keys())
                    
    def calculate_baseline_ind(self, ref_geogrid_data):
        if ref_geogrid_data is None:
            ref_scores=self.attributes_to_scores([self.base_reachable_attributes[ind] for ind in self.sim_area_geoids])
        else:
            ref_reachable_attributes=self.get_new_reachable_attributes(ref_geogrid_data)
            ref_reachable_attributes_site= [ref_reachable_attributes[ind] for ind in self.all_site_ids]
            ref_scores=self.attributes_to_scores(ref_reachable_attributes_site)
        self.ref_ind=self.normalise_ind(ref_scores)
            
    def get_reachable_geoms(self, zone, tolerance):
        """
        find all grid cells, zones and point POIs reachable from a geometry
        """
        sub_graph=self.make_ego_graph_around_geometry(zone, tolerance)
        sub_graph_nodes=sub_graph.nodes(data=False)
        reachable_zones= list(self.zones.loc[
                ((self.zones['central_node'].isin(list(sub_graph_nodes)))&
                 (self.zones['nearest_dist']<self.zone_to_node_tolerance))
                ].index.values)
        reachable_grid_cells=list(self.geogrid.loc[
                ((self.geogrid['interactive']=='Web')&
                    (self.geogrid['central_node'].isin(list(sub_graph_nodes)))&
                    (self.geogrid['nearest_dist']<self.grid_to_node_tolerance))
                ].index.values)
        reachable={'zones': reachable_zones, 'cells': reachable_grid_cells}
        if self.pois is not None:
            reachable_pois=list(self.pois.loc[
                (
                    (self.pois['central_node'].isin(list(sub_graph_nodes)))&
                    (self.pois['nearest_dist']<self.grid_to_node_tolerance))
                ].index.values)
            reachable['pois']= reachable_pois   	
        return reachable
    
    def get_reachable_geoms_from_all(self):
        """
        For every grid cell and every reference area zone:
        	find the reachable zones,  grid cells and POIs
        """
        print('\t Finding all reachable geometries from each geometry')
        self.grid_to_reachable, self.zone_to_reachable={}, {}
        for ind, row in self.geogrid.iterrows():
            self.grid_to_reachable[ind]=self.get_reachable_geoms(row, self.grid_to_node_tolerance)
        for ind, row in self.zones.loc[self.zones['reference_area']].iterrows():
            self.zone_to_reachable[ind]=self.get_reachable_geoms(row, self.zone_to_node_tolerance)
        # create a reverse lookup to map from each grid cell to the cells from which it is reachable
        self.grid_to_reverse_reachable={}
        for i, row_i in self.geogrid.iterrows():
            self.grid_to_reverse_reachable[i]={'zones': [], 'cells': []}
            for j, row_j in self.geogrid.iterrows():
                if i in self.grid_to_reachable[j]['cells']:
                    self.grid_to_reverse_reachable[i]['cells'].append(j)
            for ind_z, row_z in self.zones.loc[self.zones['reference_area']].iterrows():
                if i in self.zone_to_reachable[ind_z]['cells']:
                    self.grid_to_reverse_reachable[i]['zones'].append(ind_z)                    
            
    def attributes_to_scores(self, attributes):  
    	# TODO: more flexibly caclulate scores based on supplied data      
        totals={}
        totals['source_res']=sum([s['source_res'] for s in attributes])
        totals['source_emp']=sum([s['source_emp'] for s in attributes])

        scores={}
        for score_name in self.score_dict:
            if totals[self.score_dict[score_name]['from']]==0:
                scores[score_name]=0
            else:
                scores[score_name]=sum([s[self.score_dict[score_name]['from']]*s[self.score_dict[score_name]['col']] for s in attributes])/totals[self.score_dict[score_name]['from']]
        return scores
    
    def return_indicator(self, geogrid_data):
        start_ind_calc=datetime.datetime.now()
        # make copy of base_scores
        new_reachable_attributes=self.get_new_reachable_attributes(geogrid_data)
        new_reachable_attributes_site= [new_reachable_attributes[ind] for ind in self.all_site_ids]
        new_scores=self.attributes_to_scores(new_reachable_attributes_site)
        new_ind=self.normalise_ind(new_scores)
        outputs=[]
        for ind_name in new_scores:
            outputs.append({'name': ind_name ,
                            'value': new_ind[ind_name],
                            'raw_value': new_scores[ind_name],
                            'ref_value': self.ref_ind[ind_name],
                            'viz_type': self.viz_type})
        end_ind_calc=datetime.datetime.now()
        
        new_reachable_attributes_grid= [new_reachable_attributes[i_c] for i_c in range(len(geogrid_data))]
        heatmap=self.compute_heatmaps(new_reachable_attributes_grid)
        
        end_hm_calc=datetime.datetime.now()
        
        print('Prox Ind: {}'.format(end_ind_calc-start_ind_calc))
        print('Prox HM: {}'.format(end_hm_calc-end_ind_calc))
        
        return {'heatmap':heatmap,'numeric':outputs}
    
    
    def get_new_reachable_attributes(self, geogrid_data):
        new_reachable_attributes=copy.deepcopy(self.base_reachable_attributes) # the updated REACHABLE attributes for every zone and cell
        side_length=geogrid_data.get_geogrid_props()['header']['cellSize']
        cell_area=side_length*side_length
        type_def=geogrid_data.get_type_info()
        # For each cell
            # calculate the new attributes at this cell
            # increment the source attributes for this cell
            # increment the reachable attributes for every geometry which can reach this cell
        for i_c, cell in enumerate(geogrid_data):
            name=cell['name']
            height=cell['height']
            if isinstance(height, list):
                height=height[-1]
            type_info = type_def[name] 
            if 'sqm_pperson' in type_info:
                sqm_pperson=type_info['sqm_pperson']
            else:
                sqm_pperson=50
            total_capacity=height*cell_area/sqm_pperson
            agg_naics=aggregate_attributes_over_grid(
                {name: total_capacity}, 'NAICS', side_length, type_def, digits=2)
            agg_lbcs=aggregate_attributes_over_grid(
                {name: total_capacity}, 'LBCS', side_length, type_def, digits=1)

            added_attributes={} # the new attributes of THIS cell
            # For point POIs: 1 floor = 1 point
            if name in self.poi_names:
                added_attributes[name]=height
            cell_employment=sum(agg_naics.values())
            #print(i_c)
            #print(new_reachable_attributes.keys())
            #if i_c < 165:
            new_reachable_attributes[i_c]['source_emp']=cell_employment
            added_attributes['emp_total']=cell_employment

            if '1' in agg_lbcs:
                cell_population=agg_lbcs['1']
            else:
                cell_population=0
            #if i_c < 165:
            added_attributes['res_total']=cell_population
            new_reachable_attributes[i_c]['source_res']=cell_population

            for combined_code in self.naics_codes:
                naics_codes=combined_code.split('naics_')[1].split('-')
                for code in naics_codes:
                    if code in agg_naics:
                        if code in added_attributes:
                            added_attributes[combined_code]+=agg_naics[code]
                        else:
                            added_attributes[combined_code]=agg_naics[code]
            # the newly attributes are added to every zone and cell which can reach this cell   
            reverse_reachable=self.grid_to_reverse_reachable[i_c]
            for ind_z in reverse_reachable['zones']:
                for attr in added_attributes:
                    new_reachable_attributes[ind_z][attr]+=added_attributes[attr]
            for j_c in reverse_reachable['cells']:
                for attr in added_attributes:
                    new_reachable_attributes[j_c][attr]+=added_attributes[attr] 
        return new_reachable_attributes

    def normalise_ind(self, raw_ind):
        norm_ind={}
        for ind_name in raw_ind:
            norm_ind[ind_name]=min(1, (raw_ind[ind_name] - self.score_ranges[ind_name]['min'])/ self.score_ranges[ind_name]['range'])     
        return norm_ind

    def compute_heatmaps(self, grid_reachable_area_stats):
        features=[]
        heatmap={'type': 'FeatureCollection',
                 'properties': [c.split('_')[1] for c in self.score_dict]}
        x_centroid_list, y_centroid_list=self.geogrid['x_centroid'], self.geogrid['y_centroid']
        for i_c, cell_stats in enumerate(grid_reachable_area_stats):
            features.append({
              "type": "Feature",
              "properties": [min(1, (cell_stats[self.score_dict[s]['col']] - self.score_ranges[s]['min'])/ self.score_ranges[s]['range']) for s in self.score_dict],
              "geometry": {
                "type": "Point",
                "coordinates": [
                  x_centroid_list[i_c],
                  y_centroid_list[i_c]
                ]
              }
            })
        heatmap['features']=features
        return heatmap
      
    def get_network_around_geom_buffered(self, geom):
        """
        Creates a buffer around the geometry and downloads the graph for this area
        """
        geom_projected=osmnx.projection.project_gdf(geom)
        geom_projected_buffered=geom_projected.unary_union.buffer(self.buffer)

        geom_projected_buffered_gdf=gpd.GeoDataFrame(geometry=[geom_projected_buffered], crs=geom_projected.crs)
        geom_wgs_buffered_gdf=geom_projected_buffered_gdf.to_crs(geom.crs) 
        
        return osmnx.graph.graph_from_polygon(geom_wgs_buffered_gdf.iloc[0]['geometry'], network_type='walk')

class Density_Indicator(Indicator):
    def setup(self, zones, ref_geogrid_data=None, range_factor=2):
        self.name='Density'
        self.zones=zones
        self.sim_area_geoids=list(zones.loc[zones['sim_area']].index)
        self.grid_cell_area=None
        self.compute_base_score_distribution(range_factor)
        self.compute_base_ind(ref_geogrid_data)
                
    def compute_base_score_distribution(self, range_factor):
        """
        computes ranges of the indicators across the static region
        the ranges are later used to normalise indicators
        """
        self.score_ranges={}
        self.base_scores={}
        self.base_scores['res_density']=self.zones.apply(lambda row: self.res_density(row), axis=1)
        self.base_scores['emp_density']=self.zones.apply(lambda row: self.emp_density(row), axis=1)
        self.base_scores['live_work_score']=self.zones.apply(
            lambda row: self.get_live_work_score(row), axis=1)
        
        # Diversity
        industry_columns=[col for col in self.zones.columns if 'emp_naics' in col]
        res_income_columns=[col for col in self.zones.columns if 'res_income' in col]
        
        self.base_scores['industry_diversity']=self.zones.apply(
            lambda row: self.get_diversity(row, species_columns=industry_columns), axis=1)
        self.base_scores['income_diversity']=self.zones.apply(
            lambda row: self.get_diversity(row, species_columns=res_income_columns), axis=1)
        
        for score in self.base_scores:
            base_scores_no_nan=[x for x in self.base_scores[score] if x==x]
            min_range, max_range=min(base_scores_no_nan), max(base_scores_no_nan)
            max_range=min_range+range_factor*(max_range-min_range)
            diff=max((1, max_range-min_range))
            self.score_ranges[score]={'min': min_range,
                                     'max': max_range,
                                     'range': diff}

    def compute_base_ind(self, ref_geogrid_data):
        base_site_stats=self.combine_site_attributes(ref_geogrid_data)
        self.base_ind=self.calculate_indicators(base_site_stats)
        
    def combine_site_attributes(self, geogrid_data=None):
        """
        takes the attributes of the geogrid_data (new programming) and 
        the zones which overlap with the geogrid and (pre-existing programming)
        aggregates them together to get the updated site stats
        """
        stats_to_aggregate=['area']+[col for col in self.zones.columns if (('res_' in col) or ('emp_' in col))]
        temp_site_stats=dict(self.zones.loc[self.sim_area_geoids, 
                                                 stats_to_aggregate].sum())
        if geogrid_data is not None:
            side_length=geogrid_data.get_geogrid_props()['header']['cellSize']
            type_def=geogrid_data.get_type_info()
            agg_types=aggregate_types_over_grid(geogrid_data, side_length=side_length, type_def=type_def)
            agg_naics=aggregate_attributes_over_grid(agg_types, 'NAICS', side_length=side_length, type_def=type_def, digits=2)
            agg_res=aggregate_attributes_over_grid(agg_types, 'res_income', side_length=side_length, type_def=type_def)
            
            # update total residential and total employment
            add_emp=sum(agg_naics.values())
            add_res=sum(agg_res.values())   
            temp_site_stats['res_total']+=add_res
            temp_site_stats['emp_total']+=add_emp
            
            # update employment for each NAICS code
            for col in temp_site_stats:
                if 'naics' in col:
                    # if this is a double naics code column (eg. 44-45), add the new employment for both 44 and 45
                    col_naics_codes=col.split('naics_')[1].split('-')
                    for code in col_naics_codes:
                        if code in agg_naics:
                            temp_site_stats[col]+=agg_naics[code]  
                            
            # update residents for each income level
            for col in temp_site_stats:
                if 'res_income_' in col:
                    income_level=col.split('res_income_')[1]
                    if income_level in agg_res:
                        temp_site_stats[col]+=agg_res[income_level]                    
        return temp_site_stats

    
    def calculate_indicators(self, site_stats):
        raw_ind={}
        raw_ind['res_density']=self.res_density(site_stats)
        raw_ind['emp_density']=self.emp_density(site_stats)
        raw_ind['live_work_score']=self.get_live_work_score(site_stats)
        
        industry_columns=[col for col in self.zones.columns if 'emp_naics' in col]
        res_income_columns=[col for col in self.zones.columns if 'res_income' in col]
        
        raw_ind['industry_diversity']=self.get_diversity(site_stats, species_columns=industry_columns)
        raw_ind['income_diversity']=self.get_diversity(site_stats, species_columns=res_income_columns)
               
        norm_ind={}
        for ind_name in raw_ind:
            if 'density' in ind_name:
                norm_ind[ind_name]=min(1, (raw_ind[ind_name] - self.score_ranges[ind_name]['min'])/ self.score_ranges[ind_name]['range']) 
            else:
                norm_ind[ind_name]=raw_ind[ind_name]
        return {'raw': raw_ind, 'norm': norm_ind}
                  
    def return_indicator(self, geogrid_data):
        start_ind_calc=datetime.datetime.now()
        new_site_stats=self.combine_site_attributes(geogrid_data=geogrid_data)
        new_ind=self.calculate_indicators(new_site_stats)
        
#         base_site_stats=self.combine_site_attributes(geogrid_data=None)
#         base_ind=self.calculate_indicators(base_site_stats)
        
        outputs=[]
        for ind_name in new_ind['raw']:
            outputs.append({'name': ind_name.replace('_', ' ').title(),
                           'raw_value': new_ind['raw'][ind_name],
                           'value': new_ind['norm'][ind_name],
                           'ref_value': self.base_ind['norm'][ind_name]})
        end_ind_calc=datetime.datetime.now()
        print('Dens Ind: {}'.format(end_ind_calc-start_ind_calc))
#         print(outputs)
        return outputs
    
    @staticmethod
    def res_density(obj):
        """
        input can be a row if the baseline geodataframe
        or a dict representing a dynamic site
        """
        if obj['area']>0:
            return obj['res_total']/obj['area']
        return 0
    
    @staticmethod
    def emp_density(obj):
        """
        input can be a row if the baseline geodataframe
        or a dict representing a dynamic site
        """
        if obj['area']>0:
            return obj['emp_total']/obj['area'] 
        return 0
    
    @staticmethod
    def get_live_work_score(obj):
        if obj['emp_total']*obj['res_total']==0:
            return 0
        if obj['emp_total'] < obj['res_total']:
            return obj['emp_total']/obj['res_total']
        else:
            return obj['res_total']/obj['emp_total']
     
    @staticmethod
    def get_diversity(obj, species_columns):
        species_counts=[obj[col] for col in species_columns]
        diversity=0
        pop_size=sum(species_counts)
        if ((len(species_counts)>1) and (pop_size>0)):        
            for count in species_counts:
                pj=count/pop_size
                if not pj==0:
                    diversity+= -pj*math.log(pj)
            equitability=diversity/math.log(len(species_counts))
            return equitability
        else:
            return float('nan')


class Mobility_indicator(Indicator):
    def setup(self, zones, geogrid, table_name, simpop_df, mob_sys,
             mode_descriptions=None, profile_descriptions=None,
             simpop_sample_frac=1, N_max=250, mode_choice_model=None,
             route_lengths=None):
        self.N_max=N_max
        self.simpop_sample_frac=simpop_sample_frac
        self.name='Mobility'
        self.create_trip_descriptions(mode_descriptions, profile_descriptions)
        self.geogrid=geogrid
        self.zones=zones
        self.table_name=table_name
        # TODO: if route lengths is None (not precomputed), should get this by calling a method of the mob_sys
        self.route_lengths=route_lengths
        self.base_simpop_df=simpop_df.copy()
        print('Init simulation111')
        grid_zones=self.create_grid_zones()
        model_zones=zones.loc[zones['model_area']]
        model_zones['grid_area']=False
        combined_zones=model_zones.append(grid_zones)
        sim_geoids=list(zones.loc[zones['sim_area']].index)+list(grid_zones.index)
        self.sim=Simulation.Simulation(simpop_df, mob_sys, combined_zones, sim_geoids=sim_geoids,
            mode_descriptions=self.mode_descriptions, profile_descriptions=self.profile_descriptions)
        if mode_choice_model is not None:
            self.sim.set_choice_models(mode_chooser=mode_choice_model)
        self.create_zone_dist_mat()
    
    def create_trip_descriptions(self, mode_descriptions, profile_descriptions):
        if mode_descriptions is None:
            mode_descriptions = [{"name": 'drive',
                                'color': "#e41a1d"},
                                 {"name": 'cycle',
                                'color': "#377eb8"},
                                 {"name": 'walk',
                                'color': "#4daf4a"},
                                 {"name": 'pt',
                                'color': "#ffff33"},
                                ]
        if profile_descriptions is None:
            profile_descriptions = [{"name": 'low',
                                'color': "#7fc97f"},
                                 {"name": 'mid',
                                'color': "#beaed4"},
                                 {"name": 'high',
                                'color': "#fdc086"},
                                ]
        self.mode_descriptions=mode_descriptions
        self.profile_descriptions=profile_descriptions

    def create_grid_zones(self):
        print('Se están creando las zonas')
        grid_zones=self.geogrid.copy()
        centroids=grid_zones['geometry'].centroid
        #for c in centroids:
            #print(c)
        grid_zones['x_centroid']=[c.x for c in centroids]
        grid_zones['y_centroid']=[c.y for c in centroids]
        cols_to_init=[col for col in self.zones.columns if (
            ('emp_' in col) or ('res_' in col))]
        for col in cols_to_init:
            grid_zones[col]=0
        for area_col in ['model_area', 'sim_area', 'grid_area']:
            grid_zones[area_col]=1
        return grid_zones

        
    def create_zone_dist_mat(self):
        print('Computing zone-zone dist mat')
        zone_index=self.sim.zones.index
        self.dist_mat={}
        chosen_nodes_ls=[pn[0] for pn in self.sim.zones['possible_nodes_drive']]
        for i in range(len(self.sim.zones)):
            self.dist_mat[zone_index[i]]={}
            for j in range(len(self.sim.zones)):
                from_node=chosen_nodes_ls[i]
                to_node=chosen_nodes_ls[j]
                if from_node==to_node:
                    self.dist_mat[zone_index[i]][zone_index[j]]=50
                else:
                    self.dist_mat[zone_index[i]][zone_index[j]]=self.route_lengths[from_node][to_node] 
        
    
    def geogrid_updates(self, geogrid_data):
        new_simpop=[]
        cols_to_zero= [col for col in self.sim.zones.columns if (
            ('emp_' in col) or ('res_' in col))]
        res_cols=[col for col in self.sim.zones.columns if'res_' in col]
        # Initialise vacant residential assuming 4% vacancy rate of existing housing
        zones_copy=self.sim.zones.copy()
        zones_copy[res_cols]*=0.04
        zones_copy.loc[zones_copy.grid_area==True, cols_to_zero]=0

        side_length=geogrid_data.get_geogrid_props()['header']['cellSize']
        type_def=geogrid_data.get_type_info()
        for i_c, cell in enumerate(geogrid_data):
            name=cell['name']
            type_info=type_def[name]
            if ((not name =='None') and (not cell['interactive']==False)):
                height=cell['height']
                cell_area=side_length*side_length
                if isinstance(height, list):
                    height=height[-1]
                if 'sqm_pperson' in type_info:
                    sqm_pperson=type_info['sqm_pperson']
                else:
                    sqm_pperson=50
                total_capacity=height*cell_area/sqm_pperson
                # update where the residential capacity exists
                if ('res_income' in type_info):
                    if type_info['res_income'] is not None:
                        zones_copy.loc[i_c, 'res_total']=total_capacity
                        for income_level in type_info['res_income']:
                            zones_copy.loc[i_c, 'res_income_{}'.format(income_level)
                                              ]=type_info['res_income'][income_level]*total_capacity
        for i_c, cell in enumerate(geogrid_data):
            name=cell['name']
            type_info=type_def[name]
            if not name =='None':
                height=cell['height']
                cell_area=side_length*side_length
                if isinstance(height, list):
                    height=height[-1]
                if 'sqm_pperson' in type_info:
                    sqm_pperson=type_info['sqm_pperson']
                else:
                    sqm_pperson=50
                total_capacity=self.simpop_sample_frac*height*cell_area/sqm_pperson
                # update the agents
                if type_info['NAICS'] is not None:
                    workers={code: prob_floor(type_info['NAICS'][code]*total_capacity) for code in  type_info['NAICS']}   
                    for code in workers:
                        # TODO: choose the income level based on a mapping from NAICS to income level
                        # estimate from the available data
                        profile=random.sample([p['name'] for p in self.profile_descriptions], 1)[0]
                        home_locations=self.sample_home_locations(zones_copy, i_c, profile, n=workers[code])
                        for i_w in range(workers[code]):
                            new_simpop.append({'work_geoid': i_c,'home_geoid': home_locations[i_w],
                                               'naics': code, 'earnings': profile,
                                              'age': 40})

        return new_simpop
            
    def sample_home_locations(self, zones_copy, work_geoid, earnings, n, beta=1):
        attraction=zones_copy['res_income_{}'.format(earnings)]
        impedance=[self.dist_mat[hid][work_geoid] for hid in zones_copy.index]
        weights=np.divide(attraction,np.power(impedance, beta))
#         weights=np.array(attraction)
        return np.random.choice(
            zones_copy.index, replace=True, p=weights/sum(weights), size=n)
        
    def simulate(self, simpop_df):
        print('Schedules and locations')
        simpop_df=self.sim.create_simple_HWH_schedules(simpop_df)
        print('OD')
        simpop_df_w_nodes=simpop_df.merge(self.sim.zones['possible_nodes_drive'], left_on='home_geoid', right_index=True, how='left').rename(columns={'possible_nodes_drive': 'home_nodes'})
        simpop_df_w_nodes=simpop_df_w_nodes.merge(self.sim.zones['possible_nodes_drive'], left_on='work_geoid', right_index=True, how='left').rename(columns={'possible_nodes_drive': 'work_nodes'})
        simpop_df_w_nodes['route_distance']=simpop_df_w_nodes.apply(lambda row: self.route_lengths[row['home_nodes'][0]]
                                                                    [row['work_nodes'][0]], axis=1)
        simpop_df_w_mode=self.sim.mode_chooser.predict_modes(simpop_df_w_nodes)
        od_output=simpop_df_w_mode[self.sim.person_attributes+['home_geoid', 'work_geoid']].to_dict(orient='records')
        print('Trip table')
        all_trips_df=self.sim.create_trip_table(simpop_df_w_mode)
        # all_trips_df['route_distance']=all_trips_df.apply(lambda row: self.route_lengths[row['from_possible_nodes_drive'][0]]
        #                                                             [row['to_possible_nodes_drive'][0]], axis=1)
        # all_trips_df=self.sim.mode_chooser(all_trips_df)
        print('Route table')
        route_table=self.sim.get_routes_table(all_trips_df)
        print('DeckGL')
        deckgl_trips=self.sim.routes_to_deckgl_trip(route_table)
        by_mode=route_table.groupby('mode').size()
        active=0
        for mode in ['walk', 'cycle']:
            if mode in by_mode:
                active+=by_mode[mode]
        ind=active/by_mode.sum()
        return od_output, deckgl_trips, ind
    
    def return_indicator(self, geogrid_data):
        print('Starting MM Update')
        start_calc=datetime.datetime.now()
        new_simpop=self.geogrid_updates(geogrid_data)
        new_simpop_df=pd.DataFrame(new_simpop)
        combined_simpop=self.base_simpop_df.append(new_simpop_df)
        sample_simpop_df=combined_simpop.sample(min(self.N_max, len(combined_simpop)))
        od_output, deckgl_trips, ind=self.simulate(sample_simpop_df)
        end_calc=datetime.datetime.now()
        self.post_trips(deckgl_trips)
        self.post_od(od_output)
        end_post=datetime.datetime.now()
        print('Finished MM Update')
        print('\t calculation took {}'.format(end_calc-start_calc))
        print('\t trips post took {}'.format(end_post-end_calc))
        return {'name': 'Active Mobility',
                           'raw_value': ind,
                           'value': ind,
                           'ref_value': 0.1,
                           'viz_type':'radar'}
    
    def post_trips(self, deckgl_trips):
        post_url='https://cityscope-api.smartaraucania.org/api/table/'+self.table_name
        r = requests.post(post_url+'/ABM2', data = json.dumps(deckgl_trips),
            headers={'Content-Type': 'application/json'})
        print('Post ABM: {}'.format(r))

    def post_od(self, od_output):
        post_url='https://cityscope-api.smartaraucania.org/api/table/'+self.table_name
        r = requests.post(post_url+'/od', data = json.dumps(od_output),
            headers={'Content-Type': 'application/json'})
        print('Post OD: {}'.format(r))
        
    def get_combined_zones(self):
        comb_zones=self.sim.zones
        comb_zones.index.name='GEOID'
        return self.sim.zones
    
    def get_edges_geojson(self, network_name):
        G_gdf=osmnx.utils_graph.graph_to_gdfs(self.sim.mob_sys.networks[network_name].G, nodes=False, edges=True)
        return G_gdf[['geometry', 'speed_kph', 'from', 'to']]


if __name__ == "__main__":
    table_name=sys.argv[1]

    # Load the saved data
    geogrid=gpd.read_file('tables/{}/geogrid.geojson'.format(table_name))
    zones=gpd.read_file('tables/{}/zones.geojson'.format(table_name))
    zones['GEOID']=zones['GEOID'].astype(int)
    zones=zones.set_index('GEOID')
    simpop_df=pd.read_csv('tables/{}/simpop_df.csv'.format(table_name))


    H=Handler(table_name=table_name)
    H.reset_geogrid_data()

    d=Density_Indicator(zones=zones)
    p=Proximity_Indicator(zones=zones, geogrid=geogrid)
    m=Mobility_indicator(zones, geogrid, table_name, simpop_df)
    H.add_indicators([d, p, m])

    H.listen()