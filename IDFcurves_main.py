# -*- coding: utf-8 -*-
"""
Created on Sun May 19 08:52:22 2024

@author: webbe

This file imports .py files for IDF curves


"""

#Dependencies
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
import geopandas as gpd
pd.options.mode.chained_assignment = None  # default='warn'

t0 = time.perf_counter() #seconds

#%%
#user-defined paths
path_main = r"/Users/webbe/Box/Marissa's Research/IDF Curves/IDFcurve_code/"
path_to_save_atlas = path_main + r"atlas14/"

#pre-defined lists
RCP = ["rcp45_2020-2070", "rcp45_2050-2100", "rcp85_2020-2070", "rcp85_2050-2100"]
RP_list = ["2-yr", "5-yr", "10-yr", "25-yr", "50-yr", "100-yr"]
ptiles = ["mean", "min", "10th ptile", "25th ptile", "50th ptile", "75th ptile", "90th ptile", "max"]


#%% county centroids
#get centroids of each county in chesapeake bay watershed dataset
#doing this every time now because I want to save the duplicated list

#start by reading in all the counties for which we have data from the tool
cf = pd.read_csv(path_main + "CBP_data_" + RCP[0] +".csv")
county_centroids = pd.DataFrame()
county_centroids.insert(0, cf.columns[0], cf[cf.columns[0]]) #first column = name
county_centroids.insert(1, cf.columns[1], cf[cf.columns[1]]) #second column = type 
county_centroids.insert(2, cf.columns[2], cf[cf.columns[2]]) #third column = state

#find centroid of each county in the US in arcgis and import the file here
latlon_US = pd.read_csv(path_main + "counties_latlon.csv")
stateIDs = pd.read_csv(path_main + "state_fips.csv")
#countyfips = pd.read_csv(path + "county_fips_master.csv")

#save the state fips code for each state in our county_centroids
state_fips = []
for i in range(len(county_centroids)):
    for j in range(len(stateIDs)):
        if county_centroids["STATE"].iloc[i] == stateIDs['STATE'].iloc[j]:
            state_fips.append(stateIDs['STATEFP'].iloc[j])
county_centroids.insert(3, "stateIDs", state_fips)

#using the county name and state fips, save lat lon and GEOID for each county
lat_list = []
lon_list = []
countyID = []
listofindices = []
for i in range(len(county_centroids)):
    for j in range(len(latlon_US)):
        if county_centroids["NAME"].iloc[i] + "_" + str(county_centroids["stateIDs"].iloc[i]) == latlon_US["NAME"].iloc[j] + "_" + str(latlon_US["STATEFP"].iloc[j]):
            listofindices.append(i)
            lat_list.append(latlon_US["Latitude"][j])
            lon_list.append(latlon_US["Longitude"][j])
            countyID.append(latlon_US['GEOID'][j])
#use difference between indices to remove extra lat lons 
#so that this data matches length of county centroids from the tool (N = 326)
diff = []
for i in range(len(listofindices)-1):
    diff.append(listofindices[i+1] - listofindices[i])
count = 0
droplist = []
for d in range(len(diff)): 
    if diff[d] == 0: 
        count+=1
        droplist.append(d)
#anywhere diff = 0, drop that lat lon (it means I have a repeat or an extra)
indicesList = sorted(droplist, reverse=True)
for indx in indicesList: # Traversing in the indices list
   if indx < len(lat_list): # checking whether the corresponding iterator index is less than the list length
      lat_list.pop(indx) # removing element by index using pop() function
for indx in indicesList: # Traversing in the indices list
   if indx < len(lon_list): # checking whether the corresponding iterator index is less than the list length
      lon_list.pop(indx) # removing element by index using pop() function
for indx in indicesList: # Traversing in the indices list
   if indx < len(countyID): # checking whether the corresponding iterator index is less than the list length
      countyID.pop(indx) # removing element by index using pop() function

#save to this dataframe
county_centroids.insert(4, "lat", lat_list)
county_centroids.insert(5, "lon", lon_list)
county_centroids.insert(6, "GEOID", countyID)

#remove any duplicates based on GEOIDs
county_duplicated = county_centroids.duplicated(subset=['GEOID'])
county_centroids = county_centroids.drop_duplicates(subset=["GEOID"]).reset_index()

#save CSV
county_centroids.to_csv(path_main + "CBP_countycentroids.csv")

#%% simple map of study area
#use geopandas to plotmap of CBW counties
allUScounties = gpd.read_file(path_main + "cb_2018_us_county_20m/cb_2018_us_county_20m.shp")
#select counties with GEOID from county_centroids
geoid_list = []
for i in range(len(county_centroids)):
    geoid_list.append(str(county_centroids['GEOID'][i]))
CBWcounties = allUScounties[allUScounties['GEOID'].isin(geoid_list)]
#fig, ax = plt.subplots(1, figsize=(10, 6))
CBWcounties.plot()
plt.savefig(path_main + "CBWmap.png")
plt.close()

#%% main files

#import IDFcurves_checkPOR #TODO does not generate any figures?
import IDFcurves_atlasmedian #strategy 1
import IDFcurves_atlasCI #strategy 2
import IDFcurves_atlasupRP #strategy 3
import IDFcurves_atlasupRPCI #strategy 4
import IDFcurves_maps #maps and figures comparing four strategies

#import IDFcurves_atlasmedian_percentiles #SI figures for strategy 1
#import IDFcurves_atlasCI_percentiles #SI figures for strategy 2

t1 = time.perf_counter() #seconds

totaltime_sec = t1-t0 #3435.7 #2484.3
totaltime_min = totaltime_sec/60 #57.26 #41.41
totaltime_hr = totaltime_sec/3600 #0.954 #0.69


