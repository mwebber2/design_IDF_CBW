# -*- coding: utf-8 -*-
"""
Created on Mon May 20 15:00:17 2024

@author: webbe

IDF curves 2024 for Chesapeake Bay watershed counties (applying change factors, not depths)
using the data as directly downloaded from the tool: https://midatlantic-idf.rcc-acis.org/

"""

import pandas as pd
import numpy as np
#import math
import statistics as st
import matplotlib.pyplot as plt
#import os
import requests 
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
#import joypy
import geopandas as gpd

#to run this code, need a folder (path) with the following files:
#- counties GIS folder: cb_2018_us_county_20m
#- counties_latlon.csv (as created by author using ArcGIS)
#- state_fips.csv (as downloaded from https://www2.census.gov/geo/docs/reference/codes2020/national_state2020.txt saved as excel --> csv)
#- file from box folder representing stations used to extract change factors (rp_best_depth_index_info.csv)
#- 4 extracted files downloaded from CBP website 
#   - (CBP_data_rcp45_2020-2070, CBP_data_rcp45_2050-2100, CBP_data_rcp85_2020-2070, CBP_data_rcp85_2050-2100)

#%% methods

#method that saves change factors for all return periods in a single dataframe
def mediandf(cf, r, m_list, RP_list):
    X_avg = pd.DataFrame()
    X_avg.insert(0, cf.columns[0], cf[cf.columns[0]]) #first column = name
    X_avg.insert(1, cf.columns[1], cf[cf.columns[1]]) #second column = type 
    X_avg.insert(2, cf.columns[2], cf[cf.columns[2]]) #third column = state
    for i in range(len(m_list)):
        X_avg.insert(3+i, RP_list[i], m_list[i])
    return X_avg

#method that can webscrape for multiple lat-lon locations
def webscrape_atlas14(path, grid, lat, lon):
    filename = grid +"_atlas14.csv"
    linka = r"https://hdsc.nws.noaa.gov/cgi-bin/hdsc/new/fe_text.csv?lat="
    linkb = r"&lon="
    linkc = r"&data=depth&units=metric&series=pds" #metric!
    link = linka + str(lat) + linkb + str(lon) + linkc
    r = requests.get(link)
    with open(path + filename, 'w') as out:
        out.write(r.text)

#%% paths and other user defined variables

path = r"/Users/webbe/Box/Marissa's Research/IDF Curves/IDFcurve_code/atlasmedian_percentiles/"
path_main = r"/Users/webbe/Box/Marissa's Research/IDF Curves/IDFcurve_code/"
path_to_save = path + r"results/"
path_to_save_atlas = path_main + r"atlas14/"

RCP = ["rcp45_2020-2070", "rcp45_2050-2100", "rcp85_2020-2070", "rcp85_2050-2100"]
RP_list = ["2-yr", "5-yr", "10-yr", "25-yr", "50-yr", "100-yr"]
ptiles = ["mean", "min", "10th ptile", "25th ptile", "50th ptile", "75th ptile", "90th ptile", "max"]

#%%
#read in the change factors and save the medians in dataframe organized by RCP
#(only need to do this once and save the csvs)
for p in range(len(ptiles)):
    for r in range(len(RCP)):
        median = [] #each item in this list is median change factor for return period
        cf = pd.read_csv(path_main + "CBP_data_" + RCP[r] +".csv")
        for i in range(len(cf.columns)):
            if cf.columns[i][-10:] == ptiles[p]: #choose the percentile
                median.append(cf[cf.columns[i]])
        output = mediandf(cf, r, median, RP_list)
        output.to_csv(path_to_save + ptiles[p] + "_" + RCP[r] + ".csv")

#mean
for r in range(len(RCP)):
    median = [] #each item in this list is median change factor for return period
    cf = pd.read_csv(path_main + "CBP_data_" + RCP[r] +".csv")
    for i in range(len(cf.columns)):
        if cf.columns[i][-4:] == ptiles[0]: #median
            median.append(cf[cf.columns[i]])
    output = mediandf(cf, r, median, RP_list)
    output.to_csv(path_to_save + "mean_" + RCP[r] + ".csv")
    
#min
for r in range(len(RCP)):
    minimum = [] #each item in this list is median change factor for return period
    cf = pd.read_csv(path_main + "CBP_data_" + RCP[r] +".csv")
    for i in range(len(cf.columns)):
        if cf.columns[i][-3:] == ptiles[1]: #min
            minimum.append(cf[cf.columns[i]])
    output = mediandf(cf, r, minimum, RP_list)
    output.to_csv(path_to_save + "min_" + RCP[r] + ".csv")
    
#max
for r in range(len(RCP)):
    maximum = [] #each item in this list is median change factor for return period
    cf = pd.read_csv(path_main + "CBP_data_" + RCP[r] +".csv")
    for i in range(len(cf.columns)):
        if cf.columns[i][-3:] == ptiles[-1]: #max
            maximum.append(cf[cf.columns[i]])
    output = mediandf(cf, r, maximum, RP_list)
    output.to_csv(path_to_save + "max_" + RCP[r] + ".csv")

#%% 
#get centroids of each county in chesapeake bay watershed dataset
#doing this every time now because I want to save the duplicated list

#get the center of each county in this dataset

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
#county_centroids.to_csv(path + "CBP_countycentroids.csv")

#%%
#read in the change factors, apply the (correct) change factor to Atlas14 for each county
#clean up this section by adding more methods for repeated code

#collect change factors for all percentiles (use my ptiles list)
cf_masterlist = [] #list of lists
for p in range(len(ptiles)):
    cf_indivlist = []
    for r in range(len(RCP)):
        cf_raw = pd.read_csv(path_to_save + ptiles[p] + "_" + RCP[r] + ".csv")
        cf_raw2 = cf_raw.copy()
        for i in range(len(county_duplicated)):
            if county_duplicated[i] == True:
                cf_raw2 = cf_raw2.drop([cf_raw.index[i]])
        cf_indivlist.append(cf_raw2.reset_index()) #each list has 4 dataframes, one for each RCP
    cf_masterlist.append(cf_indivlist) #6 lists of lists, one for each percentile

#lists to save results for atlas 14 protection RQ for all counties
output_45_short = [] #output depths for all counties
output_45_long = []
output_85_short = []
output_85_long = []
output_atlas_CBW = []
output_atlas_min = []
output_atlas_max = []
output_45_short_min = [] #output depths for min and max of RCP 4.5
output_45_long_min = []
output_45_short_max = []
output_45_long_max = []
output_45_short_25 = [] #output depths for 25th ptile of RCP 4.5
output_45_long_25 = []
output_45_short_75 = [] #output depths for 75th ptile of RCP 4.5
output_45_long_75 = []
output_85_short_min = [] #output depths for min and max of RCP 4.5
output_85_long_min = []
output_85_short_max = []
output_85_long_max = []
output_85_short_25 = [] #output depths for 25th ptile of RCP 4.5
output_85_long_25 = []
output_85_short_75 = [] #output depths for 75th ptile of RCP 4.5
output_85_long_75 = []
output_45_shortcf = [] #output change factors for all counties
output_45_longcf = []
output_85_shortcf = []
output_85_longcf = []
output_diff45_short = [] #output the factional difference (y-axis of box plots for all counties
output_diff85_short = []
output_diff45_long = []
output_diff85_long = []
protect_min_master = [] #does atlas14 upper protect from at least one scenario?
protect_max_master = [] #does atlas14 upper protect from max scenario?
perc_diff_min_master_all = [] #for all counties = highest fractional difference between atlas 14 and scenario
perc_diff_max_master_all = [] #for all counties = fractional difference between atlas 14 and max scenario
perc_diff_45short_master = [] #fractional difference between atlas14 upper and RCPs - for all counties
perc_diff_45long_master = []
perc_diff_85short_master = []
perc_diff_85long_master = []
protect_45short_master = [] #does atlas14 upper protect from RCPs? for all counties
protect_45long_master = []
protect_85short_master = []
protect_85long_master = []
output_45short_master = [] #depths
output_45long_master = [] #depths
output_85short_master = [] #depths
output_85long_master = [] #depths

#for many locations 
#add some filter here = only proceed if POR_list[c] == 0 #no longer needed, we flag them later
for c in range(len(county_centroids)):
    grid_test = county_centroids['NAME'][c] + "_" + county_centroids['STATE'][c]
    try:
        atlas = pd.read_csv(path_to_save_atlas + grid_test + "_atlas14.csv", skiprows=13)
    except:
        lat_test = county_centroids['lat'][c]
        lon_test = county_centroids['lon'][c]
        webscrape_atlas14(path_to_save_atlas, grid_test, lat_test, lon_test)
        atlas = pd.read_csv(path_to_save_atlas + grid_test + "_atlas14.csv", skiprows=13)
    #read the correct line that I want from this csv
    line = []
    for i in range(len(atlas)):
        if atlas['by duration for ARI (years):'][i] == "24-hr:":
            line.append(i)
    atlas_median = atlas.iloc[line[0]][2:8]
    atlas_upper = atlas.iloc[line[1]][2:8]
    atlas_lower = atlas.iloc[line[2]][2:8]
    #find fractional difference between atlas_median and atlas_upper
    atlas_diff = np.divide(np.array(atlas_upper), np.array(atlas_median))
    output = []
    output_min = []
    output_max = []
    output_25 = []
    output_75 = []
    outputcf = []
    outputdiff = []
    
    def cf_ptile(p, j, i): #method creates the output depths for the 6 return periods
        #for percentile p, RCP j, and county i
        cf = cf_masterlist[p][j].iloc[i][5:]
        out = np.multiply(np.array(atlas_median), np.array(cf))
        return out
    
    for j in range(len(cf_masterlist[0])): #all 4 RCPs
        for i in range(len(cf_masterlist[0][j])): #321 counties
            if grid_test == cf_masterlist[0][j]['NAME'][i] + "_" + cf_masterlist[0][j]['STATE'][i]:
                cf_2 = cf_masterlist[4][j].iloc[i][5:] #median
                output_min.append(cf_ptile(1, j, i))#output_min.append(mini)
                output_25.append(cf_ptile(3, j, i)) #25th
                output.append(cf_ptile(4, j, i))#output.append(med)
                output_75.append(cf_ptile(5, j, i)) #75th
                output_max.append(cf_ptile(-1, j, i))#output_max.append(maxi)
                outputcf.append(cf_2)
                outputdiff.append(atlas_diff)
    output_atlas_CBW.append(atlas_median)
    output_atlas_min.append(atlas_lower)
    output_atlas_max.append(atlas_upper)
    output_45_short_min.append(output_min[0]) #min
    output_45_long_min.append(output_min[1])
    output_45_short_25.append(output_25[0]) #25th ptile
    output_45_long_25.append(output_25[1])
    output_45_short_75.append(output_75[0]) #75th ptile
    output_45_long_75.append(output_75[1])
    output_45_short_max.append(output_max[0]) #max
    output_45_long_max.append(output_max[1])
    output_45_short.append(output[0])
    output_45_long.append(output[1])
    output_85_short.append(output[2])
    output_85_long.append(output[3])
    output_85_short_min.append(output_min[2]) #min
    output_85_long_min.append(output_min[3])
    output_85_short_25.append(output_25[2]) #25th ptile
    output_85_long_25.append(output_25[3])
    output_85_short_75.append(output_75[2]) #75th ptile
    output_85_long_75.append(output_75[3])
    output_85_short_max.append(output_max[2]) #max
    output_85_long_max.append(output_max[3])
    output_45_shortcf.append(outputcf[0])
    output_45_longcf.append(outputcf[1])
    output_85_shortcf.append(outputcf[2])
    output_85_longcf.append(outputcf[3])
    output_diff45_short.append(outputdiff[0])
    output_diff45_long.append(outputdiff[1])
    output_diff85_short.append(outputdiff[2])
    output_diff85_long.append(outputdiff[3])
    
    #way to keep track of which is higher (atlas 14 vs projections) or how much higher they are
    #depends on what I want to plot/show
    r_master = []
    protect_min = []
    protect_max = []
    perc_diff_min_all = []
    perc_diff_max_all = []
    perc_diff_45short = []
    perc_diff_45long = []
    perc_diff_85short = []
    perc_diff_85long = []
    protect_45short = []
    protect_45long = []
    protect_85short = []
    protect_85long = []
    
    for a in range(len(atlas_median)): #for each RP
        r_list = []
        for o in range(len(output)): #save the depths for each scenario
            r_list.append(output[o][a])
        r_master.append(r_list) #6 lists, each list = 4 values (one for each scenario)
    
    for r in range(len(r_master)): #for each RP #TODO update this section for other strategies
        diff_min = (atlas_median[r] - min(r_master[r]))/atlas_median[r]
        perc_diff_min_all.append(diff_min)
        if diff_min >= 0: #check if atlas 14 >= any of the depths
            protect_min.append(1) #protected
        elif diff_min < 0: 
            protect_min.append(0)
        
        diff_max = (atlas_median[r] - max(r_master[r]))/atlas_median[r]
        perc_diff_max_all.append(diff_max)
        if diff_max >= 0: #check if atlas 14 >= max depths
            protect_max.append(1) #protected
        elif diff_max < 0: 
            protect_max.append(0)
            
        diff_45short = (atlas_median[r] - r_master[r][0])/atlas_median[r]
        perc_diff_45short.append(diff_45short)
        if ((atlas_median[r] - output_max[0][r])/atlas_median[r]) >= 0:
            pp1 = 1
        elif ((atlas_median[r] - output_75[0][r])/atlas_median[r]) >= 0:
            pp1 = 0.75
        elif diff_45short >= 0:
            pp1 = 0.5
        elif ((atlas_median[r] - output_25[0][r])/atlas_median[r]) >= 0:
            pp1 = 0.25
        #else: pp1 = 0
        elif ((atlas_median[r] - output_min[0][r])/atlas_median[r]) >= 0:
            pp1 = 0
        else: pp1 = np.nan #-1
        protect_45short.append(pp1)
        
        diff_45long = (atlas_median[r] - r_master[r][1])/atlas_median[r]
        perc_diff_45long.append(diff_45long)
        if ((atlas_median[r] - output_max[1][r])/atlas_median[r]) >= 0:
            pp2 = 1
        elif ((atlas_median[r] - output_75[1][r])/atlas_median[r]) >= 0:
            pp2 = 0.75
        elif diff_45long >= 0:
            pp2 = 0.5
        elif ((atlas_median[r] - output_25[1][r])/atlas_median[r]) >= 0:
            pp2 = 0.25
        #else: pp2 = 0
        elif ((atlas_median[r] - output_min[1][r])/atlas_median[r]) >= 0:
            pp2 = 0
        else: pp2 = np.nan #-1
        protect_45long.append(pp2)

        diff_85short = (atlas_median[r] - r_master[r][2])/atlas_median[r]
        perc_diff_85short.append(diff_85short)
        if ((atlas_median[r] - output_max[2][r])/atlas_median[r]) >= 0:
            pp3 = 1
        elif ((atlas_median[r] - output_75[2][r])/atlas_median[r]) >= 0:
            pp3 = 0.75
        elif diff_85short >= 0:
            pp3 = 0.5
        elif ((atlas_median[r] - output_25[2][r])/atlas_median[r]) >= 0:
            pp3 = 0.25
        #else: pp2 = 0
        elif ((atlas_median[r] - output_min[2][r])/atlas_median[r]) >= 0:
            pp3 = 0
        else: pp3 = np.nan #-1
        protect_85short.append(pp3)
        
        diff_85long = (atlas_median[r] - r_master[r][3])/atlas_median[r]
        perc_diff_85long.append(diff_85long)
        if ((atlas_median[r] - output_max[3][r])/atlas_median[r]) >= 0:
            pp4 = 1
        elif ((atlas_median[r] - output_75[3][r])/atlas_median[r]) >= 0:
            pp4 = 0.75
        elif diff_85long >= 0:
            pp4 = 0.5
        elif ((atlas_median[r] - output_25[3][r])/atlas_median[r]) >= 0:
            pp4 = 0.25
        #else: pp2 = 0
        elif ((atlas_median[r] - output_min[3][r])/atlas_median[r]) >= 0:
            pp4 = 0
        else: pp4 = np.nan #-1
        protect_85long.append(pp4)

    protect_min_master.append(protect_min)
    protect_max_master.append(protect_max)
    perc_diff_min_master_all.append(perc_diff_min_all)
    perc_diff_max_master_all.append(perc_diff_max_all)
    
    protect_45short_master.append(protect_45short)
    protect_45long_master.append(protect_45long)
    protect_85short_master.append(protect_85short)
    protect_85long_master.append(protect_85long)
    perc_diff_45short_master.append(perc_diff_45short)
    perc_diff_45long_master.append(perc_diff_45long)
    perc_diff_85short_master.append(perc_diff_85short)
    perc_diff_85long_master.append(perc_diff_85long)

#%%
#method rearranges a list of lists for counties[rps] i.e. 321 lists with each list 6 items
#to a list of lists for rps[counties] i.e. 6 lists with each list 321 items
def resultsbycountyandRP(listoflists):
    diff_allRP = []
    for rp in range(len(RP_list)):
        diff = []
        for p in range(len(listoflists)):
            #if len(listoflists[p])>rp: #no longer need a conditional because I have results for all return periods
            diff.append(listoflists[p][rp])
        diff_allRP.append(diff)
    return diff_allRP

#for all counties
diff_allRP_protectmin = resultsbycountyandRP(protect_min_master)
diff_allRP_protectmax = resultsbycountyandRP(protect_max_master)
diff_allRP_percmin_all = resultsbycountyandRP(perc_diff_min_master_all)
diff_allRP_percmax_all = resultsbycountyandRP(perc_diff_max_master_all)

#for all counties for different different RCPs
diff_allRP_percmin_RCP45_short = resultsbycountyandRP(perc_diff_45short_master)
diff_allRP_percmin_RCP85_short = resultsbycountyandRP(perc_diff_85short_master)
diff_allRP_percmin_RCP45_long = resultsbycountyandRP(perc_diff_45long_master)
diff_allRP_percmin_RCP85_long = resultsbycountyandRP(perc_diff_85long_master)
diff_allRP_protect_45short = resultsbycountyandRP(protect_45short_master)
diff_allRP_protect_45long = resultsbycountyandRP(protect_45long_master)
diff_allRP_protect_85short = resultsbycountyandRP(protect_85short_master)
diff_allRP_protect_85long = resultsbycountyandRP(protect_85long_master)


#%%
allUScounties = gpd.read_file(path_main + "cb_2018_us_county_20m/cb_2018_us_county_20m.shp")
#select counties with GEOID from county_centroids
geoid_list = []
for i in range(len(county_centroids)):
    geoid_list.append(str(county_centroids['GEOID'][i]))
CBWcounties = allUScounties[allUScounties['GEOID'].isin(geoid_list)]

#method to plot maps for each RP
#new color bar
N = 256
vals2 = np.ones((N, 4))
vals2[:, 0] = np.linspace(173/256, 0/256, N)
vals2[:, 1] = np.linspace(216/256, 0/256, N)
vals2[:, 2] = np.linspace(230/256, 255/256, N)
newblue = ListedColormap(vals2)

#col_dict = {1: "azure", #"silver",
#            2: "cornflowerblue",
#            3: "blue",
#            4: "midnightblue"}
#col_dict = {1: "orangered", #stoplight
#            2: "mistyrose", #"lightcoral",
#            3: "gold",
#            4: "darkgreen"}
col_dict = {1: "tomato", #"red", #red to blue
            2: "mistyrose", #"lightcoral",
            3: "lightskyblue", #"cornflowerblue",
            4: "blue"} #kinda difficult to see in monochrome #1 vs #4
col_list = ["tomato", #"red", #red to blue
            "mistyrose", #"lightcoral",
            "lightskyblue", #"cornflowerblue",
            "blue"]
col_listwgrey = ["firebrick",#"white", #"lightgrey",
                 "tomato", #"red", #red to blue
                 "mistyrose", #"lightcoral",
                 "lightskyblue", #"cornflowerblue",
                 "blue"]
newblue = ListedColormap([col_dict[x] for x in col_dict.keys()])
newblue_list = []
for i in range(len(col_list)):
    newblue_indiv = ListedColormap([col_list[x] for x in range(i+1)])
    newblue_list.append(newblue_indiv)
newbluewgrey_list = []
for i in range(len(col_listwgrey)):
    newbluewgrey_indiv = ListedColormap([col_listwgrey[x] for x in range(i+1)])
    newbluewgrey_list.append(newbluewgrey_indiv)

#method for plotting multiple maps at the same time
def multiplemaps(listname, titletosave): #with shared legend
    fig, axs = plt.subplots(2, 3, figsize = (10, 8.8))
    axs = axs.flatten()
    #vmin, vmax = 0.0, 1.0
    for i in range(len(RP_list)): #6 RPs
        CBWcounties[RP_list[i]] = listname[i]
        variable = RP_list[i]
        set_len = len(set(CBWcounties[RP_list[i]].dropna())) #if using NAs
        cmap_curr = newblue_list[set_len-1]
        #set_len = len(set(CBWcounties[RP_list[i]])) #if using -1
        #print(set(CBWcounties[RP_list[i]]))
        if min(set(CBWcounties[RP_list[i]].dropna())) == 0: 
            cmap_curr = newbluewgrey_list[set_len-1]
        if min(set(CBWcounties[RP_list[i]].dropna())) == 0.25: 
            cmap_curr = newblue_list[set_len-1]
        #CBWcounties.plot(column = variable, cmap=plt.cm.get_cmap(newblue), vmin = vmin, vmax = vmax, linewidth=0.8, ax=axs[i], edgecolor='0.8', categorical=False)#, legend = True)
        CBWcounties.plot(facecolor="white", hatch="XX", ax = axs[i])#, edgecolor = '0.3')
        if i == (len(RP_list)-1): #or i == 2:
            CBWcounties.plot(column = variable, cmap=plt.cm.get_cmap(cmap_curr), linewidth=0.8, ax=axs[i], edgecolor='0.8', categorical=True,
                             legend = True, legend_kwds={'bbox_to_anchor':(2, 1.05), 'fontsize':12, 
                                                         'title': "% of counties \n climate ready \n at different percentiles",
                                                         'title_fontsize':12})
        else:
            CBWcounties.plot(column = variable, cmap=plt.cm.get_cmap(cmap_curr), linewidth=0.8, ax=axs[i], edgecolor='0.8', categorical=True, 
                             legend = False)
        axs[i].set_title(RP_list[i])
        axs[i].set_xticks([])
        axs[i].set_yticks([])
        # add text to the figures
        t1 = str(int((sum(num >= 1 for num in list(CBWcounties[RP_list[i]]))/321)*100))
        t75 = str(int((sum(num >= 0.75 for num in list(CBWcounties[RP_list[i]]))/321)*100))
        t50 = str(int((sum(num >= 0.5 for num in list(CBWcounties[RP_list[i]]))/321)*100))
        t25 = str(int((sum(num >= 0.25 for num in list(CBWcounties[RP_list[i]]))/321)*100))
        t0 = str(int((sum(num >= 0 for num in list(CBWcounties[RP_list[i]]))/321)*100))
        #axs[i].text(0.05, 0.65, t1 + "% max \n" + t75 + "% 75th \n" + t50 + "% 50th \n" + t25 + "% 25th \n" + t0 + "% min", 
        #            transform=axs[i].transAxes)
        #axs[i].text(0.05, 0.75, t50 + "% 50th \n" + t25 + "% 25th \n" + t0 + "% min", 
        #            transform=axs[i].transAxes) #when the fig size was (10, 7)
        #axs[i].text(-60, 42, t50 + "% 50th \n" + t25 + "% 25th \n" + t0 + "% min")
        #axs[i].text(0.25, 1, t50 + "% 50th \n" + t25 + "% 25th \n" + t0 + "% min",
        #            horizontalalignment='left',
        #            verticalalignment='top',
        #            transform = ax.transAxes)
        axs[i].table([['percentile', '% counties'],['maximum', t1+"%"], ['75th', t75+"%"], ['50th', t50+"%"], ['25th', t25+"%"], ['minimum', t0+"%"]])
    #cbar_ax = fig.add_axes([0.125, 0.075, 0.78, 0.04]) #specify location of the colorbar = bottom
    #sm = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap(newblue, 4), norm=plt.Normalize(vmin=vmin, vmax=vmax))
    #sm = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap('tab20b', 4), norm=plt.Normalize(vmin=vmin, vmax=vmax))
    #sm._A = [] #empty array for the data range
    #cbar = fig.colorbar(sm, cax = cbar_ax, orientation = "horizontal", ticks = [0, 0.25, 0.5, 0.75, 1]) #add the colorbar to the figure
    plt.savefig(path + titletosave, bbox_inches='tight')
    #plt.close()

#multiplemaps(diff_allRP_protectmin, "_minprotectionmap_allRPs.png") #counties protected from min RCP and time period
#multiplemaps(diff_allRP_protectmax, "_maxprotectionmap_allRPs.png") #counties protected from max RCP and time period

#counties protected from each RP
multiplemaps(diff_allRP_protect_45short, "percentilemap_RCP45_2020-2070.png") #RCP 4.5 2020-2070
multiplemaps(diff_allRP_protect_45long, "percentilemap_RCP45_2050-2100.png") #RCP 4.5 2050-2100
multiplemaps(diff_allRP_protect_85short, "percentilenmap_RCP85_2020-2070.png") #RCP 8.5 2020-2070
multiplemaps(diff_allRP_protect_85long, "percentilenmap_RCP85_2050-2100.png") #RCP 8.5 2050-2100
