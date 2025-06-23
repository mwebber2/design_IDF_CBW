# -*- coding: utf-8 -*-
"""
Created on Sun May 19 10:31:09 2024

@author: webbe

IDF curves 2024 for Chesapeake Bay watershed counties (applying change factors, not depths)
using the data as directly downloaded from the tool: https://midatlantic-idf.rcc-acis.org/

"""

import pandas as pd
import numpy as np
#import math
#import statistics as st
import matplotlib.pyplot as plt
import os
import requests 
#import joypy
#import geopandas as gpd

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

path_main = r"/Users/webbe/Box/Marissa's Research/IDF Curves/IDFcurve_code/"
path = path_main + r"atlasCI/"
path_to_save = path + r"results/"
path_to_save_atlas = path_main + r"atlas14/"


if not os.path.exists(path):
    os.makedirs(path)
if not os.path.exists(path_to_save):
    os.makedirs(path_to_save)
if not os.path.exists(path_to_save_atlas):
    os.makedirs(path_to_save_atlas)

RCP = ["rcp45_2020-2070", "rcp45_2050-2100", "rcp85_2020-2070", "rcp85_2050-2100"]
RP_list = ["2-yr", "5-yr", "10-yr", "25-yr", "50-yr", "100-yr"]
ptiles = ["mean", "min", "10th ptile", "25th ptile", "50th ptile", "75th ptile", "90th ptile", "max"]

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

cf_list = []
cf_list_min = []
cf_list_max = []
for r in range(len(RCP)):
    cf_raw = pd.read_csv(path_main + "median_" + RCP[r] + ".csv")
    cf_raw_min = pd.read_csv(path_main + "min_" + RCP[r] + ".csv")
    cf_raw_max = pd.read_csv(path_main + "max_" + RCP[r] + ".csv")
    
    #drop the counties with repeated GEOIDs
    cf_raw2 = cf_raw.copy()
    cf_raw_min2 = cf_raw_min.copy()
    cf_raw_max2 = cf_raw_max.copy()
    for i in range(len(county_duplicated)):
        if county_duplicated[i] == True:
            cf_raw2 = cf_raw2.drop([cf_raw.index[i]])
            cf_raw_min2 = cf_raw_min2.drop([cf_raw_min.index[i]])
            cf_raw_max2 = cf_raw_max2.drop([cf_raw_max.index[i]])
    cf_list.append(cf_raw2.reset_index())
    cf_list_min.append(cf_raw_min2.reset_index())
    cf_list_max.append(cf_raw_max2.reset_index())

#lists to save results for atlas 14 protection RQ for all counties
output_45_short = [] #output depths for all counties
output_45_long = []
output_85_short = []
output_85_long = []
output_atlas_CBW = [] #output depths for all counties for atlas and upper and lower
output_atlas_min = []
output_atlas_max = []
output_45_short_min = [] #output depths for min and max of RCP 4.5
output_45_long_min = []
output_45_short_max = []
output_45_long_max = []
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
    outputcf = []
    outputdiff = []
    for j in range(len(cf_list)):
        for i in range(len(cf_list[j])):
            if grid_test == cf_list[j]['NAME'][i] + "_" + cf_list[j]['STATE'][i]:
                cf_2 = cf_list[j].iloc[i][5:]
                cf_min = cf_list_min[j].iloc[i][5:]
                cf_max = cf_list_max[j].iloc[i][5:]
                #multiply atlas median, upper lower by the change factors
                med = np.multiply(np.array(atlas_median),np.array(cf_2))
                mini = np.multiply(np.array(atlas_median),np.array(cf_min)) 
                maxi = np.multiply(np.array(atlas_median),np.array(cf_max))
                #upper = np.multiply(np.array(atlas_upper),np.array(cf_2))
                #lower = np.multiply(np.array(atlas_lower),np.array(cf_2))
                output.append(med)
                output_min.append(mini)
                output_max.append(maxi)
                outputcf.append(cf_2)
                outputdiff.append(atlas_diff)
    output_atlas_CBW.append(atlas_median)
    output_atlas_min.append(atlas_lower)
    output_atlas_max.append(atlas_upper)
    output_45_short_min.append(output_min[0])
    output_45_long_min.append(output_min[1])
    output_45_short_max.append(output_max[0])
    output_45_long_max.append(output_max[1])
    output_45_short.append(output[0])
    output_45_long.append(output[1])
    output_85_short.append(output[2])
    output_85_long.append(output[3])
    output_45_shortcf.append(outputcf[0])
    output_45_longcf.append(outputcf[1])
    output_85_shortcf.append(outputcf[2])
    output_85_longcf.append(outputcf[3])
    output_diff45_short.append(outputdiff[0])
    output_diff45_long.append(outputdiff[1])
    output_diff85_short.append(outputdiff[2])
    output_diff85_long.append(outputdiff[3])
    
    #way to keep track of variables for individual counties (and later to be aggregated to lists for all counties)
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
    
    for a in range(len(atlas_upper)): #for each RP
        r_list = []
        for o in range(len(output)): #save the atlas14 depths for each scenario
            r_list.append(output[o][a])
        r_master.append(r_list) #6 lists, each list = 4 values (one for each scenario)

    for r in range(len(r_master)): #for each RP
        diff_min = (atlas_upper[r] - min(r_master[r]))/atlas_upper[r]
        perc_diff_min_all.append(diff_min)
        if diff_min >= 0: #check if atlas 14 >= any of the depths
            protect_min.append(1) #protected
        elif diff_min < 0: 
            protect_min.append(0)
        
        diff_max = (atlas_upper[r] - max(r_master[r]))/atlas_upper[r]
        perc_diff_max_all.append(diff_max)
        if diff_max >= 0: #check if atlas 14 >= max depths
            protect_max.append(1) #protected
        elif diff_max < 0: 
            protect_max.append(0)
            
        diff_45short = (atlas_upper[r] - r_master[r][0])/atlas_upper[r]
        perc_diff_45short.append(diff_45short)
        if diff_45short >= 0:
            protect_45short.append(1)
        elif diff_45short < 0:
            protect_45short.append(0)
        
        diff_45long = (atlas_upper[r] - r_master[r][1])/atlas_upper[r]
        perc_diff_45long.append(diff_45long)
        if diff_45long >= 0:
            protect_45long.append(1)
        elif diff_45long < 0:
            protect_45long.append(0)

        diff_85short = (atlas_upper[r] - r_master[r][2])/atlas_upper[r]
        perc_diff_85short.append(diff_85short)
        if diff_85short >= 0:
            protect_85short.append(1)
        elif diff_85short < 0:
            protect_85short.append(0)
        
        diff_85long = (atlas_upper[r] - r_master[r][3])/atlas_upper[r]
        perc_diff_85long.append(diff_85long)
        if diff_85long >= 0:
            protect_85long.append(1)
        elif diff_85long < 0:
            protect_85long.append(0)

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

#export these files for future mapping
#protected 0 vs 1
#pd.DataFrame(diff_allRP_protectmin).T.to_csv(path_to_save + "output_protect_min.csv")
#pd.DataFrame(diff_allRP_protectmax).T.to_csv(path_to_save + "output_protect_max.csv")
#pd.DataFrame(diff_allRP_protect_45short).T.to_csv(path_to_save + "output_protect_45short.csv")
#pd.DataFrame(diff_allRP_protect_45long).T.to_csv(path_to_save + "output_protect_45long.csv")
#pd.DataFrame(diff_allRP_protect_85short).T.to_csv(path_to_save + "output_protect_85short.csv")
#pd.DataFrame(diff_allRP_protect_85long).T.to_csv(path_to_save + "output_protect_85long.csv")
#percentage protection/fracitonal difference
#pd.DataFrame(diff_allRP_percmin_all).T.to_csv(path_to_save + "output_percmin.csv")
#pd.DataFrame(diff_allRP_percmax_all).T.to_csv(path_to_save + "output_prercmax.csv")
pd.DataFrame(diff_allRP_percmin_RCP45_short).T.to_csv(path_to_save + "output_frac_45short.csv")
pd.DataFrame(diff_allRP_percmin_RCP45_long).T.to_csv(path_to_save + "output_frac_45long.csv")
pd.DataFrame(diff_allRP_percmin_RCP85_short).T.to_csv(path_to_save + "output_frac_85short.csv")
pd.DataFrame(diff_allRP_percmin_RCP85_long).T.to_csv(path_to_save + "output_frac_85long.csv")
#depths
pd.DataFrame(output_atlas_CBW).to_csv(path_to_save + "output_atlas_CBW.csv")
pd.DataFrame(output_atlas_min).to_csv(path_to_save + "output_atlas_min.csv")
pd.DataFrame(output_atlas_max).to_csv(path_to_save + "output_atlas_max.csv")
pd.DataFrame(output_45_short).to_csv(path_to_save + "output_45short.csv")
pd.DataFrame(output_45_long).to_csv(path_to_save + "output_45long.csv")
pd.DataFrame(output_85_short).to_csv(path_to_save + "output_85short.csv")
pd.DataFrame(output_85_long).to_csv(path_to_save + "output_85long.csv")

#separating out counties that are only protected (pos) or only unprotected (neg)
diff_allRP_percmin_pos = []
diff_allRP_percmin_neg = []
for i in range(len(diff_allRP_percmin_all)):
    RP_pos = []
    RP_neg = []
    for j in diff_allRP_percmin_all[i]:
        if j >= 0:
            RP_pos.append(j)
        elif j < 0:
            RP_neg.append(j)
    diff_allRP_percmin_pos.append(RP_pos)
    diff_allRP_percmin_neg.append(RP_neg)

diff_allRP_percmax_pos = []
diff_allRP_percmax_neg = []
for i in range(len(diff_allRP_percmax_all)):
    RP_pos = []
    RP_neg = []
    for j in diff_allRP_percmax_all[i]:
        if j >= 0:
            RP_pos.append(j)
        elif j < 0:
            RP_neg.append(j)
    diff_allRP_percmax_pos.append(RP_pos)
    diff_allRP_percmax_neg.append(RP_neg)

#convert all the fractions to factors of safety (i.e. add 1+fraction)
#with if statement that factors cannot be < 1.0
def fractiontoFS(listoflists):
    FS = []    
    for l in listoflists:        
        fs = []
        for i in l:
            if i <= 0: #if fraction < 0, FS = 1
                fs.append(1.0)
            else: #else FS = 1 + fraction
                fs.append(1+i)
        FS.append(fs) #output a new df
    return FS

#convert all fractions to a percentage (i.e. add 1+fraction)
def fractiontoCRF(listoflists):
    FS = []    
    for l in listoflists:        
        fs = []
        for i in l:
            fs.append(1+i)
        FS.append(fs) #output a new df
    return FS

FS_allRP_percmin_all = fractiontoFS(diff_allRP_percmin_all)
FS_allRP_percmax_all = fractiontoFS(diff_allRP_percmax_all)
FS_allRP_percmin_pos = fractiontoFS(diff_allRP_percmin_pos)
FS_allRP_percmax_pos = fractiontoFS(diff_allRP_percmax_pos)
#doesnt make sense to do neg = they will all have FS = 1
FS_allRP_percmin_RCP45_short = fractiontoFS(diff_allRP_percmin_RCP45_short)
FS_allRP_percmin_RCP45_long = fractiontoFS(diff_allRP_percmin_RCP45_long)
FS_allRP_percmin_RCP85_short = fractiontoFS(diff_allRP_percmin_RCP85_short)
FS_allRP_percmin_RCP85_long = fractiontoFS(diff_allRP_percmin_RCP85_long)
#also calculate for CRF < 1
CRF_allRP_percmin_RCP45_short = fractiontoCRF(diff_allRP_percmin_RCP45_short)
CRF_allRP_percmin_RCP45_long = fractiontoCRF(diff_allRP_percmin_RCP45_long)
CRF_allRP_percmin_RCP85_short = fractiontoCRF(diff_allRP_percmin_RCP85_short)
CRF_allRP_percmin_RCP85_long = fractiontoCRF(diff_allRP_percmin_RCP85_long)

#methods for plotting
def boxplots_counties(data, xlab, ylimits, hline, ylabel, titletosave):
    plt.figure()
    plt.boxplot(data)
    plt.xticks(range(1, 7), xlab)#, rotation = 25)
    plt.hlines(hline, 0.5, 6.5, color = "grey", linestyle="dotted")
    plt.ylim(ylimits)
    plt.xlabel("Return Period")
    plt.ylabel(ylabel)
    plt.savefig(path + titletosave, bbox_inches='tight')
    #plt.close()

def violinplots_counties(data, xlab, ylimits, hline, ylabel, titletosave):
    plt.figure()
    plt.violinplot(data)
    plt.xticks(range(1, 7), xlab)#, rotation = 25)
    plt.hlines(hline, 0.5, 6.5, color = "grey", linestyle="dotted")
    plt.ylim(ylimits)
    plt.xlabel("Return Period")
    plt.ylabel(ylabel)
    plt.savefig(path + titletosave, bbox_inches='tight')
    plt.close()
    
#uniform min and max for y axis
yaxislimits = (-0.4, 0.4)
yaxislimits_FS = (0.5, 2.0)

#labels for all counties - representing # of counties protected
xlabels_min = []
xlabels_max = []
for i in range(len(diff_allRP_percmin_all)): #label for protected counties
    #xlabels_min.append(RP_list[i] + " (N = " + str(len(diff_allRP_percmin_pos[i])) + ")")
    #xlabels_max.append(RP_list[i] + " (N = " + str(len(diff_allRP_percmax_pos[i])) + ")")
    xlabels_min.append(RP_list[i] + "\n(" + str(int(len(diff_allRP_percmin_pos[i])/len(county_centroids)*100)) + "%)")
    xlabels_max.append(RP_list[i] + "\n(" + str(int(len(diff_allRP_percmax_pos[i])/len(county_centroids)*100)) + "%)")

#Box plot #all counties
#fraction
#boxplots_counties(diff_allRP_percmin_all, xlabels_min, yaxislimits, 0, "Fractional Difference (-)", "minprotection_allcounties.png")
#boxplots_counties(diff_allRP_percmax_all, xlabels_max, yaxislimits, 0, "Fractional Difference (-)", "maxprotection_allcounties.png")
#FS
#boxplots_counties(FS_allRP_percmin_all, xlabels_min, yaxislimits_FS, 1, "Climate Factor of Safety (-)", "minprotection_allcounties_FS.png")
#boxplots_counties(FS_allRP_percmax_all, xlabels_max, yaxislimits_FS, 1, "Climate Factor of Safety (-)", "maxprotection_allcounties_FS.png")

#violin plot #all counties
#fraction
#violinplots_counties(diff_allRP_percmin_all, xlabels_min, yaxislimits, 0, "Fractional Difference (-)", "minprotection_violin_allcounties.png")
#violinplots_counties(diff_allRP_percmax_all, xlabels_max, yaxislimits, 0, "Fractional Difference (-)", "maxprotection_violin_allcounties.png")
#FS
#violinplots_counties(FS_allRP_percmin_all, xlabels_min, yaxislimits_FS, 1, "Climate Factor of Safety (-)", "minprotection_violin_allcounties_FS.png")
#violinplots_counties(FS_allRP_percmax_all, xlabels_max, yaxislimits_FS, 1,  "Climate Factor of Safety (-)", "maxprotection_violin_allcounties_FS.png")

#labels for positive counties (protected)
xlabels_min = []
xlabels_max = []
for i in range(len(diff_allRP_percmin_pos)):
    xlabels_min.append(RP_list[i] + "\n(" + str(int(len(diff_allRP_percmin_pos[i])/len(county_centroids)*100)) + "%)")
    xlabels_max.append(RP_list[i] + "\n(" + str(int(len(diff_allRP_percmax_pos[i])/len(county_centroids)*100)) + "%)")

#positive only
#fraction
#boxplots_counties(diff_allRP_percmin_pos, xlabels_min, yaxislimits, 0, "Fractional Difference (-)", "minprotection.png")
#boxplots_counties(diff_allRP_percmax_pos, xlabels_max, yaxislimits, 0 ,"Fractional Difference (-)", "maxprotection.png")
#violinplots_counties(diff_allRP_percmin_pos, xlabels_min, "minprotection_violin.png")
#violinplots_counties(diff_allRP_percmax_pos, xlabels_max, "maxprotection_violin.png")
#FS
#boxplots_counties(FS_allRP_percmin_pos, xlabels_min, yaxislimits_FS, 1, "Climate Factor of Safety (-)", "minprotection_FS.png")
#boxplots_counties(FS_allRP_percmax_pos, xlabels_max, yaxislimits_FS, 1, "Climate Factor of Safety (-)", "maxprotection_FS.png")

#labels for negative counties (unprotected)
xlabels_min = []
xlabels_max = []
for i in range(len(diff_allRP_percmin_neg)):
    xlabels_min.append(RP_list[i] + "\n(" + str(int(len(diff_allRP_percmin_neg[i])/len(county_centroids)*100)) + "%)")
    xlabels_max.append(RP_list[i] + "\n(" + str(int(len(diff_allRP_percmax_neg[i])/len(county_centroids)*100)) + "%)")
    
#negative only
#boxplots_counties(diff_allRP_percmin_neg, xlabels_min, yaxislimits, 0, "Fractional Difference (-)", "minprotection_unprotected.png")
#boxplots_counties(diff_allRP_percmax_neg, xlabels_max, yaxislimits, 0 , "Fractional Difference (-)", "maxprotection_unprotected.png")
#violinplots_counties(diff_allRP_percmin_neg, xlabels_min, "minprotection_violin_unprotected.png")
#violinplots_counties(diff_allRP_percmax_neg, xlabels_max, "maxprotection_violin_unprotected.png")

#box plots for each RCP and time period in turn
def countprotected(difflist):
    count = [0, 0, 0, 0, 0, 0]
    for i in range(len(difflist)):
        for j in difflist[i]:
            if j >= 0:
                count[i] += 1
    return(count)

count_45s = countprotected(diff_allRP_percmin_RCP45_short)
xlabels_45s = []
for i in range(len(RP_list)):
    #xlabels_45s.append(RP_list[i] + " (N = " + str(count_45s[i]) + ")")
    xlabels_45s.append(RP_list[i] + "\n(" + str(int((count_45s[i]/len(county_centroids))*100)) + "%)")

count_45l = countprotected(diff_allRP_percmin_RCP45_long)
xlabels_45l = []
for i in range(len(RP_list)):
    #xlabels_45l.append(RP_list[i] + " (N = " + str(count_45l[i]) + ")")
    xlabels_45l.append(RP_list[i] + "\n(" + str(int((count_45l[i]/len(county_centroids))*100)) + "%)")

count_85s = countprotected(diff_allRP_percmin_RCP85_short)
xlabels_85s = []
for i in range(len(RP_list)):
    #xlabels_45s.append(RP_list[i] + " (N = " + str(count_45s[i]) + ")")
    xlabels_85s.append(RP_list[i] + "\n(" + str(int((count_85s[i]/len(county_centroids))*100)) + "%)")

count_85l = countprotected(diff_allRP_percmin_RCP85_long)
xlabels_85l = []
for i in range(len(RP_list)):
    #xlabels_45l.append(RP_list[i] + " (N = " + str(count_45l[i]) + ")")
    xlabels_85l.append(RP_list[i] + "\n(" + str(int((count_85l[i]/len(county_centroids))*100)) + "%)")


#fraction
boxplots_counties(diff_allRP_percmin_RCP45_short, xlabels_45s, yaxislimits, 0, "Fractional Difference (-)", "RCP45_short.png")
boxplots_counties(diff_allRP_percmin_RCP45_long, xlabels_45l, yaxislimits, 0, "Fractional Difference (-)", "RCP45_long.png")
boxplots_counties(diff_allRP_percmin_RCP85_short, xlabels_85s, yaxislimits, 0, "Fractional Difference (-)", "RCP85_short.png")
boxplots_counties(diff_allRP_percmin_RCP85_long, xlabels_85l, yaxislimits, 0, "Fractional Difference (-)", "RCP85_long.png")
#fraction + 1
boxplots_counties(CRF_allRP_percmin_RCP45_short, xlabels_45s, tuple([item + 1 for item in yaxislimits]), 1, "Climate Factor of Safety (-)", "RCP45_short_CRF.png")
boxplots_counties(CRF_allRP_percmin_RCP45_long, xlabels_45l, tuple([item + 1 for item in yaxislimits]), 1, "Climate Factor of Safety (-)", "RCP45_long_CRF.png")
boxplots_counties(CRF_allRP_percmin_RCP85_short, xlabels_85s, tuple([item + 1 for item in yaxislimits]), 1, "Climate Factor of Safety (-)", "RCP85_short_CRF.png")
boxplots_counties(CRF_allRP_percmin_RCP85_long, xlabels_85l, tuple([item + 1 for item in yaxislimits]), 1, "Climate Factor of Safety (-)", "RCP85_long_CRF.png")
#FS
boxplots_counties(FS_allRP_percmin_RCP45_short, xlabels_45s, yaxislimits_FS, 1, "Climate Factor of Safety (-)", "RCP45_short_FS.png")
boxplots_counties(FS_allRP_percmin_RCP45_long, xlabels_45l, yaxislimits_FS, 1, "Climate Factor of Safety (-)", "RCP45_long_FS.png")
boxplots_counties(FS_allRP_percmin_RCP85_short, xlabels_85s, yaxislimits_FS, 1, "Climate Factor of Safety (-)", "RCP85_short_FS.png")
boxplots_counties(FS_allRP_percmin_RCP85_long, xlabels_85l, yaxislimits_FS, 1, "Climate Factor of Safety (-)", "RCP85_long_FS.png")
