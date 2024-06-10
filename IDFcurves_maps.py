# -*- coding: utf-8 -*-
"""
Created on Sun May 19 17:35:26 2024

@author: webbe

creates maps for IDF curves 2024 for Chesapeake Bay watershed counties
using the data as directly downloaded from the tool: https://midatlantic-idf.rcc-acis.org/

"""

import pandas as pd
import numpy as np
#import math
import statistics as st
import matplotlib.pyplot as plt
import scipy.stats as stats
#import os
#import requests 
import geopandas as gpd
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
pd.options.mode.chained_assignment = None  # default='warn'

#%% paths and other user defined variables

path_main = r"/Users/webbe/Box/Marissa's Research/IDF Curves/IDFcurve_code/"
path1 = path_main + r"atlasmedian/"
path2 = path_main + r"atlasCI/"
path3 = path_main + r"atlasupRP/"
path4 = path_main + r"atlasupRPCI/"
path_to_save = path_main + r"maps/"
path_to_save_atlas = path_main + r"atlas14/"

RCP = ["rcp45_2020-2070", "rcp45_2050-2100", "rcp85_2020-2070", "rcp85_2050-2100"]
RP_list = ["2-yr", "5-yr", "10-yr", "25-yr", "50-yr", "100-yr"]
ptiles = ["mean", "min", "10th ptile", "25th ptile", "50th ptile", "75th ptile", "90th ptile", "max"]

#%% introductory maps of all US and CBW within US
#use geopandas to plotmap of CBW counties
county_centroids = pd.read_csv(path_main + "CBP_countycentroids.csv")

#list to selectcounties with GEOID from county_centroids
geoid_list = []
for i in range(len(county_centroids)):
    geoid_list.append(str(county_centroids['GEOID'][i]))

#all US
allUScounties = gpd.read_file(path_main + "cb_2018_us_county_20m/cb_2018_us_county_20m.shp")
allUScounties['color'] = [1]*len(allUScounties)
allUScounties.plot(column = 'color', cmap ='Reds', edgecolor='0')
#plt.savefig(path_to_save + r"USmap.png", bbox_inches='tight')
plt.close()

#CONUS only
nonCONUSstates = ['02', '03', '07','14', '15', '43', '72'] #AL, old AR, old CT, old IA, HI, old UT, PR
CONUScounties = allUScounties[~allUScounties['STATEFP'].isin(nonCONUSstates)].reset_index()
CONUScounties.plot(color = "silver")
#plt.savefig(path_to_save + r"CONUSmap.png", bbox_inches='tight')
plt.close()

#CONUS with CBW highlighted
CBWinCONUS = []
for i in range(len(CONUScounties)):
    if CONUScounties['GEOID'][i] in geoid_list:
        CBWinCONUS.append(1)
    else:
        CBWinCONUS.append(0)
CONUScounties['CBW'] = CBWinCONUS
mycmap = ListedColormap(['lightgray', 'brown'])
CONUScounties.plot(column = 'CBW', cmap = mycmap)#, edgecolor='0.5', linewidth = 0.5)
#plt.savefig(path_to_save + r"CONUSmap+CBW.png", bbox_inches='tight')
plt.close()

#CBW only
CBWcounties = allUScounties[allUScounties['GEOID'].isin(geoid_list)] #select counties with GEOID from county_centroids
CBWcounties['color'] = [1]*len(CBWcounties)
CBWcounties.plot(color = "grey")#, edgecolor='0.5', linewidth = 0.5)
#plt.savefig(path_to_save + r"CBWmap.png", bbox_inches='tight')
plt.close()
CBWcounties = allUScounties[allUScounties['GEOID'].isin(geoid_list)]

#create an inset map of CBW within CONUS
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)
ax.set_xticks([])
ax.set_yticks([])
CONUScounties.plot(column = 'CBW', cmap = mycmap, ax = ax)
ax2 = inset_axes(ax, "30%", "40%", loc="lower left")
CBWcounties.plot(color = "brown", ax = ax2)
ax2.xaxis.set_ticks_position('top')
ax2.xaxis.set_label_position('top') 
plt.savefig(path_to_save + r"CONUSmap+CBW+inset.png", bbox_inches='tight')
plt.close()

#%% create new color maps for the 10 panel maps

#new green - that does not start from white #34,139,34 #144, 238, 144
#N = 256
#vals = np.ones((N, 4))
#vals[:, 0] = np.linspace(144/256, 34/256, N)
#vals[:, 1] = np.linspace(238/256, 139/256, N)
#vals[:, 2] = np.linspace(144/256, 34/356, N)
#newgreen = ListedColormap(vals)
N = 256 #1, 50, 32
vals = np.ones((N, 4))
vals[:, 0] = np.linspace(144/256, 1/256, N)
vals[:, 1] = np.linspace(238/256, 50/256, N)
vals[:, 2] = np.linspace(144/256, 32/356, N)
vals[:1, :] = np.linspace(128/256, 128/256, 1) #grey
#vals[:1, :] = np.linspace(1, 1, 1) #white
newgreen = ListedColormap(vals)

#new blue - that does not start from white #0,0,255 #173, 216, 230
vals2 = np.ones((N, 4))
vals2[:, 0] = np.linspace(173/256, 0/256, N)
vals2[:, 1] = np.linspace(216/256, 0/256, N)
vals2[:, 2] = np.linspace(230/256, 255/256, N)
vals2[:1, :] = np.linspace(128/256, 128/256, 1) #grey
#vals2[:1, :] = np.linspace(1, 1, 1) #white
newblue = ListedColormap(vals2)

#%% 10 panel maps
#Description
##12 maps in total: each of the pairs of RCP and future periods (RCP4.5 together, RCP 8.5 together) 
#has 6 maps for each of the 6 return periods
#each map has 10 panels: for the 4 strategies + 1 for the RCP depth-Atlas14 median x2
def multiplemaps_byRP_10paneldiff(outdf, outputrcp1, outputrcp2, t1, t2, legendtitle, titletosave):
    fig, axs = plt.subplots(2, 5, figsize = (18, 6), dpi = 300)
    axs = axs.flatten()
    vmin, vmax = np.min(outdf.values), np.max(outdf.values)
    for i in list(range(0, 4)) + list(range(5, 9)):
        if i <= 4:
            CBWcounties[outdf.columns[i]] = list(outdf[outdf.columns[i]])
            variable = outdf.columns[i]
            CBWcounties.plot(column = variable, cmap=newgreen, vmin =vmin, vmax=vmax, linewidth=0.8, ax=axs[i], edgecolor='0.8', legend = False)
            #CBWcounties.plot(column = variable, cmap='BrBG', vmin =vmin, vmax=vmax, linewidth=0.8, ax=axs[i], edgecolor='0.8', legend = False)
            #axs[i].set_facecolor("grey")
            axs[i].set_title(outdf.columns[i], fontsize = 14)
            axs[i].set_xticks([])
            axs[i].set_yticks([])
            if i == 0:
                axs[i].set_ylabel("Projected climate \n 2020-2070", fontsize = 14)
                axs[i].text(0.1, 0.8, t1[i] + "% of counties \nclimate ready", transform=axs[i].transAxes, size = 14, linespacing=1)
            else:
                axs[i].text(0.1, 0.8, t1[i] + "% \n ", transform=axs[i].transAxes, size = 14, linespacing=1)
                
        else:
            CBWcounties[outdf.columns[i-1]] = list(outdf[outdf.columns[i-1]])
            variable = outdf.columns[i-1]
            CBWcounties.plot(column = variable, cmap=newgreen, vmin =vmin, vmax=vmax, linewidth=0.8, ax=axs[i], edgecolor='0.8', legend = False)
            #CBWcounties.plot(column = variable, cmap='BrBG', vmin =vmin, vmax=vmax, linewidth=0.8, ax=axs[i], edgecolor='0.8', legend = False)
            #axs[i].set_facecolor("grey")
            axs[i].set_xticks([])
            axs[i].set_yticks([])
            if i == 5:
                axs[i].set_ylabel("Projected climate \n 2050-2100", fontsize = 14)
                axs[i].text(0.1, 0.8, t1[i] + "% of counties \nclimate ready", transform=axs[i].transAxes, size = 14, linespacing=1)
            else:
                axs[i].text(0.1, 0.8, t1[i] + "% \n ", transform=axs[i].transAxes, size = 14, linespacing=1)
    sm1 = plt.cm.ScalarMappable(cmap=newgreen, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    #sm1 = plt.cm.ScalarMappable(cmap='BrBG', norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm1._A = [] #empty array for the data range
    #cbar_ax = fig.add_axes([0.125, 0.52, 0.535, 0.04]) #specify location of the colorbar
    cbar_ax = fig.add_axes([0.125, 0, 0.535, 0.04]) #specify location of the colorbar
    cbar = fig.colorbar(sm1, cax = cbar_ax, orientation = "horizontal", extend = "min")#, extendrect = "True")
    cbar.set_label(label = "Stormwater Infrastructure Climate Factor of Safety (" + legendtitle + ")", size=16)#, weight='bold')
    cbar.ax.tick_params(labelsize=14) 
    CBWcounties.loc[:, "RCP_2020-2070"] = outputrcp1
    CBWcounties.loc[:, "RCP_2050-2100"] = outputrcp2
    #vmin2, vmax2 = min(outputrcp1), max(outputrcp2)
    vmin2, vmax2 = 0, max(outputrcp2) #ensude the difference starts from zero
    CBWcounties.plot(column = 'RCP_2020-2070', cmap=newblue, vmin=vmin2, vmax=vmax2, linewidth=0.8, ax=axs[4], edgecolor='0.8')#, legend = True)
    axs[4].set_title("Projected \nAdditional Rainfall", fontsize = 14)
    axs[4].set_xticks([])
    axs[4].set_yticks([])
    #sm = plt.cm.ScalarMappable(cmap='Blues', norm=plt.Normalize(vmin=vmin2, vmax=vmax2))
    #sm._A = [] #empty array for the data range
    #cbar = fig.colorbar(sm, ax = axs[4], label = "Depth (mm)") #add the colorbar to the figure
    CBWcounties.plot(column = 'RCP_2050-2100', cmap=newblue, vmin=vmin2, vmax=vmax2, linewidth=0.8, ax=axs[9], edgecolor='0.8')#, legend = True)
    #axs[9].set_title("RCP projection")
    axs[9].set_xticks([])
    axs[9].set_yticks([])
    sm = plt.cm.ScalarMappable(cmap=newblue, norm=plt.Normalize(vmin=vmin2, vmax=vmax2))
    sm._A = [] #empty array for the data range
    #cbar_ax2 = fig.add_axes([0.68, 0.52, 0.125, 0.04]) #specify location of the colorbar
    cbar_ax2 = fig.add_axes([0.68, 0, 0.125, 0.04]) #specify location of the colorbar
    cbar2 = fig.colorbar(sm, cax = cbar_ax2, orientation = "horizontal", extend = "min")#, extendrect = "True") #add the colorbar to the figure
    cbar2.set_label(label = "Rainfall depth (mm)", size=16)#, weight='bold')
    cbar2.ax.tick_params(labelsize=14) 
    #cbar_ax2 = fig.add_axes([0.8, 0.15, 0.02, 0.7])
    #cbar2 = fig.colorbar(sm, cax = cbar_ax2, label = "Depth (mm)", orientation = "vertical") #add the colorbar to the figure
    fig.subplots_adjust(right=0.8)
    plt.savefig(titletosave, bbox_inches='tight')
    
def multiplemaps_byRP_10paneldiff_100(outdf, outputrcp1, outputrcp2, t1, t2, legendtitle, titletosave):
    fig, axs = plt.subplots(2, 5, figsize = (18, 6), dpi = 300)
    axs = axs.flatten()
    vmin, vmax = np.min(outdf.dropna(axis = 1).values), np.max(outdf.dropna(axis = 1).values)
    #vmin, vmax = 0, np.max(outdf.values)
    for i in list(range(0, 4)) + list(range(5, 9)):
        if i <= 4:
            CBWcounties[outdf.columns[i]] = list(outdf[outdf.columns[i]])
            variable = outdf.columns[i]
            CBWcounties.plot(column = variable, cmap=newgreen, vmin=vmin, vmax=vmax, linewidth=0.8, ax=axs[i], edgecolor='0.8', legend = False)
            #CBWcounties.plot(column = variable, cmap='BrBG', vmin =vmin, vmax=vmax, linewidth=0.8, ax=axs[i], edgecolor='0.8', legend = False)
            axs[i].set_xticks([])
            axs[i].set_yticks([])
            axs[i].set_title(outdf.columns[i], fontsize = 14)
            if i == 0:
                axs[i].set_ylabel("Projected climate \n 2050-2100", fontsize = 14)
                axs[i].text(0.1, 0.8, t1[i] + "% of counties \nclimate ready", transform=axs[i].transAxes, size = 14, linespacing=1)
            else:
                axs[i].text(0.1, 0.8, t1[i] + "% \n", transform=axs[i].transAxes, size = 14, linespacing=1)
        else:
            CBWcounties[outdf.columns[i-1]] = list(outdf[outdf.columns[i-1]])
            variable = outdf.columns[i-1]
            CBWcounties.plot(column = variable, cmap=newgreen, vmin=vmin, vmax=vmax, linewidth=0.8, ax=axs[i], edgecolor='0.8', legend = False)
            #CBWcounties.plot(column = variable, cmap='BrBG', vmin =vmin, vmax=vmax, linewidth=0.8, ax=axs[i], edgecolor='0.8', legend = False)
            axs[i].set_xticks([])
            axs[i].set_yticks([])
            if i == 5:
                axs[i].set_ylabel("Projected climate \n 2020-2070", fontsize = 14)
                axs[i].text(0.1, 0.8, t1[i] + "% of counties \nclimate ready", transform=axs[i].transAxes, size = 14, linespacing=1)
            else:
                axs[i].text(0.1, 0.8, t1[i] + "% \n", transform=axs[i].transAxes, size = 14, linespacing=1)
    sm1 = plt.cm.ScalarMappable(cmap=newgreen, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    #sm1 = plt.cm.ScalarMappable(cmap='BrBG', norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm1._A = [] #empty array for the data range
    #cbar_ax = fig.add_axes([0.125, 0.52, 0.535, 0.04]) #specify location of the colorbar
    cbar_ax = fig.add_axes([0.125, 0, 0.535, 0.04]) #specify location of the colorbar
    cbar = fig.colorbar(sm1, cax = cbar_ax, orientation = "horizontal", extend = "min")#, extendrect = "True")
    cbar.set_label(label = "Stormwater Infrastructure Climate Factor of Safety (" + legendtitle + ")", size=16)#, weight='bold')
    cbar.ax.tick_params(labelsize=14) 
    CBWcounties.loc[:, "RCP_2020-2070"] = outputrcp1
    CBWcounties.loc[:, "RCP_2050-2100"] = outputrcp2
    #vmin2, vmax2 = min(outputrcp1), max(outputrcp2)
    vmin2, vmax2 = 0, max(outputrcp2) #ensure the difference starts from zero
    CBWcounties.plot(column = 'RCP_2020-2070', cmap=newblue, vmin=vmin2, vmax=vmax2, linewidth=0.8, ax=axs[4], edgecolor='0.8')#, legend = True)
    axs[4].set_title("Projected \nAdditional Rainfall", fontsize = 14)
    axs[4].set_xticks([])
    axs[4].set_yticks([])
    #sm = plt.cm.ScalarMappable(cmap='Blues', norm=plt.Normalize(vmin=vmin2, vmax=vmax2))
    #sm._A = [] #empty array for the data range
    #cbar = fig.colorbar(sm, ax = axs[4], label = "Depth (mm)") #add the colorbar to the figure
    CBWcounties.plot(column = 'RCP_2050-2100', cmap=newblue, vmin=vmin2, vmax=vmax2, linewidth=0.8, ax=axs[9], edgecolor='0.8')#, legend = True)
    #axs[9].set_title("RCP projection")
    axs[9].set_xticks([])
    axs[9].set_yticks([])
    sm = plt.cm.ScalarMappable(cmap=newblue, norm=plt.Normalize(vmin=vmin2, vmax=vmax2))
    sm._A = [] #empty array for the data range
    #cbar_ax2 = fig.add_axes([0.68, 0.52, 0.125, 0.04]) #specify location of the colorbar
    cbar_ax2 = fig.add_axes([0.68, 0, 0.125, 0.04]) #specify location of the colorbar
    cbar2 = fig.colorbar(sm, cax = cbar_ax2, orientation = "horizontal", extend = "min")#, extendrect = "True") #add the colorbar to the figure
    cbar2.set_label(label = "Rainfall depth (mm)", size=16)#, weight='bold')
    cbar2.ax.tick_params(labelsize=14)
    #cbar_ax2 = fig.add_axes([0.8, 0.15, 0.02, 0.7])
    #cbar2 = fig.colorbar(sm, cax = cbar_ax2, label = "Depth (mm)", orientation = "vertical") #add the colorbar to the figure
    axs[-2].set_visible(False) #axs[-2].axis('off')
    axs[-3].set_visible(False) #axs[-3].axis('off')
    axs[2].set_visible(False) #axs[2].axis('off')
    axs[3].set_visible(False) #axs[3].axis('off')
    fig.subplots_adjust(right=0.8)
    plt.savefig(titletosave, bbox_inches='tight')
    
def importandmapoutput_byRP_10paneldiff(path_list, filename1, filename2, filename3, filename4, titletosave_list):
    out1 = pd.read_csv(path_list[0] + r"results/" + filename1).drop("Unnamed: 0", axis = 1)
    out2 = pd.read_csv(path_list[1] + r"results/" + filename1).drop("Unnamed: 0", axis = 1)
    out3 = pd.read_csv(path_list[2] + r"results/" + filename1).drop("Unnamed: 0", axis = 1)
    out4 = pd.read_csv(path_list[3] + r"results/" + filename1).drop("Unnamed: 0", axis = 1)
    out5 = pd.read_csv(path_list[0] + r"results/" + filename2).drop("Unnamed: 0", axis = 1)
    out6 = pd.read_csv(path_list[1] + r"results/" + filename2).drop("Unnamed: 0", axis = 1)
    out7 = pd.read_csv(path_list[2] + r"results/" + filename2).drop("Unnamed: 0", axis = 1)
    out8 = pd.read_csv(path_list[3] + r"results/" + filename2).drop("Unnamed: 0", axis = 1)
    outrcp1 = pd.read_csv(path_list[0] + r"results/" + filename3).drop("Unnamed: 0", axis = 1)
    outrcp2 = pd.read_csv(path_list[0] + r"results/" + filename4).drop("Unnamed: 0", axis = 1)
    outputatlas = pd.read_csv(path_list[0]+ r"results/output_atlas_CBW.csv").drop("Unnamed: 0", axis = 1)
    outputatlasCI = pd.read_csv(path_list[0]+ r"results/output_atlas_max.csv").drop("Unnamed: 0", axis = 1)
    for i in range(len(out1.columns)-1):
        output = pd.DataFrame(
            {"Current " + RP_list[i] + " \nmedian": [1+x for x in out1[out1.columns[i]]],
             "Current " + RP_list[i] + " \nupper bound": [1+x for x in out2[out2.columns[i]]],
             "Current " + RP_list[i+1] + " \nmedian": [1+x for x in out3[out3.columns[i]]],
             "Current " + RP_list[i+1] + " \nupper bound": [1+x for x in out4[out4.columns[i]]],
             "Current " + RP_list[i] + " \nmedian ": [1+x for x in out5[out5.columns[i]]],
             "Current " + RP_list[i] + " \nupper bound ": [1+x for x in out6[out6.columns[i]]],
             "Current " + RP_list[i+1] + " \nmedian ": [1+x for x in out7[out7.columns[i]]],
             "Current " + RP_list[i+1] + " \nupper bound ": [1+x for x in out8[out8.columns[i]]]
             })
        output = output.clip(lower = 1.00) #set lower limit for factor of safety = 1
        t1 = [] #text for figures
        t2 = [] #text for figures
        tr1 = []
        tr2 = []
        tr3 = []
        tr4 = []
        tr5 = []
        tr6 = []
        tr7 = []
        tr8 = []
        for c in range(len(output.columns)):
            t1.append(str(int(((sum(x > 1.0 for x in output[output.columns[c]]))/321)*100)))
            t2.append(str(round(st.median(output[output.columns[c]]), 2)))
        tr1.append(stats.ks_2samp(outputatlas[outputatlas.columns[i]], outrcp1[outrcp1.columns[i]], alternative = "less")[1])
        tr2.append(stats.ks_2samp(outputatlas[outputatlas.columns[i]], outrcp2[outrcp2.columns[i]], alternative = "less")[1])
        tr3.append(stats.ks_2samp(outputatlasCI[outputatlasCI.columns[i]], outrcp1[outrcp1.columns[i]], alternative = "less")[1])
        tr4.append(stats.ks_2samp(outputatlasCI[outputatlasCI.columns[i]], outrcp2[outrcp2.columns[i]], alternative = "less")[1])
        tr5.append(stats.ks_2samp(outputatlas[outputatlas.columns[i+1]], outrcp1[outrcp1.columns[i]], alternative = "less")[1])
        tr6.append(stats.ks_2samp(outputatlas[outputatlas.columns[i+1]], outrcp2[outrcp2.columns[i]], alternative = "less")[1])
        tr7.append(stats.ks_2samp(outputatlasCI[outputatlasCI.columns[i+1]], outrcp1[outrcp1.columns[i]], alternative = "less")[1])
        tr8.append(stats.ks_2samp(outputatlasCI[outputatlasCI.columns[i+1]], outrcp2[outrcp2.columns[i]], alternative = "less")[1])
        #TODO repeat for the other raw depths, and save in a dataframe that will be exported as a csv
        t1.insert(4, 'NA')
        t2.insert(4, 'NA')
        t1.insert(9, 'NA')
        t2.insert(9, 'NA')
        #pd.DataFrame([tr1, tr2]).to_csv(path_to_save + RP_list[i] + "atlas14med_KS_2sampletest.csv")
        #pd.DataFrame([tr3, tr4]).to_csv(path_to_save + RP_list[i] + "atlas14CI_KS_2sampletest.csv")
        #pd.DataFrame([tr5, tr6]).to_csv(path_to_save + RP_list[i] + "atlas14upRP_KS_2sampletest.csv")
        #pd.DataFrame([tr7, tr8]).to_csv(path_to_save + RP_list[i] + "atlas14upRPCI_KS_2sampletest.csv")
        #pd.DataFrame([tr1, tr3, tr5, tr7]).to_csv(path_to_save + RP_list[i] + "KS_2sampletest_midcentury.csv")
        #pd.DataFrame([tr2, tr4, tr6, tr8]).to_csv(path_to_save + RP_list[i] + "KS_2sampletest_endcentury.csv")
        scenario = ""
        if filename1[-11:] == "45short.csv":
            scenario = "RCP45"
        elif filename1[-11:] == "85short.csv":
            scenario = "RCP85"
        multiplemaps_byRP_10paneldiff(output, 
                                  [xi - yi for xi, yi in zip(outrcp1[outrcp1.columns[i]], outputatlas[outputatlas.columns[i]])],
                                  [xi - yi for xi, yi in zip(outrcp2[outrcp2.columns[i]], outputatlas[outputatlas.columns[i]])], 
                                  t1, t2, RP_list[i], titletosave_list[i]+scenario+".png")
    output100 = pd.DataFrame( #NB no up RP for this one
        {"Current 100-yr \nmedian": [1+x for x in out1[out1.columns[-1]]],
         "Current 100-yr \nupper bound": [1+x for x in out2[out2.columns[-1]]],
         "": [-1]*321,#[np.nan]*321,
         " ": [-1]*321, #[np.nan]*321,
         "Current 100-yr \nmedian ": [1+x for x in out5[out5.columns[-1]]],
         "Current 100-yr \nupper bound ": [1+x for x in out6[out6.columns[-1]]],
         "  ": [-1]*321, #[np.nan]*321,
         "   ": [-1]*321#[np.nan]*321
         })                  
    output100 = output100.clip(lower = 1.00) #set lower limit for percentage difference = 0          
    t1 = [] #text for figures
    t2 = [] #text for figures
    tr1 = []
    tr2 = []
    tr3 = []
    tr4 = []
    for c in range(len(output100.columns)):
        t1.append(str(int(((sum(x > 1.0 for x in output100[output100.columns[c]]))/321)*100)))
        t2.append(str(round(st.median(output100[output100.columns[c]]), 2)))
    tr1.append(stats.ks_2samp(outputatlas[outputatlas.columns[-1]], outrcp1[outrcp1.columns[-1]], alternative = "less")[1])
    tr2.append(stats.ks_2samp(outputatlas[outputatlas.columns[-1]], outrcp2[outrcp2.columns[-1]], alternative = "less")[1])
    tr3.append(stats.ks_2samp(outputatlasCI[outputatlasCI.columns[-1]], outrcp1[outrcp1.columns[-1]], alternative = "less")[1])
    tr4.append(stats.ks_2samp(outputatlasCI[outputatlasCI.columns[-1]], outrcp2[outrcp2.columns[-1]], alternative = "less")[1])
    t1.insert(4, 'NA')
    t2.insert(4, 'NA')
    t1.insert(9, 'NA')
    t2.insert(9, 'NA')
    #pd.DataFrame([tr1, tr2]).to_csv(path_to_save + RP_list[-1] + "atlas14med_KS_2sampletest.csv")
    #pd.DataFrame([tr3, tr4]).to_csv(path_to_save + RP_list[-1] + "atlas14CI_KS_2sampletest.csv")
    #pd.DataFrame([tr1, tr3]).to_csv(path_to_save + RP_list[-1] + "KS_2sampletest_midcentury.csv")
    #pd.DataFrame([tr2, tr4]).to_csv(path_to_save + RP_list[-1] + "KS_2sampletest_endcentury.csv")
    scenario = ""
    if filename1[-11:] == "45short.csv":
        scenario = "RCP45"
    elif filename1[-11:] == "85short.csv":
        scenario = "RCP85"
    multiplemaps_byRP_10paneldiff_100(output100, 
                                  [xi - yi for xi, yi in zip(outrcp1[outrcp1.columns[-1]], outputatlas[outputatlas.columns[-1]])],
                                  [xi - yi for xi, yi in zip(outrcp2[outrcp2.columns[-1]], outputatlas[outputatlas.columns[-1]])],  
                                  t1, t2, RP_list[-1], titletosave_list[-1]+scenario+".png")
    
pathlist = [path1, path2, path3, path4]
titletosave_list = [path_to_save + "2yr_fractionmap+diff_",
                    path_to_save + "5yr_fractionmap+diff_",
                    path_to_save + "10yr_fractionmap+diff_",
                    path_to_save + "25yr_fractionmap+diff_",
                    path_to_save + "50yr_fractionmap+diff_",
                    path_to_save + "100yr_fractionmap+diff_"]
importandmapoutput_byRP_10paneldiff(pathlist, "output_percmin_85short.csv", "output_percmin_85long.csv", "output_85short.csv", "output_85long.csv", titletosave_list)
importandmapoutput_byRP_10paneldiff(pathlist, "output_percmin_45short.csv", "output_percmin_45long.csv", "output_45short.csv", "output_45long.csv", titletosave_list)
