# Requires geo-env (laptop)
import numpy as np
import geopandas as gp
import pandas as pd
from datetime import datetime, timedelta
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt


''' This script plot England and Wales regions '''
''' and show the ratio of climate deaths to covid deaths on them '''
''' Eunice Lo '''
''' Created 24/03/2023 '''
''' Updated with results from pop weighted tmean '''
''' Updated 18/04/2024 to save as pdf '''


if __name__ == "__main__":

    # paths
    shp_dir = "ONS/"
    out_dir = "plots/"
    
    # read England regions shapefile (BUC is 500 m resolution)
    eng_df = gp.GeoDataFrame.from_file(shp_dir+"RGN_DEC_2022_EN_BUC.shp")
    eng_df.to_crs(epsg=4326, inplace=True)
    # Wales shapefile
    wales_df = gp.GeoDataFrame.from_file(shp_dir+"CTRY_DEC_2022_GB_BUC.shp")
    wales_df = wales_df[wales_df["CTRY22NM"]=="Wales"]
    wales_df.to_crs(epsg=4326, inplace=True)
    # rename columns to match those in England
    wales_df.rename(columns={"CTRY22CD": "RGN22CD", "CTRY22NM": "RGN22NM"}, inplace=True) 
    # combined England ans Wales shapefile
    eng_wales_df = pd.concat([eng_df, wales_df])
    
    # in the order of regions listed in shapefile, i.e., NE, NW, Y&H, E Mids, W Mids, E Eng, London, SE, SW
    regnames = ["NorthEastEngland", "NorthWestEngland", "YorkshireandHumber", "EastMidlands", \
               "WestMidlands", "EastofEngland", "London", "SouthEastEngland", "SouthWestEngland", \
               "Wales"]

    # mid-2021 population by region, copied from 
    # www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland    
    reg_pops = [2646772, 7422295, 5481431, 4880094, 5954240, 6348096, 8796628, 9294023, 5712840, 3105410]

    # UKHSA-defined heatwaves in period 2020-01-30 to 2022-12-31
    # stats: 10 HWs, total = 70 days
    # first covid death in England and Wales was 2020-01-30 in SE England
    hws = [["2020-06-23", "2020-06-27"], ["2020-07-30", "2020-08-01"], ["2020-08-05", "2020-08-15"], \
           ["2021-07-16", "2021-07-23"], ["2021-09-06", "2021-09-09"], \
           ["2022-06-16", "2022-06-19"], ["2022-07-10", "2022-07-25"], ["2022-07-30", "2022-08-05"], ["2022-08-08", "2022-08-17"], ["2022-08-23", "2022-08-25"]]

    # cold snaps in period 
    # defined here as ave temp <= 2 deg C in any region for 2 consective days, between 1 Nov to 31 Mar (c.f. CWP England)
    # see outputs/daily_tas_by_region_cold_indicator_20200130-20221231.csv
    # stats: 14 CSs, total = 70 days
    #css = [["2020-12-06", "2020-12-08"], ["2020-12-24", "2020-12-25"], ["2020-12-27", "2021-01-10"], \
    #       ["2021-01-12", "2021-01-16"], ["2021-01-23", "2021-01-26"], ["2021-01-30", "2021-02-03"], ["2021-02-07", "2021-02-14"], ["2021-03-06", "2021-03-07"], ["2021-11-27", "2021-11-28"], ["2021-12-19", "2021-12-22"], \
    #       ["2022-01-04", "2022-01-07"], ["2022-01-14", "2022-01-15"], ["2022-01-20", "2022-01-21"], ["2022-12-06", "2022-12-17"]] 

    # cold waves in period
    # defined here as dates on which a L3 Cold Weather Alert was issued for any region in England (see CWP England)
    # reference: Cold and heat alerts 2016-2022.xlsx
    # stats: 8 cold waves, total = 70 days; double checked
    css = [["2020-12-29", "2021-01-18"], \
           ["2021-01-22", "2021-02-02"], ["2021-02-08", "2021-02-12"], ["2021-11-26", "2021-11-29"], ["2021-12-20", "2021-12-23"], \
           ["2022-01-04", "2022-01-10"], ["2022-01-13", "2022-01-17"], ["2022-12-07", "2022-12-18"]] 
    
    # covid and temperature-related deaths, 2020-01-30 to 2022-12-31 
    covid_all_ds = []   # all days in period
    temp_all_ds = []
    covid_hw_ds = []    # heatwaves
    temp_hw_ds = []
    covid_cs_ds = []    # cold snaps
    temp_cs_ds = []
    
    for reg in regnames:
    
        print(reg) 

        # cumulated covid deaths by region
        # these are deaths with COVID-19 on the death certificate
        # downloaded data: /Users/yl17544/OneDrive - University of Bristol/ONS_data/daily_deaths/daily_covid_deaths_on_certificate_202003_202302
        # source: https://coronavirus.data.gov.uk/details/deaths?areaType=region&areaName=East%20Midlands
        covid_df = pd.read_csv("data/covid_deaths_on_cert_data_"+reg+"_2023-Mar-23.csv", \
                   usecols=["date", "newDailyNsoDeathsByDeathDate", "cumDailyNsoDeathsByDeathDate"])
        covid_df = covid_df.loc[covid_df["date"]<="2022-12-31"]
        # reverse in time
        covid_df.sort_values(by="date", inplace=True)       
        # add dates from 2020-01-30 if doesn't exist
        date1 = datetime.strptime(covid_df["date"].iloc[0], '%Y-%m-%d').date()
        sdate = datetime.strptime("2020-01-30", '%Y-%m-%d').date()
        if date1 > sdate:
           edate = date1 - timedelta(days=1)    
           extras = pd.date_range(start=sdate, end=edate)        	    # extra dates
           extras2 = [date.strftime('%Y-%m-%d') for date in extras]
           extra_ds = [0] * len(extras2)                                    # extra daily deaths, all 0
           extra_data = {"date":extras2, "newDailyNsoDeathsByDeathDate":extra_ds, "cumDailyNsoDeathsByDeathDate":extra_ds}
           extra_df = pd.DataFrame(data=extra_data)       
           covid_df = pd.concat([extra_df, covid_df])
        # remove leap days
        covid_df = covid_df[~covid_df.date.str.endswith("02-29")]   

        # temperature-related deaths from ER modelling
        temp_df = pd.read_csv("outputs/daily_attributable_deaths_ER1981-2022_yearround_21dayslag_MMT2to98_nsim100_1981-2022_ONSdata_"+\
                              reg+".csv", usecols=["date","tmean","mmt","est"])
        temp_df = temp_df.loc[(temp_df["date"]>="2020-01-30") & (temp_df["date"]<="2022-12-31")]
        # isolate non-optimal temperature deaths
        temp_df.loc[temp_df.tmean == temp_df.mmt, "est"] = 0     # best estimate
  
        # the whole period
        covid_all_ds.append(int(covid_df["newDailyNsoDeathsByDeathDate"].sum()))
        temp_all_ds.append(int(temp_df["est"].sum()))

        # only heatwave days in period
        covid_hwd = 0
        temp_hwd = 0
        for hw in hws:
            # subset heatwave dates
            covid_hw_df = covid_df.loc[(covid_df["date"]>=hw[0]) & (covid_df["date"]<=hw[1])]
            covid_hwd += int(covid_hw_df["newDailyNsoDeathsByDeathDate"].sum())
            temp_hw_df = temp_df.loc[(temp_df["date"]>=hw[0]) & (temp_df["date"]<=hw[1])]
            temp_hwd += int(temp_hw_df["est"].sum())
        covid_hw_ds.append(covid_hwd)
        temp_hw_ds.append(temp_hwd)    

        # only coldsnaps in period
        covid_csd = 0
        temp_csd = 0
        for cs in css:
            # subset coldsnap dates
            covid_cs_df = covid_df.loc[(covid_df["date"]>=cs[0]) & (covid_df["date"]<=cs[1])]
            covid_csd += int(covid_cs_df["newDailyNsoDeathsByDeathDate"].sum())
            temp_cs_df = temp_df.loc[(temp_df["date"]>=cs[0]) & (temp_df["date"]<=cs[1])]
            temp_csd += int(temp_cs_df["est"].sum())
        covid_cs_ds.append(covid_csd)     # checked with Wales value here vs covid spreahsheet manually!
        temp_cs_ds.append(temp_csd)
    
    # put data in shapes dataframe
    # ratios
    eng_wales_df["tempcovid_death_ratio"] = np.array(temp_all_ds)/np.array(covid_all_ds)
    eng_wales_df["hw_tempcovid_death_ratio"] = np.array(temp_hw_ds)/np.array(covid_hw_ds)
    eng_wales_df["cs_tempcovid_death_ratio"] = np.array(temp_cs_ds)/np.array(covid_cs_ds)
       
    # plot regions on map     
    font = {'size' : 12}
    plt.rc('font', **font)

    # set discrete colourmap for ratios
    c_below1 = mpl.cm.Greys(np.linspace(0.2, 0.6, 5))
    c_above1 = mpl.cm.BuPu(np.linspace(0.25, 1, 9))
    c_both = np.vstack((c_below1, c_above1))
    cmap = mcolors.LinearSegmentedColormap.from_list('my_colormap', c_both)
    norm = mpl.colors.BoundaryNorm(np.arange(0,15)/5, cmap.N) 
    
    f, axs = plt.subplots(1, 3, figsize=(12,4))
 
    f.subplots_adjust(bottom=0, top=0.9, left=0, right=0.9,
                    wspace=0.01, hspace=0.01) 

    # whole period ratio
    eng_wales_df.plot(column="tempcovid_death_ratio", ax=axs[0], edgecolor="black", cmap=cmap, norm=norm, \
                      cax=axs[0], legend=False)
    axs[0].title.set_text(r"$\bf{a}$"+"  30 Jan 2020 to 31 Dec 2022")

    # heatwaves in period ratio
    eng_wales_df.plot(column="hw_tempcovid_death_ratio", ax=axs[1], edgecolor="black", cmap=cmap, norm=norm, \
                      cax=axs[1], legend=False)
    axs[1].title.set_text(r"$\bf{b}$"+"  Heatwaves (70 days)")

    # cold snaps in period ratio
    eng_wales_df.plot(column="cs_tempcovid_death_ratio", ax=axs[2], edgecolor="black", cmap=cmap, norm=norm, \
                      cax=axs[2], legend=False)
    axs[2].title.set_text(r"$\bf{c}$"+"  Cold snaps (70 days)")
    
    # add colourbar for ratio, same as above
    col = axs[0].collections[0]
    cb_ax = f.add_axes([0.9, 0.04, 0.02, 0.9])        # x0, y0, width, height
    cb = f.colorbar(col, cax=cb_ax, fraction=0.05)
    cb.ax.set_ylabel("Temperature to COVID deaths ratio")
    
    for ax in axs:     # axs.flatten():
        ax.axis("off")

    plt.savefig(out_dir+"EngWales_tempcovid_mortality_fraction_202001-202212_all_heatwaves_coldwaves_officaldef_ONSdata.pdf", format="pdf", dpi=300)
    plt.close("all")

    print("Saved graph!")
