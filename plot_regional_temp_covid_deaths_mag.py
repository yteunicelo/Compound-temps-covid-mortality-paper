# Requires conda activate geo_env (local computer)
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
''' Created 04/01/2023, from the corresponding ratios script'''
''' Updated 21/03/2024, with pop w tmean results '''
''' Updated 22/03/2024, with mid-2015 pop estimates for 2010-19 results '''
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
    # mid-2015 population by region, copied from the same source
    reg_pops_15 = [2624579, 7175178, 5390211, 4677425, 5755032, 6075970, 8666930, 8949392, 5471610, 3099086]

    # UKHSA-defined heatwaves in period 2020-01-30 to 2022-12-31
    # stats: 10 HWs, total = 70 days
    # first covid death in England and Wales was 2020-01-30 in SE England
    hws = [["2020-06-23", "2020-06-27"], ["2020-07-30", "2020-08-01"], ["2020-08-05", "2020-08-15"], \
           ["2021-07-16", "2021-07-23"], ["2021-09-06", "2021-09-09"], \
           ["2022-06-16", "2022-06-19"], ["2022-07-10", "2022-07-25"], ["2022-07-30", "2022-08-05"], ["2022-08-08", "2022-08-17"], ["2022-08-23", "2022-08-25"]]

    # cold waves in period
    # defined here as dates on which a L3 Cold Weather Alert was issued for any region in England (see CWP England)
    # reference: Cold and heat alerts 2016-2022.xlsx
    # stats: 8 cold waves, total = 70 days; double checked
    css = [["2020-12-29", "2021-01-18"], \
           ["2021-01-22", "2021-02-02"], ["2021-02-08", "2021-02-12"], ["2021-11-26", "2021-11-29"], ["2021-12-20", "2021-12-23"], \
           ["2022-01-04", "2022-01-10"], ["2022-01-13", "2022-01-17"], ["2022-12-07", "2022-12-18"]] 
  
    # PHE-defined heatwaves in period 2010-2019, with the ones before 2016 backdated by me (see Lo et al., 2022)
    # None in 2010 (see Lo_2021_Englandheatwavedeaths_outline.docx)
    hws_prev = [["2011-06-25", "2011-06-28"], ["2011-07-31", "2011-08-02"], \
                ["2012-08-17", "2012-08-20"], \
                ["2013-07-11", "2013-07-14"], ["2013-07-16", "2013-07-24"], ["2013-07-31", "2013-08-02"], \
                ["2014-07-17", "2014-07-20"], ["2014-07-22", "2014-07-27"], \
                ["2015-06-29", "2015-07-02"], ["2015-08-21", "2015-08-23"], \
                ["2016-07-19", "2016-07-21"], ["2016-08-23", "2016-08-25"], ["2016-09-14", "2016-09-16"], \
                ["2017-06-16", "2017-06-23"], ["2017-07-05", "2017-07-07"], \
                ["2018-06-25", "2018-06-27"], ["2018-06-30", "2018-07-10"], ["2018-07-21", "2018-07-29"], ["2018-08-01","2018-08-09"], \
                ["2019-06-28", "2019-06-30"], ["2019-07-21", "2019-07-28"], ["2019-08-23", "2019-08-29"]]
   
    # PHE-defined heatwaves in period 2010-2019, with the ones before 2016 backdated by me (see find_coldsnap_dates.py)
    css_prev = [["2010-01-01", "2010-01-15"], ["2010-01-19", "2010-01-20"], ["2010-01-24", "2010-01-27"], ["2010-01-29", "2010-02-04"], ["2010-02-07", "2010-02-24"], ["2010-02-27", "2010-03-05"], ["2010-03-07", "2010-03-08"], ["2010-11-24", "2010-12-09"], ["2010-12-12", "2010-12-28"], \
                ["2011-01-02", "2011-01-09"], ["2011-01-19", "2011-01-23"], ["2011-01-27", "2011-02-01"], ["2011-02-19", "2011-02-21"], ["2011-12-05", "2011-12-06"], ["2011-12-15", "2011-12-19"], \
                ["2012-01-14", "2012-01-17"], ["2012-01-28", "2012-02-12"], ["2012-02-19", "2012-02-20"], ["2012-11-29", "2012-12-03"], ["2012-12-05", "2012-12-07"], ["2012-12-11", "2012-12-14"], \
                ["2013-01-10", "2013-01-25"], ["2013-02-05", "2013-02-08"], ["2013-02-10", "2013-02-13"], ["2013-02-19", "2013-02-25"], ["2013-03-09", "2013-03-14"], ["2013-03-18", "2013-03-31"], ["2013-11-19", "2013-11-20"], \
                ["2014-12-03", "2014-12-04"], ["2014-12-11", "2014-12-13"], ["2014-12-26", "2014-12-30"], \
                ["2015-01-16", "2015-01-24"], ["2015-01-29", "2015-02-08"], ["2015-02-10", "2015-02-11"], \
                ["2016-02-13", "2016-02-16"], ["2016-02-24", "2016-02-28"], ["2016-12-27", "2016-12-30"], \
                ["2017-01-04", "2017-01-06"], ["2017-01-11", "2017-01-14"], ["2017-01-19", "2017-01-28"], ["2017-02-08", "2017-02-12"], ["2017-11-29", "2017-12-02"], ["2017-12-07", "2017-12-17"], ["2017-12-27", "2017-12-29"], \
                ["2018-01-06", "2018-01-12"], ["2018-01-15", "2018-01-21"], ["2018-02-04", "2018-02-08"], ["2018-02-23", "2018-03-04"], ["2018-03-16", "2018-03-20"], \
                ["2019-01-21", "2019-01-25"], ["2019-01-28", "2019-02-03"]]

    # covid and temperature-related deaths, 2020-01-30 to 2022-12-31 
    covid_hw_ds = []    # heatwaves
    temp_hw_ds = []
    temp_hw_lw = []
    temp_hw_up = []
    covid_cs_ds = []    # cold snaps
    temp_cs_ds = []
    temp_cs_lw = []
    temp_cs_up = []   

    # temperature-related deaths, 2010-01-01 to 2019-12-31
    temp_hw_prev = []
    temp_hw_plow = []
    temp_hw_pup = []
    temp_cs_prev = []
    temp_cs_plow = []
    temp_cs_pup = []
 
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
                              reg+".csv")
        temp_df = temp_df.loc[(temp_df["date"]>="2010-01-01") & (temp_df["date"]<="2022-12-31")]
        # deal with 95% confidence intervals
        temp_df["q2pt5"] = temp_df.loc[:, temp_df.columns.str.startswith("sim")].quantile(0.025, axis=1)
        temp_df["q97pt5"] = temp_df.loc[:, temp_df.columns.str.startswith("sim")].quantile(0.975, axis=1)
        temp_df = temp_df.loc[:, ~temp_df.columns.str.startswith("sim")]
        # isolate non-optimal temperature deaths
        temp_df.loc[temp_df.tmean == temp_df.mmt, "est"] = 0     # best estimate
        temp_df.loc[temp_df.tmean == temp_df.mmt, "q2pt5"] = 0
        temp_df.loc[temp_df.tmean == temp_df.mmt, "q97pt5"] = 0 
 
        # only heatwave days in study period
        covid_hwd = 0
        temp_hwd = 0
        temp_hwl = 0
        temp_hwu = 0
        for hw in hws:
            # subset heatwave dates
            covid_hw_df = covid_df.loc[(covid_df["date"]>=hw[0]) & (covid_df["date"]<=hw[1])]
            covid_hwd += int(covid_hw_df["newDailyNsoDeathsByDeathDate"].sum())
            temp_hw_df = temp_df.loc[(temp_df["date"]>=hw[0]) & (temp_df["date"]<=hw[1])]
            temp_hwd += int(temp_hw_df["est"].sum())
            temp_hwl += int(temp_hw_df["q2pt5"].sum())
            temp_hwu += int(temp_hw_df["q97pt5"].sum())
        covid_hw_ds.append(covid_hwd)
        temp_hw_ds.append(temp_hwd)    
        temp_hw_lw.append(temp_hwl)
        temp_hw_up.append(temp_hwu)

        # only coldsnaps in study period
        covid_csd = 0
        temp_csd = 0
        temp_csl = 0
        temp_csu = 0
        for cs in css:
            # subset coldsnap dates
            covid_cs_df = covid_df.loc[(covid_df["date"]>=cs[0]) & (covid_df["date"]<=cs[1])]
            covid_csd += int(covid_cs_df["newDailyNsoDeathsByDeathDate"].sum())
            temp_cs_df = temp_df.loc[(temp_df["date"]>=cs[0]) & (temp_df["date"]<=cs[1])]
            temp_csd += int(temp_cs_df["est"].sum())
            temp_csl += int(temp_cs_df["q2pt5"].sum())
            temp_csu += int(temp_cs_df["q97pt5"].sum())
        covid_cs_ds.append(covid_csd)     # checked with Wales value here vs covid spreahsheet manually!
        temp_cs_ds.append(temp_csd)
        temp_cs_lw.append(temp_csl)
        temp_cs_up.append(temp_csu)   
 
        # only heatwave days in previous 10 years
        temp_hwd2 = 0
        temp_hwl2 = 0
        temp_hwu2 = 0
        hwdays = 0
        for hw in hws_prev:
            # subset heatwave dates
            temp_hw_df = temp_df.loc[(temp_df["date"]>=hw[0]) & (temp_df["date"]<=hw[1])]
            hwdays += len(temp_hw_df["date"])
            temp_hwd2 += int(temp_hw_df["est"].sum())
            temp_hwl2 += int(temp_hw_df["q2pt5"].sum())
            temp_hwu2 += int(temp_hw_df["q97pt5"].sum())
        temp_hw_prev.append(temp_hwd2)
        temp_hw_plow.append(temp_hwl2)
        temp_hw_pup.append(temp_hwu2)

        # only coldsnaps in previous 10 years
        temp_csd2 = 0
        temp_csl2 = 0
        temp_csu2 = 0
        csdays = 0
        for cs in css_prev:
            # subset coldsnap dates
            temp_cs_df = temp_df.loc[(temp_df["date"]>=cs[0]) & (temp_df["date"]<=cs[1])]
            csdays += len(temp_cs_df["date"])
            temp_csd2 += int(temp_cs_df["est"].sum())
            temp_csl2 += int(temp_cs_df["q2pt5"].sum())
            temp_csu2 += int(temp_cs_df["q97pt5"].sum())
        temp_cs_prev.append(temp_csd2)
        temp_cs_plow.append(temp_csl2)
        temp_cs_pup.append(temp_csu2)
 
    # put data in shapes dataframe
    # sums (in magnitudes) in study period
    eng_wales_df["hw_tempcovid_death_sum_per100k"] = ((np.array(temp_hw_ds)+np.array(covid_hw_ds))/np.array(reg_pops))*100000
    eng_wales_df["hw_tempcovid_low_sum_per100k"] = ((np.array(temp_hw_lw)+np.array(covid_hw_ds))/np.array(reg_pops))*100000
    eng_wales_df["hw_tempcovid_hgh_sum_per100k"] = ((np.array(temp_hw_up)+np.array(covid_hw_ds))/np.array(reg_pops))*100000

    eng_wales_df["cs_tempcovid_death_sum_per100k"] = ((np.array(temp_cs_ds)+np.array(covid_cs_ds))/np.array(reg_pops))*100000    
    eng_wales_df["cs_tempcovid_low_sum_per100k"] = ((np.array(temp_cs_lw)+np.array(covid_cs_ds))/np.array(reg_pops))*100000
    eng_wales_df["cs_tempcovid_hgh_sum_per100k"] = ((np.array(temp_cs_up)+np.array(covid_cs_ds))/np.array(reg_pops))*100000
    
    # sums (in magnitudes) in previous decade
    eng_wales_df["hw_prev_death_70days_per100k"] = (np.array(temp_hw_prev)/(hwdays*np.array(reg_pops_15)))*(70*100000) # ave heatwave deaths per 100 pop per day * 70 days
    eng_wales_df["hw_prev_low_70days_per100k"] = (np.array(temp_hw_plow)/(hwdays*np.array(reg_pops_15)))*(70*100000)
    eng_wales_df["hw_prev_hgh_70days_per100k"] = (np.array(temp_hw_pup)/(hwdays*np.array(reg_pops_15)))*(70*100000)

    eng_wales_df["cs_prev_death_70days_per100k"] = (np.array(temp_cs_prev)/(csdays*np.array(reg_pops_15)))*(70*100000)
    eng_wales_df["cs_prev_low_70days_per100k"] = (np.array(temp_cs_plow)/(csdays*np.array(reg_pops_15)))*(70*100000)
    eng_wales_df["cs_prev_hgh_70days_per100k"] = (np.array(temp_cs_pup)/(csdays*np.array(reg_pops_15)))*(70*100000)

    # ratios of sums (in magnitudes: study period over previous decade
    eng_wales_df["hw_tempcovid_prev_death_ratio"] = eng_wales_df["hw_tempcovid_death_sum_per100k"]/eng_wales_df["hw_prev_death_70days_per100k"]
    eng_wales_df["cs_tempcovid_prev_death_ratio"] = eng_wales_df["cs_tempcovid_death_sum_per100k"]/eng_wales_df["cs_prev_death_70days_per100k"]
    
    # plot regions on map     
    font = {'size' : 14}
    plt.rc('font', **font)

    # set discrete colourmap for sums
    csum = mpl.cm.magma_r(np.linspace(0.1, 0.8, 15))
    cmap2 = mcolors.LinearSegmentedColormap.from_list('mycolormap2', csum)
    norm2 = mpl.colors.BoundaryNorm(np.arange(0,140,10), cmap2.N)
    
    # set discrete colourmap for sums ratios
    crat = mpl.cm.viridis_r(np.linspace(0.1, 0.7, 11))
    cmap3 = mcolors.LinearSegmentedColormap.from_list('mycolormap2', crat)
    norm3 = mpl.colors.BoundaryNorm(np.linspace(1.6, 3.4, 10), cmap3.N)
    
    f, axs = plt.subplots(3, 2, figsize=(10,14))
 
    f.subplots_adjust(bottom=0, top=0.95, left=-0.05, right=0.9,
                    wspace=0, hspace=0.07) 
    
    # heatwaves in period magnitude
    eng_wales_df.plot(column="hw_tempcovid_death_sum_per100k", ax=axs[0,0], edgecolor="black", cmap=cmap2, \
                      norm=norm2, cax=axs[0,0], legend=False)        
    axs[0,0].title.set_text(r"$\bf{a}$"+"  Study period\nHeatwaves (70 days)")   
 
    # cold snaps in period magnitude
    eng_wales_df.plot(column="cs_tempcovid_death_sum_per100k", ax=axs[0,1], edgecolor="black", cmap=cmap2, \
                      norm=norm2, cax=axs[0,1], legend=False)
    axs[0,1].title.set_text(r"$\bf{b}$"+"  Study period\nCold snaps (70 days)")   

    # heatwaves in previous 10 years
    eng_wales_df.plot(column="hw_prev_death_70days_per100k", ax=axs[1,0], edgecolor="black", cmap=cmap2, \
                      norm=norm2, cax=axs[1,0], legend=False)
    axs[1,0].title.set_text(r"$\bf{c}$"+"  Period 2010-2019\nHeatwaves (70 days)")

    # cold snaps in previous 10 years
    eng_wales_df.plot(column="cs_prev_death_70days_per100k", ax=axs[1,1], edgecolor="black", cmap=cmap2, \
                      norm=norm2, cax=axs[1,1], legend=False)
    axs[1,1].title.set_text(r"$\bf{d}$"+"  Period 2010-2019\nCold snaps (70 days)")
    
    # heatwaves in period magnitude over heatwaves in previous 10 years magnitude
    eng_wales_df.plot(column="hw_tempcovid_prev_death_ratio", ax=axs[2,0], edgecolor="black", cmap=cmap3, \
                      norm=norm3, cax=axs[2,0], legend=False)
    axs[2,0].title.set_text(r"$\bf{e}$"+"  Study period vs 2010-2019\nHeatwaves (70 days)")
    
    # cold snaps in period magnitude over cold snaps in previous 10 years magnitude
    eng_wales_df.plot(column="cs_tempcovid_prev_death_ratio", ax=axs[2,1], edgecolor="black", cmap=cmap3, \
                      norm=norm3, cax=axs[2,1], legend=False)
    axs[2,1].title.set_text(r"$\bf{f}$"+"  Study period vs 2010-2019\nCold snaps (70 days)")
    
    # add colourbar for sum
    col = axs[0,1].collections[0]
    cb_ax = f.add_axes([0.88, 0.68, 0.025, 0.25])        # x0, y0, width, height
    cb = f.colorbar(col, cax=cb_ax, fraction=0.05)
    cb.ax.set_ylabel("Temperature plus COVID deaths\nper 100,000 population")

    # add colourbar for previous 10 years
    cb_ax2 = f.add_axes([0.88, 0.355, 0.025, 0.25])
    cb2 = f.colorbar(col, cax=cb_ax2, fraction=0.05)
    cb2.ax.set_ylabel("Temperature deaths\nper 100,000 population")
    
    # add colourbar for study period vs previous 10 years
    col3 = axs[2,1].collections[0]
    cb_ax3 = f.add_axes([0.88, 0.03, 0.025, 0.25])
    cb3 = f.colorbar(col3, cax=cb_ax3, fraction=0.05)
    cb3.ax.set_ylabel("Study to previous period deaths ratio")
    
    for ax in axs.flatten():
        ax.axis("off")

    plt.savefig(out_dir+"EngWales_tempcovid_mortality_sumper100k_202001-202212_all_heatwaves_coldwaves_officaldef_prev_ONSdata_3rows_pop15.pdf", format="pdf", dpi=300)
    plt.close("all")

    print("Saved graph!")