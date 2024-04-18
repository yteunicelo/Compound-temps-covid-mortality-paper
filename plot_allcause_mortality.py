# requires source activate pyn_env
from __future__ import division
import numpy as np
import pandas as pd
from datetime import datetime, timedelta, date
from functools import reduce
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.dates as mdates


''' This script plots the time series of all cause mortality '''
''' between 2020-01-30 and 2022-12-31 '''
''' vs the same calendar days but from the 10 years before that '''
''' Eunice Lo '''
''' Created 01/12/2023 '''
''' Updated 21/03/2024 with population weighted temperature results (in SI) '''
''' Updated 17/04/2024 with means for error bars and saving figure as pdf '''


if __name__ == "__main__":

    in_dir = "/home/bridge/yl17544/papers_repos/EnglandWales_climate_vs_covid_paper/"
    ons_dir = "/export/anthropocene/array-01/yl17544/ONS/Daily_deaths_by_region_England_Wales/"
    out_dir = "/home/bridge/yl17544/plots/EnglandWales_climate_vs_covid/"

    regions = ["NorthEastEngland", "NorthWestEngland", "YorkshireandHumber", "EastMidlands", \
                "WestMidlands", "EastofEngland", "London", "SouthEastEngland", "SouthWestEngland", \
                "Wales"] 
     
    # within study period
    regional_allcause_nocovid = []
    regional_covid_mort = []
    regional_heat_mort = []    # each entry is a region. each region has columns: date, est, low, hgh; with low-hgh being 2.5-97.5th percentiles
    regional_cold_mort = []
    regional_heat_prev = []    # this will be a nested list: nregs x 10
    regional_cold_prev = []

    for reg in regions:

        print("======"+reg+"======")
        
        #### within study period
        # all-cause mortality without covid mortality
        allcause = pd.read_csv(in_dir+"data/ONS_HadUK-Grid_2020_2022_without_covidcert_noleap_"+reg+".csv", \
                               usecols=[0,6])
        allcause = allcause.loc[allcause["date"]>="2020-01-30"]
        regional_allcause_nocovid.append(allcause)
        
        # covid mortality (on death certs)
        covid = pd.read_csv(in_dir+"data/covid_deaths_on_cert_data_"+reg+"_2023-Mar-23.csv", usecols=[3,4])
        covid = covid.loc[covid["date"]<="2022-12-31"]
        # reverse in time
        covid.sort_values(by="date", inplace=True)       
        # add dates from 2020-01-30 if doesn't exist
        date1 = datetime.strptime(covid["date"].iloc[0], '%Y-%m-%d').date()
        sdate = datetime.strptime("2020-01-30", '%Y-%m-%d').date()
        if date1 > sdate:
           edate = date1 - timedelta(days=1)    
           extras = pd.date_range(start=sdate, end=edate)        	    # extra dates
           extras2 = [date.strftime('%Y-%m-%d') for date in extras]
           extra_ds = [0] * len(extras2)                                    # extra daily deaths, all 0
           extra_data = {"date":extras2, "newDailyNsoDeathsByDeathDate":extra_ds}
           extra_df = pd.DataFrame(data=extra_data)       
           covid = pd.concat([extra_df, covid])
        # remove leap days
        covid = covid[~covid.date.str.endswith("02-29")]
        regional_covid_mort.append(covid)
        
        # temperature-related mortality
        temp_mort = pd.read_csv(in_dir+"outputs/daily_attributable_deaths_ER1981-2022_yearround_21dayslag_MMT2to98_nsim100_1981-2022_ONSdata_"+reg+".csv") 
        temp_period = temp_mort.loc[temp_mort["date"]>="2020-01-30"]
        # heat est and sims
        heat_mask = temp_period["tmean"] < temp_period["mmt"]    # make heat deaths 0 if too cold
        heat_period = temp_period.copy(deep=True) 
        heat_period.iloc[heat_mask, 3:] = 0
        regional_heat_mort.append(heat_period)
        # cold est and sims
        cold_mask = temp_period["tmean"] > temp_period["mmt"]    # make cold deaths 0 if too hot
        cold_period = temp_period.copy(deep=True)
        cold_period.iloc[cold_mask, 3:] = 0
        regional_cold_mort.append(cold_period) 
           
        #### previous 10 years: 2010-01-01 to 2019-12-31
        sdate = "2010-01-01"
        edate = "2019-12-31"

        # temperature-related mortality
        temp_prev = temp_mort.loc[(temp_mort["date"]>=sdate) & (temp_mort["date"]<=edate)]     
        # heat est and sims
        heat_mask2 = temp_prev["tmean"] < temp_prev["mmt"]    # make heat deaths 0 if too cold
        heat_prev = temp_prev.copy(deep=True)
        heat_prev.iloc[heat_mask2, 3:] = 0
        regional_heat_prev.append(heat_prev)
        # cold est and sims
        cold_mask2 = temp_prev["tmean"] > temp_prev["mmt"]    # make heat deaths 0 if too cold
        cold_prev = temp_prev.copy(deep=True)
        cold_prev.iloc[cold_mask2, 3:] = 0
        regional_cold_prev.append(cold_prev)          
 
    # do some plotting!
    font = {'size' : 12}
    plt.rc('font', **font)

    f, axs = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]}, figsize=(12, 6))

    # England and Wales totals
    allcause_tots = np.zeros(1066)
    covid_tots = np.zeros(1066)
    heat_tots = np.zeros(1066)
    heat_tots_sims = np.zeros((1066,100))
    cold_tots = np.zeros(1066)
    cold_tots_sims = np.zeros((1066,100))

    heat_prev_tots = np.zeros(3650)
    heat_prev_tots_sims = np.zeros((3650,100))
    cold_prev_tots = np.zeros(3650)
    cold_prev_tots_sims = np.zeros((3650,100))

    for nr in range(10):
        allcause_tots += np.array(regional_allcause_nocovid[nr]["death_allage"])
        covid_tots += np.array(regional_covid_mort[nr]["newDailyNsoDeathsByDeathDate"])
        heat_tots += np.array(regional_heat_mort[nr]["est"])
        cold_tots += np.array(regional_cold_mort[nr]["est"])
        
        heat_tots_sims += np.array(regional_heat_mort[nr].iloc[:,4:])
        cold_tots_sims += np.array(regional_cold_mort[nr].iloc[:,4:])

        heat_prev_tots += np.array(regional_heat_prev[nr]["est"])
        cold_prev_tots += np.array(regional_cold_prev[nr]["est"])

        heat_prev_tots_sims += np.array(regional_heat_prev[nr].iloc[:,4:])
        cold_prev_tots_sims += np.array(regional_cold_prev[nr].iloc[:,4:])

    # previous ten years but ONS all-cause, still England and Wales tots
    ons_ds = pd.read_excel(ons_dir+"dailydeaths19812020adhocfinal.xlsx", \
                           sheet_name="Table 2", header=4, skipfooter=0, usecols="A:E,AI:AR", engine="openpyxl")
    # remove leap days
    ons_ds.drop(ons_ds[(ons_ds.Month==2) & (ons_ds.Day==29)].index, inplace=True)
    allcause_prev_tots = ons_ds.groupby(["Month","Day"]).sum() 

    # main plot over the study period
    tots = pd.DataFrame({"date":regional_allcause_nocovid[0]["date"], \
                         "allcause_nocovid":reduce(lambda a, b: a.add(b, fill_value=0), regional_allcause_nocovid)["death_allage"], \
                         "covid":covid_tots, \
                         "heat":np.array(reduce(lambda a, b: a.add(b, fill_value=0), regional_heat_mort)["est"]), \
                         "cold":np.array(reduce(lambda a, b: a.add(b, fill_value=0), regional_cold_mort)["est"])})
    
    dates = [datetime.strptime(d, '%Y-%m-%d') for d in regional_heat_mort[0]["date"].values]

    axs[0].plot(dates, allcause_tots+covid_tots, color="black", linewidth=2, label="All-cause deaths", zorder=1)
    axs[0].plot(dates, covid_tots, color="tab:purple", linewidth=2, label="COVID-19 deaths", zorder=2)
    # 95% eCI!
    axs[0].fill_between(dates, np.percentile(heat_tots_sims,q=2.5,axis=1), np.percentile(heat_tots_sims,q=97.5,axis=1), color="pink", alpha=0.5)
    axs[0].plot(dates, heat_tots, color="tab:red", linewidth=2, label="Heat deaths", zorder=3)
    axs[0].fill_between(dates, np.percentile(cold_tots_sims,q=2.5,axis=1), np.percentile(cold_tots_sims,q=97.5,axis=1), color="lightskyblue", alpha=0.5)
    axs[0].plot(dates, cold_tots, color="tab:blue", linewidth=2, label="Cold deaths", zorder=4)    
    axs[0].legend(loc="upper right", fontsize=12)  

    # covid variants
    axs[0].axvline(datetime.strptime("2020-12-18", '%Y-%m-%d'), color="grey", linestyle=":", zorder=0)
    axs[0].text(datetime.strptime("2020-12-18", '%Y-%m-%d'), 3300, "Alpha", rotation=90)
    axs[0].axvline(datetime.strptime("2021-05-17", '%Y-%m-%d'), color="grey", linestyle=":", zorder=0)
    axs[0].text(datetime.strptime("2021-05-17", '%Y-%m-%d'), 3300, "Delta", rotation=90)
    axs[0].axvline(datetime.strptime("2021-12-10", '%Y-%m-%d'), color="grey", linestyle=":", zorder=0)
    axs[0].text(datetime.strptime("2021-12-10", '%Y-%m-%d'), 3300, "Omicron", rotation=90)    

    # covid vaccination
    axs[0].axvline(datetime.strptime("2021-07-01", '%Y-%m-%d'), color="grey", linestyle="-", alpha=0.7, zorder=0)
    axs[0].text(datetime.strptime("2021-07-01", '%Y-%m-%d'), 3300, "~50% vaccination", rotation=90)

    axs[0].set_title("Study period")
    axs[0].set_title("a", loc='left', fontweight="bold")
 
    # axes
    axs[0].set_ylim(0, 3500)
    axs[0].set_ylabel("Number of deaths in England and Wales")
    axs[0].yaxis.set_minor_locator(MultipleLocator(100))

    axs[0].set_xlim(datetime.strptime("2020-01-01",'%Y-%m-%d'), datetime.strptime("2023-05-01",'%Y-%m-%d')) 
    xt = ["2020-01-01","2020-03-01","2020-05-01","2020-07-01","2020-09-01","2020-11-01", \
          "2021-01-01","2021-03-01","2021-05-01","2021-07-01","2021-09-01","2021-11-01", \
          "2022-01-01","2022-03-01","2022-05-01","2022-07-01","2022-09-01","2022-11-01", \
          "2023-01-01"]
    xticks = [datetime.strptime(t, '%Y-%m-%d') for t in xt]
    axs[0].set_xticks(xticks)
    axs[0].set_xticklabels(xt, rotation=90, ha="center")
    axs[0].set_xlabel("Date")

    # Minor ticks every month
    fmt_month = mdates.MonthLocator()
    axs[0].xaxis.set_minor_locator(fmt_month)

    axs[0].tick_params('both', length=8, width=1, which='major')
    axs[0].tick_params('both', length=4, width=1, which='minor')

    # remove extra xticks on the right
    xticks_min = axs[0].xaxis.get_minor_ticks()
    xticks_min[-1].set_visible(False)
    xticks_min[-2].set_visible(False)
    xticks_min[-3].set_visible(False)
    xticks_min[-4].set_visible(False)
    xticks_min[-5].set_visible(False)

    # previous ten years' range in the same plot
    allcause_ave = np.mean(allcause_prev_tots.values.ravel())
    allcause_ave_low = np.min(allcause_prev_tots.values.ravel())
    allcause_ave_hgh = np.max(allcause_prev_tots.values.ravel())

    heat_ave_est = np.mean(heat_prev_tots)     
    heat_ave_low = np.min(heat_prev_tots)
    heat_ave_hgh = np.max(heat_prev_tots)
    cold_ave_est = np.mean(cold_prev_tots)
    cold_ave_low = np.min(cold_prev_tots)
    cold_ave_hgh = np.max(cold_prev_tots)

    axs[0].errorbar(datetime.strptime("2023-03-15", '%Y-%m-%d'), allcause_ave, yerr=[[allcause_ave-allcause_ave_low], [allcause_ave_hgh-allcause_ave]], \
                    fmt=".", ms=10, color="k", ecolor="k", elinewidth=3, capsize=4) 
    axs[0].errorbar(datetime.strptime("2023-03-01", '%Y-%m-%d'), heat_ave_est, yerr=[[heat_ave_est-heat_ave_low], [heat_ave_hgh-heat_ave_est]], \
                    fmt=".", ms=10, color="tab:red", ecolor="tab:red", elinewidth=3, capsize=4)
    axs[0].errorbar(datetime.strptime("2023-04-01", '%Y-%m-%d'), cold_ave_est, yerr=[[cold_ave_est-cold_ave_low], [cold_ave_hgh-cold_ave_est]], \
                    fmt=".", ms=10, color="tab:blue", ecolor="tab:blue", elinewidth=3, capsize=4)     
 
    # annotate the period
    axs[0].annotate("Period\n2010-2019", (datetime.strptime("2023-04-01", '%Y-%m-%d'), -400), rotation=90, ha="center", \
                    style="italic", annotation_clip=False)

    # second panel
    # zoomed in over a heatwave/cold snap
    axs[1].plot(dates, allcause_tots+covid_tots, color="black", linewidth=2)
    axs[1].plot(dates, covid_tots, color="tab:purple", linewidth=2)
    # 95% eCI!
    axs[1].fill_between(dates, np.percentile(heat_tots_sims,q=2.5,axis=1), np.percentile(heat_tots_sims,q=97.5,axis=1), color="pink", alpha=0.5)
    axs[1].plot(dates, heat_tots, color="tab:red", linewidth=2)
    axs[1].fill_between(dates, np.percentile(cold_tots_sims,q=2.5,axis=1), np.percentile(cold_tots_sims,q=97.5,axis=1), color="lightskyblue", alpha=0.5)
    axs[1].plot(dates, cold_tots, color="tab:blue", linewidth=2)

    # peak temp date
    axs[1].annotate(r"40.3$^\circ$C", xy=(pd.Timestamp("2022-07-19"),700), xytext=(pd.Timestamp("2022-07-19"),950), ha="center", \
                    arrowprops=dict(facecolor="tab:red", shrink=0.05))
    
    axs[1].set_title("July 2022 heatwave")
    axs[1].set_title("b", loc='left', fontweight="bold")
    
    # axes, panel 2
    axs[1].set_ylim(0, 2200)
    axs[1].yaxis.set_major_locator(MultipleLocator(500))
    axs[1].yaxis.set_minor_locator(MultipleLocator(100))

    axs[1].set_xlim(datetime.strptime("2022-07-10",'%Y-%m-%d'), datetime.strptime("2022-07-25",'%Y-%m-%d'))
    xt2 = ["2022-07-10","2022-07-15", "2022-07-20","2022-07-25"]
    xticks2 = [datetime.strptime(t, '%Y-%m-%d') for t in xt2]
    axs[1].set_xticks(xticks2)
    axs[1].set_xticklabels(xt2, rotation=90, ha="center")
    axs[1].set_xlabel("Date")

    # Minor ticks every day
    fmt_day = mdates.DayLocator()
    axs[1].xaxis.set_minor_locator(fmt_day) 

    axs[1].tick_params('both', length=8, width=1, which='major')
    axs[1].tick_params('both', length=4, width=1, which='minor')

    # save
    plt.tight_layout()
    plt.savefig(out_dir+"EngWales_mortality_types_timeseries_202001-202212_prevtenyrs_hwzoom_ONSdata_final.pdf", format="pdf", dpi=300)
    plt.close("all")

    print("Saved graph!") 
