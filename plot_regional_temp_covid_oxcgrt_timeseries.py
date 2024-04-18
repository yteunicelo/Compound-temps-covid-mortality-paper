# requires source activate pyn_env
import numpy as np
import geopandas as gp
import pandas as pd
from datetime import datetime, timedelta
import matplotlib as mpl
mpl.use('Agg')
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.dates as mdates


''' This script plots time series of temp and covid deaths '''
''' in England and Wales regions '''
''' Eunice Lo '''
''' Created 30/03/2023 '''


if __name__ == "__main__":

    # paths
    in_dir = "/home/bridge/yl17544/papers_repos/EnglandWales_climate_vs_covid_paper/outputs/"
    in_dir2 = "/home/bridge/yl17544/papers_repos/EnglandWales_climate_vs_covid_paper/data/"
    out_dir = "/home/bridge/yl17544/plots/EnglandWales_climate_vs_covid/"

    eng_regnames = ["NorthEastEngland" , "NorthWestEngland", "YorkshireandHumber", "EastMidlands", \
                    "WestMidlands", "EastofEngland", "London", "SouthEastEngland", "SouthWestEngland", \
                    "Wales"]
 
    # mid-2021 population by region, copied from 
    # www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland    
    reg_pops = [2646772, 7422295, 5481431, 4880094, 5954240, 6348096, 8796628, 9294023, 5712840, 3105410]

    # set up figure
    font = {'size' : 14}
    plt.rc('font', **font)
    
    f, axs = plt.subplots(5, 2, sharex="col", sharey="row", figsize=(15, 12))

    f.subplots_adjust(bottom=0.1, top=0.96, left=0.08, right=0.99,
                      wspace=0.05, hspace=0.2) 

    # loop regions in England
    for reg, ax in zip(eng_regnames, axs.flatten()):
        # covid deaths on certificate
        covid_df = pd.read_csv(in_dir2+"covid_deaths_on_cert_data_"+reg+"_2023-Mar-23.csv", \
                   usecols=["date", "newDailyNsoDeathsByDeathDate"])
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
           extra_data = {"date":extras2, "newDailyNsoDeathsByDeathDate":extra_ds}
           extra_df = pd.DataFrame(data=extra_data)       
           covid_df = pd.concat([extra_df, covid_df])
        # remove leap days
        covid_df = covid_df[~covid_df.date.str.endswith("02-29")]
        # new daily deaths per 100,000 people
        covid_df["daily_per100k"] = (covid_df["newDailyNsoDeathsByDeathDate"]/reg_pops[eng_regnames.index(reg)])*100000

	# temperature
        temp_df = pd.read_csv(in_dir+"daily_attributable_deaths_ER1981-2022_yearround_21dayslag_MMT2to98_nsim100_1981-2022_"+\
                              reg+".csv", usecols=["date","tmean"])
        # subset to 30/01/2020 to 31/12/2022
        temp_df = temp_df.loc[(temp_df["date"]>="2020-01-30") & (temp_df["date"]<="2022-12-31")]
       
        # Vaccination rate data, note column names are different for England vs Wales
        if reg!="Wales":
            vac_df = pd.read_csv(in_dir2+"covid_vacc_percentage_uptake_"+reg+"_2023-May-11.csv", \
                                 usecols=["date","cumVaccinationSecondDoseUptakeByVaccinationDatePercentage"])
        else:
            vac_df = pd.read_csv(in_dir2+"covid_vacc_percentage_uptake_"+reg+"_2023-May-11.csv", \
                                 usecols=["date","cumVaccinationSecondDoseUptakeByPublishDatePercentage"])
            vac_df.rename(columns={"cumVaccinationSecondDoseUptakeByPublishDatePercentage":"cumVaccinationSecondDoseUptakeByVaccinationDatePercentage"}, \
                          inplace=True)
        # subset to 30/01/2020 to 31/12/2022
        vac_df = vac_df.loc[(vac_df["date"]>="2020-01-30") & (vac_df["date"]<="2022-12-31")]
        # reverse in time
        vac_df.sort_values(by="date", inplace=True)
        # below 0.1% are NaN in data. Treat as 0.
        vac_df.fillna(0, inplace=True)
        # add dates from 2020-01-30 if doesn't exist
        date1v = datetime.strptime(vac_df["date"].iloc[0], '%Y-%m-%d').date()
        if date1v > sdate:
           edatev = date1v - timedelta(days=1)
           extrasv = pd.date_range(start=sdate, end=edatev)                   # extra dates
           extras2v = [date.strftime('%Y-%m-%d') for date in extrasv]
           extra_dsv = [0] * len(extras2v)                                    # extra daily deaths, all 0
           extra_datav = {"date":extras2v, "cumVaccinationSecondDoseUptakeByVaccinationDatePercentage":extra_dsv}
           extra_dfv = pd.DataFrame(data=extra_datav)
           vac_df = pd.concat([extra_dfv, vac_df])
        # remove leap days
        vac_df = vac_df[~vac_df.date.str.endswith("02-29")]   

        '''     
        # Oxford Covid Gov Response Tracker (OxCGRT)
        if reg != "Wales":
             ox_df = pd.read_csv(in_dir2+"OxCGRT_GBR_latest.csv", usecols=["RegionName","Date","GovernmentResponseIndex_Average"]).loc[oxcgrt_df["RegionName"]=="England"]
        else:
             ox_df = pd.read_csv(in_dir2+"OxCGRT_GBR_latest.csv", usecols=["RegionName","Date","GovernmentResponseIndex_Average"]).loc[oxcgrt_df["RegionName"]=="Wales"]
        # subset to 30/01/2020 to 31/12/2022
        ox_df["date2"] = ox_df["Date"].apply(lambda x: pd.to_datetime(str(x), format='%Y%m%d')).astype(str)
        ox_df = ox_df.loc[(ox_df["date2"]>="2020-01-30") & (ox_df["date2"]<="2022-12-31")] 
        # remove leap days
        ox_df = ox_df[~ox_df.date2.str.endswith("02-29")]
        '''
        # plot time series
        covid_df.plot(x="date", y="daily_per100k", ax=ax, title=reg, color="black", linewidth=2, ylabel="", x_compat=True, xlabel="", legend=False, rot=30, zorder=10)
        vac_df.plot(x="date", y="cumVaccinationSecondDoseUptakeByVaccinationDatePercentage", ax=ax, secondary_y=True, color="tab:blue", linewidth=2, ylabel="", x_compat=True, \
                    xlabel="", legend=False, rot=30, zorder=5)
        #ox_df.plot(x="date2", y="GovernmentResponseIndex_Average", ax=ax, secondary_y=True, color="tab:blue", linewidth=2, ylabel="", x_compat=True, \
        #           xlabel="", legend=False, rot=30, zorder=5)   

        # third y axis
        ax3 = ax.twinx()
        rspine = ax3.spines['right']
        rspine.set_position(('axes', 1.15))
        ax3.set_frame_on(True)
        ax3.patch.set_visible(False)
        ax3.spines["right"].set_color("tab:red")
        ax3.tick_params(axis="y", colors="tab:red")

        temp_df.plot(x="date", y="tmean", ax=ax3, color="tab:red", linewidth=1, ylim=(0,25), ylabel="", x_compat=True, xlabel="", legend=False, rot=30, zorder=0)

        # set first 2 y axes and ylim
        ax.set_ylim(0,4)
        ax.right_ax.set_ylim(0,100)
        ax.right_ax.spines["right"].set_color("tab:blue")
        ax.right_ax.tick_params(axis="y", colors="tab:blue")
     
        # set x, y ticks
        ax.tick_params('both', length=8, width=1, which='major')
        ax.right_ax.tick_params('both', length=8, width=1, which='major')
        ax3.tick_params('both', length=8, width=1, which='major')     
 
    # common x and y axes
    f.text(0.5, 0.01, "Date", ha='center', fontsize=16)
    f.text(0.01, 0.5, "Daily COVID (certificate) deaths per 100,000 population", va='center', rotation='vertical', fontsize=16)
    f.text(0.94, 0.5, "Percentage of population aged 12+ vaccinated with 2 doses (%)", color="tab:blue", va="center", rotation=270, fontsize=16)
    f.text(0.96, 0.5, "Daily mean temperature ($^{\circ}$C)", color="tab:red", va="center", rotation=270, fontsize=16)    

    '''
    # set legend
    colours = ["black", "tab:blue"]
    lines = [Line2D([0], [0], color=c, linewidth=2) for c in colours]
    labels = ["COVID deaths", "GRI"]
    axs[4,0].legend(lines, labels, ncol=2, loc="upper left") 
    '''
    # save
    plt.tight_layout(rect=[0.04,0.01,0.96,1])
    plt.savefig(out_dir+"EngWales_mortality_dailyper100k_timeseries_vaccine2dose_202001-202212.png", format="png", dpi=300)
    plt.close("all")

    print("Saved graph!")     
