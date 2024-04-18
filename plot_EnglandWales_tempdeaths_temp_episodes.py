# requires source activate pyn_env
from __future__ import division
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt


''' This script plots scatter plots '''
''' of temp deaths vs mean temp in heatwave/coldsnap episodes '''
''' across England '''
''' indicating COVID years '''
''' Eunice Lo '''
''' Created 18/04/2023 '''

if __name__ == "__main__":

    # paths
    in_dir = "/home/bridge/yl17544/papers_repos/EnglandWales_climate_vs_covid_paper/outputs/"
    out_dir = "/home/bridge/yl17544/plots/EnglandWales_climate_vs_covid/"
 
    # UKHSA-defined heatwaves in period 2016 to 2022, according to PhE reports. Dates are inclusive
    hws = [["2016-07-18", "2016-07-22"], ["2016-08-22", "2016-08-26"], ["2016-09-12", "2016-09-17"], \
           ["2017-06-16", "2017-06-23"], ["2017-07-05", "2017-07-07"], \
           ["2018-06-25", "2018-06-27"], ["2018-06-30", "2018-07-10"], ["2018-07-21", "2018-07-29"], ["2018-08-01", "2018-08-09"], \
           ["2019-06-28", "2019-06-30"], ["2019-07-21", "2019-07-28"], ["2019-08-23", "2019-08-29"], \
           ["2020-06-23", "2020-06-27"], ["2020-07-30", "2020-08-01"], ["2020-08-05", "2020-08-15"], \
           ["2021-07-16", "2021-07-23"], ["2021-09-06", "2021-09-09"], \
           ["2022-06-16", "2022-06-19"], ["2022-07-10", "2022-07-25"], ["2022-07-30", "2022-08-05"], ["2022-08-08", "2022-08-17"], ["2022-08-23", "2022-08-25"]]

    # cold waves in period
    # defined here as dates on which a L3 Cold Weather Alert was issued for any region in England (see CWP England)
    # excludes single-day alerts
    css = [["2016-02-13", "2016-02-16"], ["2016-02-24", "2016-02-28"], ["2016-12-27", "2016-12-30"], \
           ["2017-01-04", "2017-01-06"], ["2017-01-11", "2017-01-14"], ["2017-01-19", "2017-01-28"], ["2017-02-08", "2017-02-12"], ["2017-11-29", "2017-12-02"], ["2017-12-07", "2017-12-17"], \
           ["2018-01-06", "2018-01-12"], ["2018-01-15", "2018-01-21"], ["2018-02-04", "2018-02-08"], ["2018-02-23", "2018-03-04"], ["2018-03-16", "2018-03-20"], \
           ["2019-01-21", "2019-01-25"], ["2019-01-28", "2019-02-03"], \
           ["2020-12-29", "2021-01-18"], \
           ["2021-01-22", "2021-02-02"], ["2021-02-08", "2021-02-12"], ["2021-11-26", "2021-11-29"], ["2021-12-20", "2021-12-23"], \
           ["2022-01-04", "2022-01-10"], ["2022-01-13", "2022-01-17"], ["2022-12-07", "2022-12-18"]]

    # England and Wales mortality
    regnames = ["NorthEastEngland", "NorthWestEngland", "YorkshireandHumber", "EastMidlands", \
                "WestMidlands", "EastofEngland", "London", "SouthEastEngland", "SouthWestEngland", \
                "Wales"]

    temp_dflist = []

    for reg in regnames:

        temp_df= pd.read_csv(in_dir+"daily_attributable_deaths_ER1981-2022_yearround_21dayslag_MMT2to98_nsim100_1981-2022_"+\
                             reg+".csv", usecols=["date","tmean","mmt","est"])
        # isolate non-optimal temperature deaths
        temp_df.loc[temp_df.tmean == temp_df.mmt, "est"] = 0     # best estimate

        # put in list
        temp_dflist.append(temp_df.loc[:,("date", "tmean", "est")])

    # mean Tmean across England and Wales
    engwales_df = temp_dflist[0].copy()
    engwales_df["tmean"] = sum([df["tmean"] for df in temp_dflist])/len(temp_dflist)

    # summed est deaths across England and Wales
    engwales_df["est"] = sum([df["est"] for df in temp_dflist])

    # lists for storing heatwave results
    hws_avetemp = []
    hws_culdeaths = []

    for hw in hws:
        # isolate heatwave
        hw_df = engwales_df.loc[(engwales_df["date"]>=hw[0]) & (engwales_df["date"]<=hw[1])]
        # heatwave mean temp
        hws_avetemp.append(hw_df["tmean"].mean())
        # heatwave deaths
        hws_culdeaths.append(hw_df["est"].mean())
       
    # lists for storing coldsnaps results
    css_avetemp = []
    css_culdeaths = []

    for cs in css:
        # isolate heatwave
        cs_df = engwales_df.loc[(engwales_df["date"]>=cs[0]) & (engwales_df["date"]<=cs[1])]
        # coldsnap mean temp
        css_avetemp.append(cs_df["tmean"].mean())
        # coldsnap deaths
        css_culdeaths.append(cs_df["est"].mean())

    # plot graph
    font = {'size' : 10}
    plt.rc('font', **font)

    f, axs = plt.subplots(1, 2, figsize=(10,5))

    # heatwaves
    # colour code by year, grey for 2016-2019, other colours for 2020, 2021, 2022 (done manually)
    cmaph = ["grey"]*12 + ["pink"]*3 + ["red"]*2 + ["darkred"]*5
    axs[0].scatter(hws_avetemp, hws_culdeaths, c=cmaph, marker="d", s=50, alpha=0.9)
    axs[0].set_title("Heatwaves")
    axs[0].set_ylabel("Average temperature-related deaths (/day)")
    axs[0].set_xlabel("Average temperature ($^{\circ}$C)")
    # legend
    hleg = [Line2D([0],[0], marker='d', label='2016-19', color="w", markerfacecolor='grey', markersize=10), \
            Line2D([0],[0], marker='d', label='2020', color="w", markerfacecolor='pink', markersize=10), \
            Line2D([0],[0], marker='d', label='2021', color="w", markerfacecolor='red', markersize=10), \
            Line2D([0],[0], marker='d', label='2022', color="w", markerfacecolor='darkred', markersize=10)]
    axs[0].legend(handles=hleg, loc='upper left')    

    # coldsnaps
    # colour code by year
    cmapc = ["grey"]*16 + ["lightskyblue"]*1 + ["blue"]*4 + ["darkblue"]*3    
    axs[1].scatter(css_avetemp, css_culdeaths, c=cmapc, marker="d", s=50, alpha=0.9)
    axs[1].set_title("Cold snaps")
    axs[1].set_ylabel("Average temperature-related deaths (/day)")
    axs[1].set_xlabel("Average temperature ($^{\circ}$C)")
    # legend
    cleg = [Line2D([0],[0], marker='d', label='2016-19', color="w", markerfacecolor='grey', markersize=10), \
            Line2D([0],[0], marker='d', label='2020', color="w", markerfacecolor='lightskyblue', markersize=10), \
            Line2D([0],[0], marker='d', label='2021', color="w", markerfacecolor='blue', markersize=10), \
            Line2D([0],[0], marker='d', label='2022', color="w", markerfacecolor='darkblue', markersize=10)]
    axs[1].legend(handles=cleg, loc=0)

    f.suptitle("England and Wales-wide results")

    # save graph
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(out_dir+"EngWales_avetempmortality_vs_avetemp_all_heatwaves_coldwaves_officaldef_2016_2022.png", format="png", dpi=300)
    plt.close("all")

    print("Saved graph!") 
