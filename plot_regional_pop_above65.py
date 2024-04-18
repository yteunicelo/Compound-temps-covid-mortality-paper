# Requires conda activate geo_env
import numpy as np
import geopandas as gp
import pandas as pd
from datetime import datetime, timedelta
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt


''' This script plot England and Wales regions '''
''' and show the proportion of population above the age of 65 '''
''' based on 2021 census data, unrounded '''
''' Eunice Lo '''
''' Created 08/03/2024 '''


if __name__ == "__main__":

    # paths
    shp_dir = "ONS/"
    in_dir = "data/"
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
    regnames = ["North East", "North West", "Yorkshire and The Humber", "East Midlands", \
               "West Midlands", "East of England", "London", "South East", "South West", \
               "Wales"]

    # 2021 census population data, unrounded
    # https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/bulletins/populationandhouseholdestimatesenglandandwales/census2021unroundeddata
    pop_df = pd.read_csv(in_dir+"Census2021_unrounded_mar2023release_population_byage_regions.csv", header=0)

    p65 = np.zeros(10)

    for rn in regnames:
        pop_rn = pop_df[pop_df["Regions"]==rn]
        pop_all = pop_rn["Observation"].sum()
        pop_65plus = pop_rn[pop_rn["Age (101 categories) Code"]>65]["Observation"].sum()
        p65[regnames.index(rn)] = (pop_65plus/pop_all)*100

    eng_wales_df["percentage_age65+"] = p65

    '''
    # 2021 census population data, first release 
    # https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationandhouseholdestimatesenglandandwalescensus2021
    skiprows = set(list(np.arange(383))) - set([7, 10, 24, 68, 93, 133, 168, 219, 255, 326, 360])
    pop_df = pd.read_excel(in_dir+"census2021firstresultsenglandwales1.xlsx", sheet_name="P02", \
                           header=0, skiprows=list(skiprows), usecols="A:C,Q:V") 
    
    # percentage of population above 65, in the same order as regnames, checked
    eng_wales_df["percentage_age65+"] = ((pop_df["Aged 65 to 69 years\n[note 12]"] + pop_df["Aged 70 to 74 years\n[note 12]"] +\
                                          pop_df["Aged 75 to 79 years\n[note 12]"] + pop_df["Aged 80 to 84 years\n[note 12]"] +\
                                          pop_df["Aged 85 to 89 years\n[note 12]"] + pop_df["Aged 90 years and over\n[note 12]"])/pop_df["All persons"])*100
    '''
    
    # plot regions on map     
    font = {'size' : 12}
    plt.rc('font', **font)

    # set discrete colourmap for percentages
    cpop = mpl.cm.magma(np.linspace(0.2, 0.95, 11))
    cmap = mcolors.LinearSegmentedColormap.from_list('my_colormap', cpop)
    norm = mpl.colors.BoundaryNorm(np.linspace(11,22,12), cmap.N) 
    
    f, ax = plt.subplots(1, 1, figsize=(4,4))
    
    f.subplots_adjust(bottom=0.05, top=0.9, left=0.05, right=0.8,
                      wspace=0.02, hspace=0.02) 

    eng_wales_df.plot(column="percentage_age65+", ax=ax, edgecolor="black", cmap=cmap, norm=norm, \
                      cax=ax, legend=False)
    ax.title.set_text("Census 2021 population")
    
    # add colourbar
    col = ax.collections[0]
    cb_ax = f.add_axes([0.8, 0.06, 0.03, 0.8])        # x0, y0, width, height
    cb = f.colorbar(col, cax=cb_ax, fraction=0.05)
    cb.ax.set_ylabel("% of population aged above 65")

    ax.axis("off")
 
    plt.savefig(out_dir+"EngWales_percentage_pop_above_age65_census2021unrounded.png", format="png", dpi=300)
    plt.close("all")

    print("Saved graph!")
