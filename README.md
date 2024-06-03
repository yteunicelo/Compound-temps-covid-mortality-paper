# Compound-temps-covid-mortality-paper
Code associated with Lo et al., 2024, Compound mortality impacts from extreme temperatures and the COVID-19 pandemic. Nature Communications. https://www.nature.com/articles/s41467-024-48207-2. 

data/ contains COVID-19 mortality and vaccination rate that were open source (links in manuscript) but are no longer available online. Other open souce data that are still available online are not included here, but their links are listed in the manuscript.

ONS/ contains the adminstrative shapfiles for this analysis.

All temperature-mortality modelling is done in 01_02_04_05_06_MortalitySeriesUncert_heatcold.r, with open source HadUK-Grid and all-cause mortality data linked in the manuscript. This script also plots supplementary Figure S3.

Figure 1 is plotted in plot_allcause_mortality.py

Figure 2 is plotted in plot_regional_temp_covid_deaths_ratio.py

Figure 3 is plotted in plot_regional_temp_covid_deaths_mag.py

Figure 4 is plotted in plot_EnglandWales_tempdeaths_temp_episodes.py

Supplementary figures:

Figure S1 is plotted in plot_regional_temp_covid_oxcgrt_timeseries.py

Figure S2 is plotted in plot_regional_pop_above65.py
