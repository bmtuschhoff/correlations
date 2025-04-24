Below are brief descriptions of the code. These scripts are used to analyze the simulated data and generate figures.

File names that contain the phrase "HeatMap" were used to generate all heat map figures, and those that contain the phrase "TimeSeries" were used to generate all time series figures.

SIRdynams_TimeSeries.R was used to generate the figure for the SIR dynamics and probability of a major epidemic with different levels of heterogeneity and correlations.

ProbPeakSizePeakTimeFES_HeatMap.R was used to generate the figures for the probability of a major epidemic and the peak size, peak time, and final epidemic size given a major epidemic.

TimetoXpercent_HeatMap.R was used to generate the figures for the time at which the epidemic reaches 50% of its final size.

TimeJthInfReff_TimeSeries_HeatMap.R was used to generate the figures for the time to the jth infection, time to the 50th infection, and R effective vs the number of of susceptible individuals.

mpox.R was used to generate the figures for the daily case counts of mpox in New York City from May 19, 2022 to March 8, 2025 and the SEIR dynamics simulated from our model with positive, zero, or negative correlations. The mpox case counts were downloaded on March 17, 2025 from GitHub (https://github.com/nychealth/monkeypox-data) and the New York City mpox information website (https://www.nyc.gov/site/doh/health/health-topics/monkeypox.page).