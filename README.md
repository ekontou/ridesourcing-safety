# Ridesourcing & Road Safety Associations Empirical Research

[![DOI](https://zenodo.org/badge/197442332.svg)](https://zenodo.org/badge/latestdoi/197442332)

This repo is hosting R code for the ridesourcing use - road safety outcomes spatial modeling portion of the Kontou and McDonald paper (https://arxiv.org/abs/2001.03461). You can find the list of libraries needed to run the spatial models in the install.R file. 

The data csv file is an join outcome of several datasets described hereinafter. Detailed historical road safety measures stem from the Texas Department of Transportation Crash Records Information System ([Texas Department of Transportation 2019](https://www.txdot.gov/government/enforcement/crash-statistics.html)) and DWI offenses from the Austin Police Department Crime Reports ([Austin Police Department 2019](https://data.austintexas.gov/Public-Safety/Crime-Reports-2018/vmn9-3bvu)). For both road safety and DWI offense records, apart from temporal variability, we have access to their exact locations as longitude and latitude coordinates. Ridesourcing comprehensive trip-level data from RideAustin’s open record (see [here](https://data.world/ride-austin)) are leveraged. We aim to control for overall traffic in the Travis county region, after conducting monthly travel demand analyses using the StreetLight Data platform (StreetLight Data 2019). Normalized trip counts for all census tracts within Travis County over the period of interest are used. The vehicular normalized traffic counts data stem from the [StreetLight Data](https://www.streetlightdata.com/?utm_source=Google-Adwords&utm_medium=Paid-Search&utm_campaign=StreetLight-Brand&utm_term=%2Bstreetlight%20%2Bdata&creative=372389731084&keyword=%2Bstreetlight%20%2Bdata&matchtype=b&network=g&device=c&utm_term=%2Bstreetlight%20%2Bdata&utm_campaign=StreetLight-Data-Brand&utm_source=adwords&utm_medium=ppc&hsa_acc=7146595976&hsa_cam=1079169723&hsa_grp=51692947334&hsa_ad=372389731084&hsa_src=g&hsa_tgt=kwd-419587122414&hsa_kw=%2Bstreetlight%20%2Bdata&hsa_mt=b&hsa_net=adwords&hsa_ver=3&gclid=EAIaIQobChMI-MyD5Pq76gIVhYbACh2y7wZtEAAYASAAEgJCWPD_BwE) platform.


