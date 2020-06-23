#Libraries
library(readr)
library(ggplot2)
library(lmtest) #langrage multiplier test
library(splm)
library(naniar)
library(gplots)
library(spdep)
library(maptools)
library(sf)
library(splm) #spatial panel data analysis package
library(tmap)
library(leaflet)
library(dplyr)
library(corrplot)

#Read the 2 databases-make sure that these files are in your working directory 
mydata <- read.csv('mydata.csv')
mydata_robust2 <- read.csv('mydata_robust2.csv')

#check missing variables plot

gg_miss_var(mydata) #no missing variables in our analysis
head(mydata)

#Data defined as panel data by determining time and spatial IDs
pdata<-pdata.frame(mydata,index=c("Geoid","Time"))
pdata_rob2<-pdata.frame(mydata_robust2,index=c("Geoid","Time"))
#Plot means and check heterogeneity

plotmeans(Crashes ~ Geoid, main="Heterogeineity across census tracts", data=mydata) 
plotmeans(Crashes ~ Time, main="Heterogeineity across month-year units", data=mydata) #useful to track differences in trajectories

#Descriptive stats of the panel data
summary(pdata)

#Spatial models testing 

#check the whole Travis County region of the base analysis
auf.poly <-readShapePoly("auf.shp", IDvar="Geoid")
str(slot(auf.poly,"data"))
summary(auf.poly@data)
leaflet(auf.poly) %>%
  addPolygons(stroke = FALSE, fillOpacity = 0.5, smoothFactor = 0.5) %>%
  addTiles()

#check the region of analysis for the robustness check without suburban extenal census tracts
auf2.poly <-readShapePoly("auf_rob2Final.shp", IDvar="Geoid")
plot(auf2.poly)
str(slot(auf2.poly,"data"))
summary(auf2.poly@data)
leaflet(auf2.poly) %>%
  addPolygons(stroke = FALSE, fillOpacity = 0.5, smoothFactor = 0.5) %>%
  addTiles()

require(RColorBrewer)

#Figure Maps included in Kontou and McDonald Paper comparing rates of crashes and ridesourcing trips before and after the RA introduction
facets = c("Crashes_B", "Crashes_A")
breaks = c(0, 0.5, 1, 2, 5, 10) 
tm_shape(auf.poly) + tm_polygons(facets, breaks=breaks, palette = "YlOrBr", free.scales=FALSE) + 
  tm_facets(nrow = 2, ncol=1, sync = TRUE) +
  tm_layout(panel.labels = c("Before: Avg. Crashes/1K people ", "After: Avg.Crashes/1K people"),  legend.position=c("left", "bottom")) 
tmap_save(filename = "CrashesRateBreaks.png", dpi = 300)
#Figure Maps Ridesourcing
facets = c("TripsRA_B","TripsRA_A")
breaks=c(0,0.05,0.1,0.3,5)
tm_shape(auf.poly) + tm_polygons(facets, breaks=breaks, palette = "YlOrBr") + 
  tm_facets(nrow = 2, ncol=1, sync = TRUE) +
  tm_layout(panel.labels = c("Before: Avg. RA Trips/Population", "After: Avg. RA Trips/Population" ), legend.position=c("left", "bottom"), legend.stack="horizontal") 
tmap_save(filename = "RA_Br.png", dpi = 300)


#Definition of spatial weights
#First order Queen contiguity 
nb.FOQ <- poly2nb(auf.poly, queen=TRUE, row.names=auf.poly$Geoid)
#row.names refers to the unique names of each polygon
nb.FOQ
#First order rook contiguity 
nb.RK <- poly2nb(auf.poly, queen=FALSE, row.names=auf.poly$Geoid)
nb.RK
#Create lists for contiguity

W.FOQ<-nb2listw(nb.FOQ, style="W", zero.policy=NULL)
W.FOQ
W.RK <- nb2listw(nb.RK, style="W", zero.policy=NULL)
#weight matrix based on distances NOT NEEDED HERE 
coords<-coordinates(auf.poly)
W_dist<-dnearneigh(coords,0,1,longlat = FALSE)

#Robustness Checks Contiguity: exclusing outside, suburban census tracts
#First order Queen contiguity 
nb2.FOQ <- poly2nb(auf2.poly, queen=TRUE, row.names=auf2.poly$Geoid)
#row.names refers to the unique names of each polygon
nb2.FOQ
#First order rook contiguity 
nb2.RK <- poly2nb(auf2.poly, queen=FALSE, row.names=auf2.poly$Geoid)
nb2.RK
#Create lists for contiguity

W.FOQ2<-nb2listw(nb2.FOQ, style="W", zero.policy=NULL)
W.FOQ2
W.RK2 <- nb2listw(nb2.RK, style="W", zero.policy=NULL)
#weight matrix based on distances NOT NEEDED HERE 
coords<-coordinates(auf2.poly)
W_dist<-dnearneigh(coords,0,1,longlat = FALSE)


##Spatial Analysis Panel - Crashes
#If p<0.001 reject the hypothesis of absence of correlation between individual effects and explanatory variables
fmc<-LogCrashes ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin

#SARAR for both spatial terms integration
sarar.c1<-spml(fmc, data = pdata, listw = W.FOQ2, model = "within", lag = TRUE, spatial.error = "b", effect = "twoways")
summary(sarar.c1)
sink("C_sarar.txt")
summary(sarar.c1)
sink()

#Spatial error
summary(sc_fe_e<-spml(fmc, pdata, listw=W.FOQ,
                      model=c("within"),
                      effect=c("twoways"),
                      lag=FALSE, spatial.error=c("b")))
sink("C_error.txt")
summary(sc_fe_e)
sink()
#spatial lag 
C_lag<-summary(sc_fe_lag<-spml(fmc, pdata, listw=W.FOQ,
                               model=c("within"),
                               effect=c("twoways"),
                               lag=TRUE, spatial.error=c("none")))
sink("C_lag.txt")
summary(sc_fe_lag)
sink()
#Tests
print(hausman_panel_c<-phtest(fmc, data = pdata))
#Hauman test robust to spatial autocorrelation, if p significant chose fixed effect specification and not random effects
print(spat_hausman_c<-sphtest(fmc,data=pdata,
                              listw =W.FOQ, spatial.model = "error", method="ML"))
#Check if the following tests confirm the rejection of the spatial lag and spatial error terms being nu;;
#Test for spatial lag dependence
c_lml<-slmtest(fmc, data=pdata, listw =W.FOQ, test="lml",model="within", effect='twoways')
#Test for spatial error dependence
c_lme<-slmtest(fmc, data=pdata, listw =W.FOQ, test="lme",model="within", effect='twoways')
#Check for spatial lag dependence sub spatial error
c_rlml<-slmtest(fmc, data=pdata, listw =W.FOQ, test="rlml",model="within", effect='twoways')
#Spatial error dependence sub spatial lag
c_rlme<-slmtest(fmc, data=pdata, listw =W.FOQ, test="rlme",model="within", effect='twoways')
#spatial lag model higher statistic

sink("C_tests")
hausman_panel_c
spat_hausman_c
c_lml
c_lme
c_rlml
c_rlme
sink()

#Injuries
fmi<-LogInjuries ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin

#SARAR for both spatial terms integration
sarar.i<-spml(fmi, data = pdata, listw = W.FOQ, model = "within", lag = TRUE, spatial.error = "b", effect = "twoways")
summary(sarar.i)
sink("I_sarar.txt")
summary(sarar.i)
sink()

#Spatial error
summary(si_fe_e<-spml(fmi, pdata, listw=W.FOQ,
                      model=c("within"),
                      effect=c("twoways"),
                      lag=FALSE, spatial.error=c("b")))
sink("I_error.txt")
summary(si_fe_e)
sink()
#spatial lag 
I_lag<-summary(si_fe_lag<-spml(fmi, pdata, listw=W.FOQ,
                               model=c("within"),
                               effect=c("twoways"),
                               lag=TRUE, spatial.error=c("none")))
sink("I_lag.txt")
summary(si_fe_lag)
sink()
#Tests
print(hausman_panel_i<-phtest(fmi, data = pdata))
#Hauman test robust to spatial autocorrelation, if p significant chose fixed effect specification and not random effects
print(spat_hausman_i<-sphtest(fmi,data=pdata,
                              listw =W.FOQ, spatial.model = "error", method="ML"))
#Check if the following tests confirm the rejection of the spatial lag and spatial error terms being nu;;
#Test for spatial lag dependence
i_lml<-slmtest(fmi, data=pdata, listw =W.FOQ, test="lml",model="within", effect='twoways')
#Test for spatial error dependence
i_lme<-slmtest(fmi, data=pdata, listw =W.FOQ, test="lme",model="within", effect='twoways')
#Check for spatial lag dependence sub spatial error
i_rlml<-slmtest(fmi, data=pdata, listw =W.FOQ, test="rlml",model="within", effect='twoways')
#Spatial error dependence sub spatial lag
i_rlme<-slmtest(fmi, data=pdata, listw =W.FOQ, test="rlme",model="within", effect='twoways')
#spatial lag model higher statistic

sink("I_tests")
hausman_panel_i
spat_hausman_i
i_lml
i_lme
i_rlml
i_rlme
sink()

#Fatalities 
fmf<-LogFatalities ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin

#SARAR for both spatial terms integration
sarar.f<-spml(fmf, data = pdata, listw = W.FOQ, model = "within", lag = TRUE, spatial.error = "b", effect = "twoways")
summary(sarar.f)
sink("F_sarar.txt")
summary(sarar.f)
sink()

#Spatial error
summary(sf_fe_e<-spml(fmf, pdata, listw=W.FOQ,
                      model=c("within"),
                      effect=c("twoways"),
                      lag=FALSE, spatial.error=c("b")))
sink("F_error.txt")
summary(sf_fe_e)
sink()
#spatial lag 
F_lag<-summary(sf_fe_lag<-spml(fmf, pdata, listw=W.FOQ,
                               model=c("within"),
                               effect=c("twoways"),
                               lag=TRUE, spatial.error=c("none")))
sink("F_lag.txt")
summary(sf_fe_lag)
sink()
#Tests
print(hausman_panel_f<-phtest(fmf, data = pdata))
#Hauman test robust to spatial autocorrelation, if p significant chose fixed effect specification and not random effects
print(spat_hausman_f<-sphtest(fmf,data=pdata,
                              listw =W.FOQ, spatial.model = "error", method="ML"))
#Check if the following tests confirm the rejection of the spatial lag and spatial error terms being nu;;
#Test for spatial lag dependence
f_lml<-slmtest(fmf, data=pdata, listw =W.FOQ, test="lml",model="within", effect='twoways')
#Test for spatial error dependence
f_lme<-slmtest(fmf, data=pdata, listw =W.FOQ, test="lme",model="within", effect='twoways')
#Check for spatial lag dependence sub spatial error
f_rlml<-slmtest(fmf, data=pdata, listw =W.FOQ, test="rlml",model="within", effect='twoways')
#Spatial error dependence sub spatial lag
f_rlme<-slmtest(fmf, data=pdata, listw =W.FOQ, test="rlme",model="within", effect='twoways')
#spatial lag model higher statistic

sink("F_tests")
hausman_panel_f
spat_hausman_f
f_lml
f_lme
f_rlml
f_rlme
sink()

#DWIs
fmd<-LogDWI~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin

#SARAR for both spatial terms integration
sarar.d<-spml(fmd, data = pdata, listw = W.FOQ, model = "within", lag = TRUE, spatial.error = "b", effect = "twoways")
summary(sarar.d)
sink("D_sarar.txt")
summary(sarar.d)
sink()

#Spatial error
summary(sd_fe_e<-spml(fmd, pdata, listw=W.FOQ,
                      model=c("within"),
                      effect=c("twoways"),
                      lag=FALSE, spatial.error=c("b")))
sink("D_error.txt")
summary(sd_fe_e)
sink()
#spatial lag 
#file_C_lag <- tempfile(pattern = "table_C_lag") # or 
D_lag<-summary(sd_fe_lag<-spml(fmd, pdata, listw=W.FOQ,
                               model=c("within"),
                               effect=c("twoways"),
                               lag=TRUE, spatial.error=c("none")))
sink("D_lag.txt")
summary(sd_fe_lag)
sink()
#Tests
print(hausman_panel_d<-phtest(fmd, data = pdata))
#Hauman test robust to spatial autocorrelation, if p significant chose fixed effect specification and not random effects
print(spat_hausman_d<-sphtest(fmd,data=pdata,
                              listw =W.FOQ, spatial.model = "error", method="ML"))
#Check if the following tests confirm the rejection of the spatial lag and spatial error terms being nu;;
#Test for spatial lag dependence
d_lml<-slmtest(fmd, data=pdata, listw =W.FOQ, test="lml",model="within", effect='twoways')
#Test for spatial error dependence
d_lme<-slmtest(fmd, data=pdata, listw =W.FOQ, test="lme",model="within", effect='twoways')
#Check for spatial lag dependence sub spatial error
d_rlml<-slmtest(fmd, data=pdata, listw =W.FOQ, test="rlml",model="within", effect='twoways')
#Spatial error dependence sub spatial lag
d_rlme<-slmtest(fmd, data=pdata, listw =W.FOQ, test="rlme",model="within", effect='twoways')
#spatial lag model higher statistic

sink("D_tests")
hausman_panel_d
spat_hausman_d
d_lml
d_lme
d_rlml
d_rlme
sink()


#Robust 2 data (excluding external suburban census tracts from the analysis)
#If p<0.001 reject the hypothesis of absence of correlation between individual effects and explanatory variables
fmc<-LogCrashes ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin

#SARAR for both spatial terms integration
sarar.c1<-spml(fmc, data = pdata_rob2, listw = W.FOQ2, model = "within", lag = TRUE, spatial.error = "b", effect = "twoways")
summary(sarar.c1)
sink("C_sarar_rob2.txt")
summary(sarar.c1)
sink()

#Spatial error
summary(sc_fe_e<-spml(fmc, pdata_rob2, listw=W.FOQ2,
                    model=c("within"),
                    effect=c("twoways"),
                    lag=FALSE, spatial.error=c("b")))
sink("C_error_rob2.txt")
summary(sc_fe_e)
sink()
#spatial lag 
#file_C_lag <- tempfile(pattern = "table_C_lag") # or 
C_lag<-summary(sc_fe_lag<-spml(fmc, pdata_rob2, listw=W.FOQ2,
                    model=c("within"),
                    effect=c("twoways"),
                    lag=TRUE, spatial.error=c("none")))
#the method not working with splm package
#table2spreadsheet(x=C_lag,file=file_C_lag, digits = 3, digitspvals = 4, type=c("csv"))
#write.table(C_lag, file = "C_lag.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
sink("C_lag_rob2.txt")
summary(sc_fe_lag)
sink()
#Tests
print(hausman_panel_c<-phtest(fmc, data = pdata_rob2))
#Hauman test robust to spatial autocorrelation, if p significant chose fixed effect specification and not random effects
print(spat_hausman_c<-sphtest(fmc,data=pdata_rob2,
                              listw =W.FOQ2, spatial.model = "error", method="ML"))
#Check if the following tests confirm the rejection of the spatial lag and spatial error terms being nu;;
#Test for spatial lag dependence
c_lml<-slmtest(fmc, data=pdata_rob2, listw =W.FOQ2, test="lml",model="within", effect='twoways')
#Test for spatial error dependence
c_lme<-slmtest(fmc, data=pdata_rob2, listw =W.FOQ2, test="lme",model="within", effect='twoways')
#Check for spatial lag dependence sub spatial error
c_rlml<-slmtest(fmc, data=pdata_rob2, listw =W.FOQ2, test="rlml",model="within", effect='twoways')
#Spatial error dependence sub spatial lag
c_rlme<-slmtest(fmc, data=pdata_rob2, listw =W.FOQ2, test="rlme",model="within", effect='twoways')
#spatial lag model higher statistic

sink("C_tests_rob2")
hausman_panel_c
spat_hausman_c
c_lml
c_lme
c_rlml
c_rlme
sink()

#Injuries
fmi<-LogInjuries ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin

#SARAR for both spatial terms integration
sarar.i<-spml(fmi, data = pdata_rob2, listw = W.FOQ2, model = "within", lag = TRUE, spatial.error = "b", effect = "twoways")
summary(sarar.i)
sink("I_sarar_rob2.txt")
summary(sarar.i)
sink()

#Spatial error
summary(si_fe_e<-spml(fmi, pdata_rob2, listw=W.FOQ2,
                      model=c("within"),
                      effect=c("twoways"),
                      lag=FALSE, spatial.error=c("b")))
sink("I_error_rob2.txt")
summary(si_fe_e)
sink()
#spatial lag 
#file_C_lag <- tempfile(pattern = "table_C_lag") # or 
I_lag<-summary(si_fe_lag<-spml(fmi, pdata_rob2, listw=W.FOQ2,
                               model=c("within"),
                               effect=c("twoways"),
                               lag=TRUE, spatial.error=c("none")))
#the method not working with splm package
#table2spreadsheet(x=C_lag,file=file_C_lag, digits = 3, digitspvals = 4, type=c("csv"))
#write.table(C_lag, file = "C_lag.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
sink("I_lag_rob2.txt")
summary(si_fe_lag)
sink()
#Tests
print(hausman_panel_i<-phtest(fmi, data = pdata_rob2))
#Hauman test robust to spatial autocorrelation, if p significant chose fixed effect specification and not random effects
print(spat_hausman_i<-sphtest(fmi,data=pdata_rob2,
                              listw =W.FOQ2, spatial.model = "error", method="ML"))
#Check if the following tests confirm the rejection of the spatial lag and spatial error terms being nu;;
#Test for spatial lag dependence
i_lml<-slmtest(fmi, data=pdata_rob2, listw =W.FOQ2, test="lml",model="within", effect='twoways')
#Test for spatial error dependence
i_lme<-slmtest(fmi, data=pdata_rob2, listw =W.FOQ2, test="lme",model="within", effect='twoways')
#Check for spatial lag dependence sub spatial error
i_rlml<-slmtest(fmi, data=pdata_rob2, listw =W.FOQ2, test="rlml",model="within", effect='twoways')
#Spatial error dependence sub spatial lag
i_rlme<-slmtest(fmi, data=pdata_rob2, listw =W.FOQ2, test="rlme",model="within", effect='twoways')
#spatial lag model higher statistic

sink("I_tests_rob2")
hausman_panel_i
spat_hausman_i
i_lml
i_lme
i_rlml
i_rlme
sink()

#Fatalities 
fmf<-LogFatalities ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin

#SARAR for both spatial terms integration
sarar.f<-spml(fmf, data = pdata_rob2, listw = W.FOQ2, model = "within", lag = TRUE, spatial.error = "b", effect = "twoways")
summary(sarar.f)
sink("F_sarar_rob2.txt")
summary(sarar.f)
sink()

#Spatial error
summary(sf_fe_e<-spml(fmf, pdata_rob2, listw=W.FOQ2,
                      model=c("within"),
                      effect=c("twoways"),
                      lag=FALSE, spatial.error=c("b")))
sink("F_error_rob2.txt")
summary(sf_fe_e)
sink()
#spatial lag 
#file_C_lag <- tempfile(pattern = "table_C_lag") # or 
F_lag<-summary(sf_fe_lag<-spml(fmf, pdata_rob2, listw=W.FOQ2,
                               model=c("within"),
                               effect=c("twoways"),
                               lag=TRUE, spatial.error=c("none")))
#the method not working with splm package
#table2spreadsheet(x=C_lag,file=file_C_lag, digits = 3, digitspvals = 4, type=c("csv"))
#write.table(C_lag, file = "C_lag.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
sink("F_lag_rob2.txt")
summary(sf_fe_lag)
sink()
#Tests
print(hausman_panel_f<-phtest(fmf, data_rob2 = pdata))
#Hauman test robust to spatial autocorrelation, if p significant chose fixed effect specification and not random effects
print(spat_hausman_f<-sphtest(fmf,data=pdata_rob2,
                              listw =W.FOQ2, spatial.model = "error", method="ML"))
#Check if the following tests confirm the rejection of the spatial lag and spatial error terms being nu;;
#Test for spatial lag dependence
f_lml<-slmtest(fmf, data=pdata_rob2, listw =W.FOQ2, test="lml",model="within", effect='twoways')
#Test for spatial error dependence
f_lme<-slmtest(fmf, data=pdata_rob2, listw =W.FOQ2, test="lme",model="within", effect='twoways')
#Check for spatial lag dependence sub spatial error
f_rlml<-slmtest(fmf, data=pdata_rob2, listw =W.FOQ2, test="rlml",model="within", effect='twoways')
#Spatial error dependence sub spatial lag
f_rlme<-slmtest(fmf, data=pdata_rob2, listw =W.FOQ2, test="rlme",model="within", effect='twoways')
#spatial lag model higher statistic

sink("F_tests_rob2")
hausman_panel_f
spat_hausman_f
f_lml
f_lme
f_rlml
f_rlme
sink()

#DWIs
fmd<-LogDWI~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin

#SARAR for both spatial terms integration
sarar.d<-spml(fmd, data = pdata_rob2, listw = W.FOQ2, model = "within", lag = TRUE, spatial.error = "b", effect = "twoways")
summary(sarar.d)
sink("D_sarar_rob2.txt")
summary(sarar.d)
sink()

#Spatial error
summary(sd_fe_e<-spml(fmd, pdata_rob2, listw=W.FOQ2,
                      model=c("within"),
                      effect=c("twoways"),
                      lag=FALSE, spatial.error=c("b")))
sink("D_error_rob2.txt")
summary(sd_fe_e)
sink()
#spatial lag 
#file_C_lag <- tempfile(pattern = "table_C_lag") # or 
D_lag<-summary(sd_fe_lag<-spml(fmd, pdata_rob2, listw=W.FOQ2,
                               model=c("within"),
                               effect=c("twoways"),
                               lag=TRUE, spatial.error=c("none")))
#the method not working with splm package
#table2spreadsheet(x=C_lag,file=file_C_lag, digits = 3, digitspvals = 4, type=c("csv"))
#write.table(C_lag, file = "C_lag.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
sink("D_lag_rob2.txt")
summary(sd_fe_lag)
sink()
#Tests
print(hausman_panel_d<-phtest(fmd, data = pdata_rob2))
#Hauman test robust to spatial autocorrelation, if p significant chose fixed effect specification and not random effects
print(spat_hausman_d<-sphtest(fmd,data=pdata_rob2,
                              listw =W.FOQ2, spatial.model = "error", method="ML"))
#Check if the following tests confirm the rejection of the spatial lag and spatial error terms being nu;;
#Test for spatial lag dependence
d_lml<-slmtest(fmd, data=pdata_rob2, listw =W.FOQ2, test="lml",model="within", effect='twoways')
#Test for spatial error dependence
d_lme<-slmtest(fmd, data=pdata_rob2, listw =W.FOQ2, test="lme",model="within", effect='twoways')
#Check for spatial lag dependence sub spatial error
d_rlml<-slmtest(fmd, data=pdata_rob2, listw =W.FOQ2, test="rlml",model="within", effect='twoways')
#Spatial error dependence sub spatial lag
d_rlme<-slmtest(fmd, data=pdata_rob2, listw =W.FOQ2, test="rlme",model="within", effect='twoways')
#spatial lag model higher statistic

sink("D_tests_rob2")
hausman_panel_d
spat_hausman_d
d_lml
d_lme
d_rlml
d_rlme
sink()

#Moran's I testing
crashesnum <- aggregate(Crashes ~ Geoid, mydata, sum)
injuriesnum <- aggregate (Injuries ~ Geoid, mydata, sum)
fatalitiesnum <- aggregate (Fatalities ~ Geoid, mydata, sum)
DWInum <- aggregate (DWI ~ Geoid, mydata, sum)
#Moran I test for cross-sectional data #null hypothesis of spatial randomness
moran.test(crashesnum$Crashes, listw = W.FOQ)
moran.test(injuriesnum$Injuries, listw=W.FOQ)
moran.test(fatalitiesnum$Fatalities, listw=W.FOQ)
moran.test(DWInum$DWI, listw=W.FOQ)
moran.plot(crashesnum$Crashes, listw = W.FOQ) #certain shapes of the cloud affect the selection of the fitted line (the slope of it correspond's to Moran's I statistic)
geary.test(crashesnum$Crashes, listw = W.FOQ) #rejecting the null hypothesis for spatial autocorrelation
plot(moran.mc(crashesnum$Crashes, listw = W.FOQ, nsim=1000))

#Correlation Plots Outcome - included in the paper's appendix
df1<-mydata %>%
select( Crashes, Fatalities, Injuries, DWI, MedHHInc, TotalTrips, EmpPrcnt, MedHHInc, VehOwn0, PopDensity, TripsRideAustin)
#Appendix Correlation Plot
corrplot.mixed(cor(df1), lower = "ellipse", upper = "number", tl.pos = "lt", diag = "u", tl.col="black")