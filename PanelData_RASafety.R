library(plm)
library(splm)
library(readr)
library(clubSandwich)
library(ggplot2)
library(naniar)
library(AER)
library(MASS)
library(lmtest)
library(pglm)
library(readstata13)

#Read Data
mydata <- read.csv("Documents/CrashDataAustin/ODs_all/Rinputs_Jun152019_v3.csv")
attach(mydata)
gg_miss_var(mydata) #check missing variables plot
head(mydata)

#Set data as panel data by determining time and spatial ids
pdata<-pdata.frame(mydata,index=c("Geoid","Time"))
library(gplots)

plotmeans(Crashes ~ Geoid, main="Heterogeineity across census tracts", data=mydata) 
plotmeans(Crashes ~ Time, main="Heterogeineity across month-year units", data=mydata) #useful to track differences in trajectories

#Descriptve stats
summary(pdata)
coplot(Crashes ~ Time|Geoid==484530000000) #plots crashes over time for a given census tract

##Pooled OLS estimator
pooling_C<-plm(LogCrashes ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin, data=pdata, model='pooling') #crashes
summary(pooling_C)
plmtest(pooling_C, effect = "twoways", type = "ghm")
pooling_Inj<-plm(LogInjuries ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin, data=pdata, model='pooling') #injuries 
summary(pooling_Inj)
pooling_Fat<-plm(LogFatalities ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin, data=pdata, model='pooling') #fatalities 
summary(pooling_Fat)
pooling_DWI<-plm(LogDWI ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin, data=pdata, model='pooling') #DWI
summary(pooling_DWI)

# Fixed effects estimator
#Crashes
wi_Crashes_ind <- plm(LogCrashes ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin, data=pdata, model= 'within',effect='individual')
summary(wi_Crashes_ind)
wi_Crashes <- plm(LogCrashes ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+TripsRideAustin, data=pdata, model= 'within',effect='twoways')
summary(wi_Crashes)
coef_test(wi_Crashes, vcov = "CR1", cluster = "individual", test = "naive-t")
coef_test(wi_Crashes, vcov = "CR2", cluster = "individual", test = "Satterthwaite")
coef_test(wi_Crashes, vcov = "CR2", cluster = "individual", test = "Satterthwaite")
coeftest(wi_Crashes, vcov=vcovHC(wi_Crashes, cluster="group"))
qqnorm(residuals(wi_Crashes), ylab = 'Residuals')
qqline(residuals(wi_Crashes))
hist(residuals(wi_Crashes), xlab = 'Residuals')
fixef(wi_Crashes)

#Fatalities
wi_Fatalities <- plm(LogFatalities ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin, data=pdata, model= 'within',effect='twoways')
summary(wi_Fatalities)
coef_test(wi_Fatalities, vcov = "CR1", cluster = "individual", test = "naive-t")
coef_test(wi_Fatalities, vcov = "CR2", cluster = "individual", test = "Satterthwaite")
coeftest(wi_Fatalities, vcov=vcovHC(wi_Fatalities, cluster="group"))
qqnorm(residuals(wi_Fatalities), ylab = 'Residuals')
qqline(residuals(wi_Fatalities))
hist(residuals(wi_Fatalities), xlab = 'Residuals')

#Injuries
wi_Injuries <- plm(LogInjuries ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin , data=pdata, model= 'within',effect='twoways')
summary(wi_Injuries)
coef_test(wi_Injuries, vcov = "CR1", cluster = "individual", test = "naive-t")
coef_test(wi_Injuries, vcov = "CR2", cluster = "individual", test = "Satterthwaite")
coeftest(wi_Injuries, vcov=vcovHC(wi_Injuries, cluster="group"))
qqnorm(residuals(wi_Injuries), ylab = 'Residuals')
qqline(residuals(wi_Injuries))
hist(residuals(wi_Injuries), xlab = 'Residuals')

#DWIs offences
wi_DWI <- plm(LogDWI ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin, data=pdata, model= 'within',effect='twoways')
summary(wi_DWI)
coef_test(wi_DWI, vcov = "CR1", cluster = "individual", test = "naive-t")
coef_test(wi_DWI, vcov = "CR2", cluster = "individual", test = "Satterthwaite")
coeftest(wi_DWI, vcov=vcovHC(wi_DWI, cluster="group"))
qqnorm(residuals(wi_DWI), ylab = 'Residuals')
qqline(residuals(wi_DWI))
hist(residuals(wi_DWI), xlab = 'Residuals')


## Random effects estimator
r_Crashes <- plm(LogCrashes ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin, data=pdata, model= 'random', effect='individual')
summary(r_Crashes)
r_Injuries <- plm(LogInjuries ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin, data=pdata, model= 'random', effect='individual')
summary(r_Injuries)
r_Fatalities <- plm(LogFatalities ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin, data=pdata, model= 'random', effect='individual')
summary(r_Injuries)
r_DWI <- plm(LogDWI ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin, data=pdata, model= 'random', effect='individual')
summary(r_DWI)

#Poisson model for count panel data
#Crashes
#dispersiontest(poisson_c) #test only for linear poisson
poisson_c <- pglm(Crashes ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+TripsRideAustin, pdata,
                  family = poisson, model = "within", effect='twoways', index=c("Geoid","Time"))
summary(poisson_c)
standard_se<-ginv(-poisson_c$hessian)
coeftest(poisson_c,standard_se)
#By exponeting the coefficients we get a ratio of sample means
exp(coef(poisson_c))

#Injuries
#dispersiontest(poisson_i)
poisson_i <- pglm(Injuries ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+TripsRideAustin, pdata,
                  family = poisson, model = "within", effect='twoways', index=c("Geoid","Time"))
#summary(poisson_i)
standard_se<-ginv(-poisson_i$hessian)
coeftest(poisson_i,standard_se)
#By exponeting the coefficients we get a ratio of sample means
exp(coef(poisson_i))

#Fatalities
#dispersiontest(poisson_f)
poisson_f <- pglm(Fatalities ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+TripsRideAustin, pdata,
                  family = poisson, model = "within", effect='twoways', index=c("Geoid","Time"))
summary(poisson_f)
standard_se<-ginv(-poisson_f$hessian)
coeftest(poisson_f,standard_se)
#By exponeting the coefficients we get a ratio of sample means
exp(coef(poisson_f))

#DWI offences
#dispersiontest(poisson_d)
poisson_d <- pglm(DWI ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+TripsRideAustin, pdata,
                  family = poisson, model = "within", effect='twoways', index=c("Geoid","Time"))
summary(poisson_d)
standard_se<-ginv(-poisson_d$hessian)
coeftest(poisson_d,standard_se)


# LM test for fixed effects versus OLS (null hypothesis OLS better than fixed effects) if pvalue<0.01 reject null
pFtest(wi_Crashes, pooling_C)
pFtest(wi_Injuries, pooling_Inj)
pFtest(wi_Fatalities, pooling_Fat)
pFtest(wi_DWI, pooling_DWI)

# Hausman test for fixed versus random effects model
#null hypothesis error not correlated with regressors, if supported use random effect
#p<0.05 thus use fixed effects
phtest(r_Crashes, wi_Crashes)
phtest(r_Injuries, wi_Injuries)
phtest(r_Fatalities, wi_Fatalities)
phtest(r_DWI, wi_DWI)

#No Controls - Only RideSourcing Exposure Variable
wi_Crashes_nc <- plm(LogCrashes ~ LogTripsRideAustin, data=pdata, model= 'within',effect='twoways')
summary(wi_Crashes_nc)
coef_test(wi_Crashes_nc, vcov = "CR1", cluster = "individual", test = "naive-t")
coef_test(wi_Crashes_nc, vcov = "CR2", cluster = "individual", test = "Satterthwaite")


wi_Fatalities_nc <- plm(LogFatalities ~ LogTripsRideAustin, data=pdata, model= 'within',effect='twoways')
summary(wi_Fatalities_nc)
coef_test(wi_Fatalities_nc, vcov = "CR1", cluster = "individual", test = "naive-t")
coef_test(wi_Fatalities_nc, vcov = "CR2", cluster = "individual", test = "Satterthwaite")

wi_Injuries_nc <- plm(LogInjuries ~ LogTripsRideAustin, data=pdata, model= 'within',effect='twoways')
summary(wi_Injuries_nc)
coef_test(wi_Injuries_nc, vcov = "CR1", cluster = "individual", test = "naive-t")
coef_test(wi_Injuries_nc, vcov = "CR2", cluster = "individual", test = "Satterthwaite")

wi_DWI_nc <- plm(LogDWI ~ LogTripsRideAustin, data=pdata, model= 'within',effect='twoways')
summary(wi_DWI_nc)
coeftest(wi_DWI_nc, vcov=vcovHC(wi_DWI, cluster="group"))
coef_test(wi_DWI_nc, vcov = "CR1", cluster = "individual", test = "naive-t")
coef_test(wi_DWI_nc, vcov = "CR2", cluster = "individual", test = "Satterthwaite")

#Spatial models testing 
library(spdep)
library(maptools)
library(sf)


plot(auf.poly)
str(slot(auf.poly,"data"))
summary(auf.poly@data$AreaSqFt2)
auf.poly <-readShapePoly("Documents/CrashDataAustin/ODs_all/auf.shp", IDvar="Geoid")
plot(auf.poly)
str(slot(auf.poly,"data"))
summary(auf.poly@data)
library(leaflet)
leaflet(auf.poly) %>%
  addPolygons(stroke = FALSE, fillOpacity = 0.5, smoothFactor = 0.5) %>%
  addTiles()

require(RColorBrewer)
## Loading required package: RColorBrewer
#Crashes
Npal <- colorNumeric(palette = "OrRd", domain = auf.poly@data$Crashes_B)
leaflet(auf.poly) %>%
  addPolygons(stroke = FALSE, fillOpacity = .9, smoothFactor = 0.9, color = ~Npal(Crashes_B)) %>%
  addTiles() %>%
  addLegend("topright", pal = Npal, values = ~Crashes_A, title = "Avg. Monthly Crashes Before",opacity = 0.9)
Npala <- colorNumeric(palette = "OrRd", domain = auf.poly@data$Crashes_A)
leaflet(auf.poly) %>%
  addPolygons(stroke = FALSE, fillOpacity = .9, smoothFactor = 0.9, color = ~Npal(Crashes_A)) %>%
  addTiles() %>%
  addLegend("topright", pal = Npala, values = ~Crashes_A, title = "Avg. Monthly Crashes After",opacity = 0.9)
#Inj
Npal <- colorNumeric(palette = "OrRd", domain = auf.poly@data$Injuries_B)
leaflet(auf.poly) %>%
  addPolygons(stroke = FALSE, fillOpacity = .9, smoothFactor = 0.9, color = ~Npal(Injuries_B)) %>%
  addTiles() %>%
  addLegend("topright", pal = Npal, values = ~Injuries_B, title = "Avg. Monthly Injuries Before",opacity = 0.9)
Npala <- colorNumeric(palette = "OrRd", domain = auf.poly@data$Inj_A)
leaflet(auf.poly) %>%
  addPolygons(stroke = FALSE, fillOpacity = .9, smoothFactor = 0.9, color = ~Npal(Injuries_A)) %>%
  addTiles() %>%
  addLegend("topright", pal = Npala, values = ~Injuries_A, title = "Avg. Monthly Injuries After",opacity = 0.9)

#Figure Maps Crashes in the paper
facets = c("Crashes_B", "Crashes_A" )
breaks = c(0, 0.5, 1, 2, 5, 10) 
tm_shape(auf.poly) + tm_polygons(facets, breaks=breaks, palette = "YlOrBr", free.scales=FALSE) + 
  tm_facets(nrow = 2, ncol=1, sync = TRUE) +
  tm_layout(panel.labels = c("Before: Avg. Crashes/1K people ", "After: Avg.Crashes/1K people"),  legend.position=c("left", "bottom")) 
tmap_save(filename = "CrashesRateBreaks.png", dpi = 300)
#Figure Maps Ridesourcing
facets = c("TripsRA_B","TripsRA_A" )
breaks=c(0,0.05,0.1,0.3,5)
tm_shape(auf.poly) + tm_polygons(facets, breaks=breaks, palette = "YlOrBr") + 
  tm_facets(nrow = 2, ncol=1, sync = TRUE) +
  tm_layout(panel.labels = c("Before: Avg. RA Trips/Population", "After: Avg. RA Trips/Population" ), legend.position=c("left", "bottom"), legend.stack="horizontal") 
tmap_save(filename = "RA_Br.png", dpi = 300)


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

#OLS
library(spdep)
fm<-LogCrashes ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin
Austin.ols <- lm(fm, data = auf.poly@data)
#Moran's I test for spatial autoregression
matrix.FOQ <- nb2mat(nb.FOQ, style="W")

fixed_s <- spml(fm,data=pdata,listw=mat2listw(matrix.FOQ), model = "within", spatial.error = "b")
summary(fixed_s)
re <- spreml(fm, data = pdata , listw = W.FOQ, effects = "spfe", method = "eigen")
summary(re)
#Random effects 
re.mod<-spreml(fm, data=pdata, w=matrix.FOQ, error="re")
summary(re.mod)
#Random effects and spatial dependence 
semre.mod<-spreml(fm, data=pdata, w=matrix.FOQ, error="semre")
summary(semre.mod)
#Sarar Model combining both lag and spatial error
sarar.mod<-spml(fm, data = pdata, listw = W.FOQ, model = "within", lag = TRUE, spatial.error = "b", effect = "individual")
summary(sarar.mod)
#sarfemod Spatial panel fixed effects error model
sarfemod<-spml(fm, data = pdata, index = NULL, listw = W.FOQ, model="within", effect="individual")
summary(sarfemod)
eff <- effects(sarfemod)
#semfemod Spatial panel fixed effects error model
semfedmod<-spml(fm, data = pdata, index = NULL, listw = W.FOQ, model="within", effect="time")
summary(semfedmod)
eff <- effects(semfedmod)


#Crashes
##Tests
#Hausman test
print(hausman_panel<-phtest(fm, data = pdata)) # result suggest fixed effects
# Hausman test robust to spatial autocorrelation (splm)
print(spat_hausman_ML_SEM<-sphtest(fm,data=pdata,listw =W.FOQ, spatial.model = "error", method="ML"))
#Hausman test for spatial models
print(spat_hausman_ML_SAR<-sphtest(fm,data=pdata,listw =W.FOQ,spatial.model = "lag", method="ML"))

#Hausman test robust to spatial autocorrelation of errors leads to rejection of the null hypothesis on absence of correlation between individual effects and explanatory variables. For the rest of the empirical analysis, a fixed effects model is thus chosen.

# Fixed effects model
# Test 1
slmtest(fm, data=pdata, listw = W.FOQ, test="lml",model="within")

slmtest(fm, data=pdata, listw = W.FOQ, test="lme",model="within")

#confirm the rejection of the hypothesis that these two terms (taken independently) are null
#Robust tests
slmtest(fm, data=pdata, listw = W.FOQ, test="rlml", model="within")
slmtest(fm, data=pdata, listw = W.FOQ, test="rlme",model="within")

# Likelihood Maximum estimation
summary(crashes_SEM_pool <- spml(fm, data=pdata, listw = W.FOQ, lag=TRUE,model="pooling"))
# Fixed-effect SEM
summary(crashes_SEM_FE<- spml(fm, data=pdata, listw = W.FOQ, lag=TRUE,model="within", effect="individual", spatial.error="b"))
summary(crashes_SEM_FE<- spml(fm, data=pdata, listw = W.FOQ, lag=TRUE,model="within", effect="individual", spatial.error="kkp"))
# Generalised moments method estimation
summary(crashes_SEM_FE_GM <- spgm(fm, data=pdata, listw = W.FOQ, model="within", moments="fullweights",spatial.error = TRUE))

#Injuries
fm_In<-LogInjuries ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin
##Tests
#Hausman test
print(hausman_panel_In<-phtest(fm_In, data = pdata)) # result suggest fixed effects

# Hausman test robust to spatial autocorrelation (splm)
print(spat_hausman_ML_SEM_In<-sphtest(fm_In,data=pdata,listw =W.FOQ, spatial.model = "error", method="ML"))
#Hausman test for spatial models
print(spat_hausman_ML_SAR_In<-sphtest(fm_In,data=pdata,listw =W.FOQ,spatial.model = "lag", method="ML"))

#Hausman test robust to spatial autocorrelation of errors leads to rejection of the null hypothesis on absence of correlation between individual effects and explanatory variables. For the rest of the empirical analysis, a fixed effects model is thus chosen.

# Fixed effects model
# Test 1
slmtest(fm_In, data=pdata, listw = W.FOQ, test="lml",model="within")

slmtest(fm_In, data=pdata, listw = W.FOQ, test="lme",model="within")

#confirm the rejection of the hypothesis that these two terms (taken independently) are null
#Robust tests
slmtest(fm_In, data=pdata, listw = W.FOQ, test="rlml", model="within")
slmtest(fm_In, data=pdata, listw = W.FOQ, test="rlme",model="within")

# Likelihood Maximum estimation
summary(inj_SEM_pool <- spml(fm_In, data=pdata, listw = W.FOQ, lag=TRUE,model="pooling"))
# Fixed-effect SEM
summary(inj_SEM_FE_In<- spml(fm_In, data=pdata, listw = W.FOQ, lag=TRUE,model="within", effect="individual", spatial.error="b"))
summary(inj_SEM_FE<- spml(fm_In, data=pdata, listw = W.FOQ, lag=TRUE,model="within", effect="individual", spatial.error="kkp"))
# Generalised moments method estimation
summary(inj_SEM_FE_GM <- spgm(fm_In, data=pdata, listw = W.FOQ, model="within", moments="fullweights",spatial.error = TRUE))

#Fatalities

fm_f<-LogFatalities ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin
##Tests
#Hausman test
print(hausman_panel_f<-phtest(fm_f, data = pdata)) # result suggest fixed effects

# Hausman test robust to spatial autocorrelation (splm)
print(spat_hausman_ML_SEM_f<-sphtest(fm_f,data=pdata,listw =W.FOQ, spatial.model = "error", method="ML"))
#Hausman test for spatial models
print(spat_hausman_ML_SAR_f<-sphtest(fm_f,data=pdata,listw =W.FOQ,spatial.model = "lag", method="ML"))

#Hausman test robust to spatial autocorrelation of errors leads to rejection of the null hypothesis on absence of correlation between individual effects and explanatory variables. For the rest of the empirical analysis, a fixed effects model is thus chosen.

# Fixed effects model
# Test 1
slmtest(fm_f, data=pdata, listw = W.FOQ, test="lml",model="within")

slmtest(fm_f, data=pdata, listw = W.FOQ, test="lme",model="within")

#confirm the rejection of the hypothesis that these two terms (taken independently) are null
#Robust tests
slmtest(fm_f, data=pdata, listw = W.FOQ, test="rlml", model="within")
slmtest(fm_f, data=pdata, listw = W.FOQ, test="rlme",model="within")

# Likelihood Maximum estimation
summary(f_SEM_pool <- spml(fm_f, data=pdata, listw = W.FOQ, lag=TRUE,model="pooling"))
# Fixed-effect SEM
summary(f_SEM_FE_b<- spml(fm_f, data=pdata, listw = W.FOQ, lag=TRUE,model="within", effect="individual", spatial.error="b"))
summary(f_SEM_FE_kkp<- spml(fm_f, data=pdata, listw = W.FOQ, lag=TRUE,model="within", effect="individual", spatial.error="kkp"))
# Generalised moments method estimation
summary(f_SEM_FE_GM <- spgm(fm_f, data=pdata, listw = W.FOQ, model="within", moments="fullweights",spatial.error = TRUE))

#DWI RESULTS
fm_D<-LogDWI ~ EmpPrcnt+MedHHInc+VehOwn0+PopDensity+TotalTrips+LogTripsRideAustin
##Tests
#Hausman test
print(hausman_panel_D<-phtest(fm_D, data = pdata)) # result suggest fixed effects

# Hausman test robust to spatial autocorrelation (splm)
print(spat_hausman_ML_SEM_D<-sphtest(fm_D,data=pdata,listw =W.FOQ, spatial.model = "error", method="ML"))
#Hausman test for spatial models
print(spat_hausman_ML_SAR_D<-sphtest(fm_D,data=pdata,listw =W.FOQ,spatial.model = "lag", method="ML"))

#Hausman test robust to spatial autocorrelation of errors leads to rejection of the null hypothesis on absence of correlation between individual effects and explanatory variables. For the rest of the empirical analysis, a fixed effects model is thus chosen.

# Fixed effects model
# Test 1
slmtest(fm_D, data=pdata, listw = W.FOQ, test="lml",model="within")
#coeff_plot
#
slmtest(fm_D, data=pdata, listw = W.FOQ, test="lme",model="within")

#confirm the rejection of the hypothesis that these two terms (taken independently) are null
#Robust tests
slmtest(fm_D, data=pdata, listw = W.FOQ, test="rlml", model="within")
slmtest(fm_D, data=pdata, listw = W.FOQ, test="rlme",model="within")

# Likelihood Maximum estimation
summary(D_SEM_pool <- spml(fm_D, data=pdata, listw = W.FOQ, lag=TRUE,model="pooling"))
# Fixed-effect SEM
summary(D_SEM_FE_b<- spml(fm_D, data=pdata, listw = W.FOQ, lag=TRUE,model="within", effect="individual", spatial.error="b"))
summary(D_SEM_FE_kkp<- spml(fm_D, data=pdata, listw = W.FOQ, lag=TRUE,model="within", effect="individual", spatial.error="kkp"))
# Generalised moments method estimation
summary(D_SEM_FE_GM <- spgm(fm_D, data=pdata, listw = W.FOQ, model="within", moments="fullweights",spatial.error = TRUE))



#Alternative hypothesis no random effects
test1 <- bsktest(x = fm, data = pdata, listw = W.FOQ, test = "LM1")
print(class(test1))
test1
#check for spatial autocorrelation on the errors of the model
test2 <- bsktest(x = fm, data = pdata, listw = W.FOQ, test = "LM2")
test2



#Fixed effects and spatial error B
fespaterr <- spml(fm, data = pdata, index = c("Geoid", "Time"), listw = W.FOQ, model="within", spatial.error="b")
summary(fespaterr)
#Fixed effects and spatial lagged error
fespatlag<- spml(fm_f, data = pdata, index = c("Geoid", "Time"), listw = W.FOQ, model = "within", lag = TRUE, spatial.error = "b")
summary(fespatlag)

#no panel considered
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

#Correlation Plots Outcome
library(corrplot)
df1<-mydata %>%
select(Crashes, Fatalities, Injuries, DWI, MedHHInc, TotalTrips, EmpPrcnt, MedHHInc, VehOwn0,PopDensity,TripsRideAustin)

corrplot.mixed(cor(df1), order="hclust", tl.col="black")
#Nicest plot to add in appendix
ggsave("Corellation_F.tiff", dpi=300)

tiff(file = "c_f_v2.tiff", res = 300)
corrplot.mixed(cor(df1), lower = "ellipse", upper = "number", tl.pos = "lt", diag = "u", tl.col="black")
dev.off()
ggcorr(df1, nbreaks=8, palette='RdGy', label=TRUE, label_size=4, label_color='black')
library("PerformanceAnalytics")

chart.Correlation(df1, histogram=TRUE, pch=19)