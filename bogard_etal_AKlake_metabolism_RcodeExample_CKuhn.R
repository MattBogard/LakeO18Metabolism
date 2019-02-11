##################################################################################
## O18 Metabolic estimates in lakes - mass balance steps
##################################################################################
#importing field and lab data for O18 mass balance metabolism modelling#######
met = read.csv("aklakes_met_input_terms.csv", 
              header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
summary(met)

#notes:
#Wind speeds measured as hourly averages on Canvasback lake (USGS weather station)
#Wind speed data for trip #3 : 3-day Average (Jun 28 to 30, 2016) of daily average (2.88 m/s) and mean daily STDEV (+/- 1.24 m/s) 
#Wind speed data for trip #4 : 3-day Average (Sept 09-12, 2016) of daily average (1.77 m/s) and mean daily STDEV (+/- 1.18 m/s) 
#Wind data in Sept do not overlap with lake sampling b/c weather station was removed from lake



#function for oxygen isotope mass balance calculations of metabolic rates and balances###########
met.fun <- function(args){
  
  #all input terms added and renamed here
  wind.ms <- met$wind.ms
  wind.height <- met$wind.height.ms
  area <- met$lake.area.km2
  temp <- met$h2o.temp.c
  zmix <- met$zmix
  do <- met$o2.g.m3
  do.pct.sat <- met$do.pct.sat
  delo18.o2 <- met$delo18.o2
  delo18.h2o <- met$delo18.h2o
  
  #calculating the gas exchange coefficient for O2 empirically from lake area and wind speed:
  u10 <- wind.ms * (1+ (((0.0013^(0.5))/0.41) * (log(10/wind.height)))) #converting wind speed from 3m to 10m height following equation 3 Vachon & Prairie (2013)
  k600cmh <- 2.51 + 1.48*u10 + 0.39*u10*(log10(area)) #k600 in cm/h from table 2 equation B vachon & prairie 2013
  k600md <- k600cmh * 24/100 #converting k600 to m/d
  sco2 <- 1800.6 - (120.1*temp) + (3.7818 * (temp^2)) - (0.047608*(temp^3))#calculating schmidt number for oxygen from Jahne et al (1987)
  ko2md <- k600md * ((sco2/600)^(-2/3)) #converting k600 to ko2 in m/d for use in mass balance
  #2/3 power used for wind speed less than 3.7 m/s following Vachon et al. (2010) and Guerin et al. (2007)
  k.z <- ko2md / zmix #input term for volumetric gas exchange
  
  #atom fractions calculated for input into mass balance:
  dosat <- do * 100/do.pct.sat #calculate dissolved oxygen concentration at equilibrium with atmosphere 
  ro18o2<-((delo18.o2/1000)*(0.0020052))+(0.0020052)#converting del value of DO-O18 back to O18:O16 ratio
  ro18h2o <- ((delo18.h2o/1000)*(0.0020052))+(0.0020052) #converting del value of H2O-O18 back to O18:O16 ratio
  ro18air = ((23.88/1000)*(0.0020052))+(0.0020052) #converting del value of O2-atmosphere back to O18:O16 ratio
  #switched to d18O-air value = 23.88 permil following Barkan & Luz 2005 
     o18o2 <- ro18o2 / (1 + ro18o2) #converting ratio of O18:O16 in DO to atom fraction following Bogard et al. (2017)
  o18h2o <- ro18h2o / (1 + ro18h2o) #converting ratio of O18:O16 in H2O to atom fraction following Bogard et al. (2017)
  o18air <- ro18air/(1 + ro18air) #fixed value: atom fraction of O18 as AF=R/(R+1)

    #generic fractionation factors
  ffp <- 1           # fractionation factor associated with photosynthesis
  ffg <- 0.9972      # fractionation factor associated with gas exchange
  ffs <- 1.0007      # fractionation factor associated with DO solubility
  ffr <- 0.985        # fractionation factor associated with ecosystem respiration - a weighted approximate average of all DO consumption processes
  
  # summarized equation terms as presented by Bocaniov et al. 2012
  a <- ffg * ffs * o18air
  b <- ffg * o18o2
  c <- ffr * o18o2
  d <- ffp * o18h2o
  
  
  gppv <- (k.z * (((do * (b - c)) - (dosat * (a - c))) / (d - c))) #volumetric GPP
  rv <- (k.z * (((do * (b - d)) - (dosat * (a - d))) / (d - c))) #volumetric R
  nepv <- gppv - rv #volumetric NEP
  gppa <- zmix * (k.z * (((do * (b - c)) - (dosat * (a - c))) / (d - c))) #areal GPP
  ra <- zmix * (k.z * (((do * (b - d)) - (dosat * (a - d))) / (d - c))) #areal R
  nepa <- gppa - ra #areal NEP
  # P:R (from Quay 1995 L&O) 
  g <- (ffg*((ffs * o18air) - ((do / dosat) * o18o2))) / (1 - (do / dosat)) #Term G directly from Quay 1995 equation
  #PtoR ratios calculated using a fractionation value of 0.985
  gpptor <- ((o18o2 * ffr) - g) / ((o18h2o * ffp) - g)
  
  output<- data.frame(gppv,gppa,rv,ra,nepv,nepa,gpptor,a,b,c,d,g)
  return(output)
}

#run the metabolism function on the input data frame#####
met.result = met.fun(met)
#print out of metabolism summary stats ####### 
summary(met.result)
#mergining metabolism results into the original input data frame #######
met.final <- cbind(met, met.result)
#exporting metabolism data as csv file#####
write.csv(met, file = "AK_MET_R_EXPORT.csv")

