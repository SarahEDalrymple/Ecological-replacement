# DESCRIPTION
# Species Distribution Modeling with biomod2.
# analysis aims to use occurence records to generate projections of climate niche,
# uses current and future climate projections to calculate species range change

# Sarah Dalrymple
# for sources, see end of script

# Script written in R version 3.5.1 (2018-07-02) -- "Feather Spray"
# Copyright (C) 2018 The R Foundation for Statistical Computing


# Notes on use of script:
# Script should be adapted for different species, geographic areas and climate projections but,
# save the script with an appropriate file name and do not change this master version.


### Step 0: preparing R for your analysis
#########################################

# for first time use of code install all necessary packages


# load libraries needed for analysis
# if not yet installed, R will return an error message in the console,
# use install.packages() inserting the name of the package in quotation marks in the brackets,
# and then try loading the libraries again

library(ade4)
library(sp)#
library(maps)
library(maptools)
library(rgdal)#
library(raster)#
library(rasterVis)
library(dismo)
library(ecodist)
library(biomod2)#
library(rgeos)#
library(usdm)
library(abind)
library(gridExtra)
library(lattice)
library(reshape2)#
library(ggplot2)#
library(dplyr)
library(markdown)
library(ggmap)
library(rgbif)

# insert file pathway with setwd() i.e. something like this: setwd("C:\\Users\\Joe\\Documents\\R")
# *IMPORTANT* if working with various species, it's a good idea to set the 
# working directory to a specific folder for each species

file.choose()
setwd("C:\\Users\\sarah\\Dropbox\\R materials\\Biomod\\biomod2_workshop")

### step 1: read distribution data
##################################


# Read species occurrences into R from your own text file
# to use code below unchanged, your text file should have x,
# or longitudinal coordinates in the first column,
# and y or latitudinal coordinates in the second column

SpOcc <-read.table("Melsyl_coordinates.txt", 
                   header = TRUE)

head(SpOcc) # displays the first five rows of object 'SpOcc'

# DATA PREPARATION
# this tells the analysis what projection system we're using for spatial data
ProjW = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs" 
#map projection (EPSG: 4326)

#occurences in SpatialPointsDataFrame() format required by function 'over' and biomod2
# square brackets below specify which rows (before comma) and columns (after comma) refer to 'xy' object
# if x and y coordinates are in different rows, change these numbers accordingly
xy <- SpOcc[,1:2]
df <- SpOcc
Occ <- SpatialPointsDataFrame(coords=xy, data=df, proj4string = CRS(ProjW))

plot(Occ)

# plot() above, will plot your datapoints spatially to check they look okay,
# but this may not mean much without a map to set them against. If you want to use a map,
# download world map files ('countries.shp') to the working directory
# from http://www.diva-gis.org/Data and use readOGR() to tell R which shapefile to read

global_map <-readOGR("countries.shp")

# use the two lines below to check that the shapefiles are on the same projection
proj4string(global_map)
proj4string(Occ)

# if they're not, use the spTransform() to change the projections to match,
# in this case, change the base map because 'Occ' is already in the required format for biomod
global_map <- spTransform(global_map, CRSobj = CRS(proj4string(Occ)))
plot(global_map)
points(Occ, pch = 16, col = "red")

# given that there are some odd occurrences, and we want the run time to be short,
# we'll actually work with a reduced dataset from this point on

SpOcc <-read.table("Melsyl_coordinates_GB_IE.txt", 
                   header = TRUE)

head(SpOcc) # displays the first five rows of object 'SpOcc'

# DATA PREPARATION
# this tells the analysis what projection system we're using for spatial data
ProjW = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs" 
#map projection (EPSG: 4326)

#occurences in SpatialPointsDataFrame() format required by function 'over' and biomod2
# square brackets below specify which rows (before comma) and columns (after comma) refer to 'xy' object
# if x and y coordinates are in different rows, change these numbers accordingly
xy <- SpOcc[,1:2]
df <- SpOcc
Occ <- SpatialPointsDataFrame(coords=xy, data=df, proj4string = CRS(ProjW))

plot(Occ)
ext <- extent(-11, 2, 50, 61)
GB_IE <- crop(global_map, ext)
plot(GB_IE)
points(Occ, pch = 16, col = "red")


#### Step 2:  create raster stack of current climate layers from WorldClim
###########################################################################


# step 2 should generate an object called 'bioclim_world' 
# Use the `getData()` function from the 'raster' package to access climate data

# to retreive global current bioclimatic variables at 2.5' resolution
bioclim_world <- getData("worldclim", res = 2.5, var = "bio")

bioclim_world #gives summary of object, in this case 'bioclim_world', including extent and class
              # A raster stack is a collection of many raster layers with the same projection,
              # spatial extent and resolution.
plot(bioclim_world) 
# plot lets you check that they have loaded properly but is not necessary to the analysis


### Step 3: create stack for restricted area
#############################################

# this section of code restricts the extent of your climate data
# it allows maps to be resized appropriately and processing time is reduced
# you should generate an object (in this case a raster stack) called 'bioclim_cropped'

# limit area by longitudinal and latitudinal extent 
# e.g. for UK and Ireland
# limit the map by lon/lat, figures given below are for UK and Ireland (-11, 2, 50, 61)
bioclim_cropped <- crop(bioclim_world, ext)
plot(bioclim_cropped) #plots the bioclim variables only for UK and Ireland


### Step 4: selecting subset of variables to avoid correlated variables
########################################################################

# at the end of this step, you should have a raster stack of bioclim variables that 
# avoid colinearity and is cropped to the area you need.
# the raster stack is called 'clim'

#convert stack into dataframe
current_df <- as.data.frame(bioclim_cropped)
current_df <- na.omit(current_df) #removes NAs

##calculate Pearson correlations between pairs of variables
cor_current <- cor(current_df)

##reformat correlation table for graphical analyses
cor_current [upper.tri(cor_current, diag = TRUE)] <- NA
cor_current_resh <- na.omit(melt(cor_current))
colnames(cor_current_resh) <- c("var1", "var2", "correlation")

##only consider the absolute value of correlations
cor_current_resh$correlation <- abs(cor_current_resh$correlation)

##correlation plot
gg_cor <- ggplot(cor_current_resh, aes(x = var1, y = var2, fill = correlation))
gg_cor <- gg_cor +geom_tile() + xlab("") + ylab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

#darker colours = lower the correlation
print(gg_cor)

## select bio17, bio3, bio9
selected_vars <-c("bio17", "bio3", "bio9")

##check correlation between selected variables
(cor_sel <- cor(current_df[,selected_vars]))

#Stack variables - only include bioclim variables that were selected in colinearity checks
clim <-stack(bioclim_cropped$bio3, bioclim_cropped$bio17, bioclim_cropped$bio9)
plot(clim)


### step 5: format data for biomod
###################################

myBiomodData <- BIOMOD_FormatingData(resp.var = rep(1, nrow( Occ )),
                                   expl.var = clim,
                                   resp.xy = xy,
                                   resp.name = "Melampyrum.sylvaticum",
                                   PA.nb.rep = 3,
                                   PA.nb.absences = 500,
                                   PA.strategy = 'random',
                                   na.rm = TRUE)


## plot of selected pseudo-absences
plot(myBiomodData)


myBiomodModelOut <- BIOMOD_Modeling(data = myBiomodData,
                                         models = c('GAM','GBM','RF'),
                                         models.options = BIOMOD_ModelingOptions(),
                                         NbRunEval=3,
                                         DataSplit=70,
                                         Yweights=NULL,
                                         VarImport=3, 
                                         models.eval.meth = c('KAPPA', 'TSS', 'ROC'),
                                         SaveObj = TRUE,
                                         rescal.all.models = FALSE, 
                                         do.full.models = FALSE)
## get models evaluation scores
myBiomodModelOut_scores <- get_evaluations(myBiomodModelOut)
## myBiomodModelOut_scores is a 5 dimension array containing the scores of the models
dim(myBiomodModelOut_scores)
dimnames(myBiomodModelOut_scores)
models_scores_graph(myBiomodModelOut, by = "models" , metrics = c("ROC","TSS"), 
                    xlim = c(0.5,1), ylim = c(0.5,1))
models_scores_graph(myBiomodModelOut, by = "cv_run" , metrics = c("ROC","TSS"), 
                    xlim = c(0.5,1), ylim = c(0.5,1))
models_scores_graph(myBiomodModelOut, by = "data_set" , metrics = c("ROC","TSS"), 
                    xlim = c(0.5,1), ylim = c(0.5,1))
(myBiomodModelOut_var_import <- get_variables_importance(myBiomodModelOut))
## make the mean of variable importance by algorithm
apply(myBiomodModelOut_var_import, c(1,2), mean)

meanVarImport_gbm <- BIOMOD_LoadModels(myBiomodModelOut, models='GBM')
meanVarImport_rf <- BIOMOD_LoadModels(myBiomodModelOut, models='RF')
meanVarImport_gam <- BIOMOD_LoadModels(myBiomodModelOut, models='GAM')

gbm_eval_strip <- biomod2::response.plot2(
  models  = meanVarImport_gbm,
  Data = get_formal_data(myBiomodModelOut,'expl.var'), 
  show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric = 'median',
  legend = FALSE,
  display_title = FALSE,
  data_species = get_formal_data(myBiomodModelOut,'resp.var'))
rf_eval_strip <- biomod2::response.plot2(
  models  = meanVarImport_rf,
  Data = get_formal_data(myBiomodModelOut,'expl.var'), 
  show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric = 'median',
  legend = FALSE,
  display_title = FALSE,
  data_species = get_formal_data(myBiomodModelOut,'resp.var'))
gam_eval_strip <- biomod2::response.plot2(
  models  = meanVarImport_gam,
  Data = get_formal_data(myBiomodModelOut,'expl.var'), 
  show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric = 'median',
  legend = FALSE,
  display_title = FALSE,
  data_species = get_formal_data(myBiomodModelOut,'resp.var'))


#### step 6: BIOMOD Ensemble models
###################################

# chosen.models = models kept for building ensemble models
# em.by = Defines the way the models will be combined to build the ensemble models. See vignette or this link for more
#     info on this one...  https://rstudio-pubs-static.s3.amazonaws.com/38564_747d4bbf87704f0394734977bd4905c4.html
# eval.metric = evaluation metric(s) used to build ensemble models
# eval.metric.quality.threshold = If not NULL, then the minimum values required for models to be included in the ensemble forecast
# prob.mean = Logical. Estimate the mean probabilities across predictions
# prob.cv = Logical. Estimate the coefficient of variation across predictions. The lower is the score, the better are the 
#       models. CV is a nice complement to the mean probability.
# prob.ci = Logical . Estimate the confidence interval around the prob.mean
# prob.ci.alpha = Numeric. Significance level for estimating the confidence interval. Default = 0.05
# prob.median = Logical. Estimate the mediane of probabilities
# committee.averaging = Logical. Estimate the committee averaging across predictions
# prob.mean.weight = Logical. Estimate the weighted sum of probabilities
# prob.mean.weight.decay = Defines the relative importance of the weights. A high value will strongly discriminate the 'good' 
#     models from the 'bad' ones (see the details section). If the value of this parameter is set to 'proportional' (default), 
#     then the attributed weights are proportional to the evaluation scores given by 'weight.method'(eval.metric)


ensemble_models <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut,
                                            chosen.models = 'all',       
                                            em.by = 'all',
                                            eval.metric = c('KAPPA', 'TSS', 'ROC'),
                                            eval.metric.quality.threshold = c(0.3, 0.5, 0.7),
                                            models.eval.meth = c('KAPPA','TSS','ROC'),
                                            prob.mean = FALSE,
                                            prob.cv = TRUE, 
                                            committee.averaging = TRUE,
                                            prob.mean.weight = TRUE,
                                            VarImport = 0 )
(ensemble_models_scores <- get_evaluations(ensemble_models))


#### step 7: Current projections
################################

#individual projections (to ecoregion-analog climates in occupied biomes)
# modeling.output = "BIOMOD.models.out" object produced by a BIOMOD_Modeling run
# new.env = A set of explanatory variables onto which models will be projected. 
# proj.name = a new folder will be created with this name
# Make sure the column names (data.frame or matrix) or layer Names (rasterStack),
# perfectly match with the names of variables used to build the models in the previous steps.

models_proj_current <- BIOMOD_Projection( modeling.output = myBiomodModelOut,
                                                 new.env = clim,
                                                 proj.name = "current", #projection to analogous climates in occupied range
                                                 binary.meth = "TSS",
                                                 output.format = ".img",
                                                 do.stack = FALSE )

# myCurrentProj <-  get_predictions(models_proj_current)#predictions for each run (maps)

ensemble_models_proj_current <- 
  BIOMOD_EnsembleForecasting( EM.output = ensemble_models,
                              projection.output = models_proj_current,
                              binary.meth = "TSS",
                              output.format = ".img",
                              do.stack = FALSE )


### step 8: get future climate projections
##########################################

# load 2050 and 2070 bioclim variables
# correspond to step 2 but generates a raster stack of bioclim variables based on future projections
# resolution, models, RCPs (IPCC sceanrios) and time slices are optional and can be altered
# if these are altered, the object names should be altered accordingly
# all steps should create a stack named 'clim_world_[add year and model/RCP]'

# to retrieve future climate projections, CMIP5, at 2.5' resolution for 2050 and 2070,
# additional arguments specify RCP 8.5, CCSM4 model
?getData
bioclim_world_50 <- getData('CMIP5', res = 2.5, var = "bio", rcp=85, model="CC", year=50)
bioclim_world_70 <- getData('CMIP5', res = 2.5, var = "bio", rcp=85, model="CC", year=70)

# select same subset of bioclim variables as used in current projections

# 2050
clim_world_50_CC85 <-stack(c(bio3 = bioclim_world_50$cc85bi503, bio17 = bioclim_world_50$cc85bi5017, bio9 = bioclim_world_50$cc85bi509))
plot(clim_world_50_CC85)

# 2070
clim_world_70_CC85 <-stack(c(bio3 = bioclim_world_70$cc85bi703, bio17 = bioclim_world_70$cc85bi7017, bio9 = bioclim_world_70$cc85bi709))
plot(clim_world_70_CC85)


#### step 9: crop future bioclim variables
###########################################

# the same method of cropping should be used as in step 3 to ensure that projections are working to the same area
# step 9 should result in an object or objects which is a raster stack called 'clim_[year]_[model/RCP]'

# limit area by longitudinal and latitudinal extent 
# 'ext' was created earlier in step 1,
# it should encompass the spatial extent of the occurence records
clim_50_CC85 <- crop(clim_world_50_CC85, ext)
clim_50_CC85 <- stack(clim_50_CC85)
plot(clim_50_CC85)

clim_70_CC85 <- crop(clim_world_70_CC85, ext)
clim_70_CC85 <- stack( clim_70_CC85 )
plot(clim_70_CC85)



#### step 10: future projections
################################

models_proj_2050_CC85 <- BIOMOD_Projection( modeling.output = myBiomodModelOut,
                                                   new.env = clim_50_CC85,
                                                   proj.name = "2050_CC85",
                                                   binary.meth = "TSS",
                                                   output.format = ".img",
                                                   do.stack = FALSE )
ensemble_models_proj_2050_CC85 <- 
  BIOMOD_EnsembleForecasting( EM.output = ensemble_models,
                              projection.output = models_proj_2050_CC85,
                              binary.meth = "TSS",
                              output.format = ".img",
                              do.stack = FALSE )
plot(ensemble_models_proj_2050_CC85, 
     str.grep = "EMwmean")

models_proj_2070_CC85 <- BIOMOD_Projection( modeling.output = myBiomodModelOut,
                                                   new.env = clim_70_CC85,
                                                   proj.name = "2070_CC85",
                                                   binary.meth = "TSS",
                                                   output.format = ".img",
                                                   do.stack = FALSE )
ensemble_models_proj_2070_CC85 <- 
  BIOMOD_EnsembleForecasting( EM.output = ensemble_models,
                              projection.output = models_proj_2070_CC85,
                              binary.meth = "TSS",
                              output.format = ".img",
                              do.stack = FALSE )
plot(ensemble_models_proj_2070_CC85, 
     str.grep = "EMwmean")


#### step 11: calculate projected species range change
######################################################

## load binary projections

bin_proj_current <- stack( 
  c( wm = "C:/Users/sarah/Dropbox/Melampyrum/Melampyrum_biomod/Melampyrum.sylvaticum/proj_current/individual_projections/Melampyrum.sylvaticum_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img") )
bin_proj_2050_CC85 <- stack( 
  c( wm = "C:/Users/sarah/Dropbox/Melampyrum/Melampyrum_biomod/Melampyrum.sylvaticum/proj_2050_CC85/individual_projections/Melampyrum.sylvaticum_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img") )
bin_proj_2070_CC85 <- stack( 
  c( wm = "C:/Users/sarah/Dropbox/Melampyrum/Melampyrum_biomod/Melampyrum.sylvaticum/proj_2070_CC85/individual_projections/Melampyrum.sylvaticum_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img") )

## SRC current -> 2050
SRC_current_2050_CC85 <- BIOMOD_RangeSize( bin_proj_current,
                                           bin_proj_2050_CC85 )
SRC_current_2050_CC85$Compt.By.Models

## SRC current -> 2070
SRC_current_2070_CC85 <- BIOMOD_RangeSize( bin_proj_current,
                                           bin_proj_2070_CC85 )
SRC_current_2070_CC85$Compt.By.Models

src_map <- stack(SRC_current_2050_CC85$Diff.By.Pixel, SRC_current_2070_CC85$Diff.By.Pixel)
names(src_map) <- c("wm cur-2050", "wm cur-2070")

my.at <- seq(-2.5,1.5,1)
myColorkey <- list(at=my.at, ## where the colors change
                   labels=list(
                     labels=c("lost", "pres", "abs","gain"), ## labels
                     at=my.at[-1]-0.5 ## where to print labels
                   ))
rasterVis::levelplot( src_map, 
                      main = "Melampyrum sylvaticum range change",
                      colorkey = myColorkey,
                      layout = c(2,1) )
#Sources:

#R
citation()

*Citation:* 
  @book{
    title={Habitat Suitability and Distribution Models: With Applications in R},
    author={Guisan, A. and Thuiller, W. and Zimmermann, N.E.},
    isbn={9780521758369},
    series={Ecology, Biodiversity and Conservation},
    year={2017},
    publisher={Cambridge University Press}
  