# Load tdr and solaR packages
library(tdr) ## The following code needs version >= 0.13
library(solaR) ## Required by the function clearSky

## Set the working directory, using the folder where the repository has been cloned.
## Example: 
## setwd('~/GitHub/pcsol/')

## Load time series of Cabauw and Carpentras stations
load('data_example/cabauw.RData')
load('data_example/carpentras.RData')
## Information about the BSRN stations
BSRNcc <- read.csv('data_example/stations.csv')

## Load code
source('R_code/clearSky.R')
source('R_code/csMother.R')

##################################################################
## Auxiliary functions
##################################################################
## Evaluate a collection of models using meteorological data from a station
evalModels <- function(meteo, station,
                       G0models, D0models, Bnmodels)
{
  ## Estimate G0 with the collection of models
    estimG0 <- lapply(G0models, FUN = function(model)
    {
        estim <- clearSky(meteo, station, model = model)
        estim$G0
    })
    estimG0 <- do.call(cbind, estimG0)
    names(estimG0) <- G0models
    
    ## Estimate D0 with the collection of models
    estimD0 <- lapply(D0models, FUN = function(model)
    {
        estim <- clearSky(meteo, station, model = model)
        estim$D0
    })
    estimD0 <- do.call(cbind, estimD0)
    names(estimD0) <- D0models
    
    ## Bn
    estimBn <- lapply(Bnmodels, FUN = function(model)
    {
        estim <- clearSky(meteo, station, model = model)
        estim$Bn
    })
    estimBn <- do.call(cbind, estimBn)
    names(estimBn) <- Bnmodels

    ## Output
    list(G0 = estimG0,
         D0 = estimD0,
         Bn = estimBn)    
}

## Compute model errors using tdr::applyStats
radStats <- function(model, vals)
{
  errG0 <- applyStats(model$G0, vals$G0)
  errD0 <- applyStats(model$D0, vals$D0)
  errBn <- applyStats(model$Bn, vals$Bn)
  list(G0 = errG0,
       D0 = errD0,
       Bn = errBn)
}

##################################################################
## Load models
##################################################################
## The whole list of models is:
## print(csModels)

## Define which models estimate each component

myModelsG0 <- c("A2005", "ABCGS", "ASHRAE1972", "B1",  "B1979", "B2", "BCLS1979", "BD", "BD2014",
                "BDG2007","BL1994", "BR1979","C1978","C1985","DF2014","DJW2012","DMLA1988","DPP","DPPLT","ESRA","FR1999", "H1978",
                "HE1","HE2","HLJ","HS","I1981","IA1983", "IB1983", "IC1983", "IP2002","Josefsson","JSP2011","K1979","K1980","KC1980",
                "KSK1997","MAC","METSTAT","MI1985", "mKasten","mP1982", "mP2000", "MRM4", "MRM5", "NAD1997",
                "PP1976","PS2003","PSA2000","PSI","PV1982", "REST2","RS2000", "RSP1978","S1976","S1994","SH1977", "SP1965","sSOLIS",
                "W1978","WL1967","WMO1","WMO2","YH2001","Z2006")

myModelsD0 <- c("A2005", "AB1978", "ABCGS", "ASHRAE1972", "B1", "B1977", "B2", "BCLS1979",
                "BDG2007","BL1994", "BR1979","C1978","C1985","DF2014","DMLA1988","DPP","DPPLT","ESRA","FR1999", "H1978",
                "HE1","HE2","HLJ","HS", "I1981","IA1983","IB1983","IC1983","IP2002","Josefsson","JSP2011","K1979","K1980","KC1980",
                "KSK1997","MAC","METSTAT","MI1985", "mKasten","mP1982", "mP2000", "MRM4", "MRM5", "NAD1997",
                "PP1976","PS2003","PSA2000","PSI", "REST2", "RSP1978","S1976","S1994","SH1977", "SP1965","sSOLIS",
                "W1978","WL1967","YH2001","Z2006")


myModelsBn <- c("A2005", "AB1978", "ABCGS", "ASHRAE1972", "B1",  "B2", "BCLS1979",
                "BDG2007","BL1994", "BR1979","C1978","C1985","DF2014","DMLA1988","DPP","DPPLT","ESRA","FR1999", "H1978",
                "HE1","HE2","HLJ","HS", "I1981","IA1983","IB1983","IC1983","IP2002","Josefsson","JSP2011","K1979","K1980",
                "KSK1997","L1970","M1972","M1976","MAC","METSTAT","MI1985", "mKasten","mP1982", "mP2000", "MRM4", "MRM5", "NAD1997",
                "PP1976","PS2003","PSA2000","PSI", "REST2","RSP1978","S1976","S1994","SH1977", "SP1965","sSOLIS",
                "W1978","WL1967","YH2001","Z2006")

##################################################################
## Evaluate models
##################################################################
## Cabauw
estimCabauw <- evalModels(cabauw, BSRNcc[1,],
                          myModelsG0, myModelsD0, myModelsBn)
statsCabauw <- radStats(estimCabauw, cabauw)

## Carpentras
estimCarpentras <- evalModels(carpentras, BSRNcc[2,],
                          myModelsG0, myModelsD0, myModelsBn)

statsCarpentras <- radStats(estimCarpentras, carpentras)

##################################################################
## Target Diagrams
##################################################################
## Example for Cabauw station

## Define the graphical parameters
myTheme <- tdTheme(pch = c(21:25),
                   col.points = 'black',
                   fill = brewer.pal(n= 12, 'Paired'),
                   cex = 0.75)
## Define the legend
myLegend <- list(columns = 2,
                            cex = 0.8,
                            type = 'l',
                            space = 'right')

## G0: Display only those stations whose nRMSE is below 0.1
idxG0 <- which(statsCabauw$G0$nrmse < 0.1)

targetDiagram(statsCabauw$G0[idxG0, ],
              groups = model,
              type='at',
              cuts = seq(0, .1, .025),
              par.settings = myTheme,
              auto.key = myLegend)

## D0: Display only those stations whose nRMSE is below 0.1
idxD0 <- which(statsCabauw$D0$nrmse < 0.1)

targetDiagram(statsCabauw$D0[idxD0, ],
              groups = model,
              type = 'at',
              cuts = seq(0, .1, .025),
              par.settings = myTheme,
              auto.key = myLegend)

## Bn: Display only those stations whose nRMSE is below 0.105
idxBn <- which(statsCabauw$Bn$nrmse < 0.105)

targetDiagram(statsCabauw$Bn[idxBn,],
              groups = model,
              type = 'at',
              cuts = seq(0, .1, .025),
              par.settings = myTheme,
              auto.key = myLegend)


