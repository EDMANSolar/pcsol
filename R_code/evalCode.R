## Load of input and output datasets for Cabauw and Carpentras
setwd('~/GitHub/pcsol/data_example')
load('data.RData')

# Load tdr and solaR libraries
library(tdr)
library(solaR)

setwd('~/GitHub/pcsol/R_code')
source('clearSky.R')
source('csMother.R')

BSRNcc <- read.csv('~/GitHub/pcsol/data_example/stations.csv')

##################################################################
## Load of models
##################################################################

# define which models estimate each component

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

Rad_est <- lapply(c(1:2), function(st){
  Scc <- BSRNcc[st,]
  Smeteo <- L[[st]]
  
  # G0
  G0models <- lapply(myModelsG0, FUN = function(model){
    estim <- clearSky(Smeteo, Scc, model = model)
    G0estim <- estim$G0
  })
  G0models <- do.call(cbind, G0models)
  names(G0models) <- myModelsG0
  
  # D0
  D0models <- lapply(myModelsD0, FUN = function(model){
    estim <- clearSky(Smeteo, Scc, model = model)
    D0estim <- estim$D0
  })
  D0models <- do.call(cbind, D0models)
  names(D0models) <- myModelsD0
  
  # Bn
  Bnmodels <- lapply(myModelsBn, FUN = function(model){
    estim <- clearSky(Smeteo, Scc, model = model)
    Bnestim <- estim$Bn
  })
  Bnmodels <- do.call(cbind, Bnmodels)
  names(Bnmodels) <- myModelsBn
  
  Est <- list(G0=G0models, D0=D0models, Bn=Bnmodels)    
  Est
})

save(Rad_est, file='Rad_est.RData')


Rad_stats <- lapply(c(1:2), function(st){
  errModel_G0 <- applyStats(Rad_est[[st]]$G0, L[[st]]$G0)
  errModel_D0 <- applyStats(Rad_est[[st]]$D0, L[[st]]$D0)
  errModel_Bn <- applyStats(Rad_est[[st]]$Bn, L[[st]]$Bn)
  Stats <- list(G0=errModel_G0, D0=errModel_D0, Bn=errModel_Bn)
  Stats
})

save(Rad_stats, file='Rad_stats.RData')

## Example for station 1

foo <- which(Rad_stats[[1]]$G0$nrmse<0.1)


targetDiagram(Rad_stats[[1]]$G0[foo,], groups = model,type='at',
              cuts=seq(0, .1, .025),  par.settings=tdTheme(pch = c(21:25),
                                                           col.points = 'black',fill = brewer.pal(n= 12, 'Paired'),cex = 0.75), 
              auto.key=list(columns=2, cex=0.8, type='l', space='right'))



foo <- which(Rad_stats[[1]]$D0$nrmse<0.1)

targetDiagram(Rad_stats[[1]]$D0[foo,], groups = model,type='at',
              cuts=seq(0, .1, .025),  par.settings=tdTheme(pch = c(21:25),
                                                           col.points = 'black',fill = brewer.pal(n= 12, 'Paired'),cex = 0.75), 
              auto.key=list(columns=2, cex=0.8, type='l', space='right'))


foo <- which(Rad_stats[[1]]$Bn$nrmse<0.105)

targetDiagram(Rad_stats[[1]]$Bn[foo,], groups = model,type='at',
              cuts=seq(0, .1, .025),  par.settings=tdTheme(pch = c(21:25),
                                                           col.points = 'black',fill = brewer.pal(n= 12, 'Paired'),cex = 0.75), 
              auto.key=list(columns=2, cex=0.8, type='l', space='right'))


