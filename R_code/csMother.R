##################################################################
## Clear Sky User Interface
##################################################################

## Clear Sky Models implemented
csModels <- c(
  "A2005",
  "AB1978",
  "ABCGS",
  "ASHRAE1972",
  "B1",
  "B1977",
  "B1979",
  "B2",
  "BCLS1979",
  "BD",
  "BD2014",
  "BDG2007",
  "BL1994",
  "BR1979",
  "C1978",
  "C1985",
  "DF2014",
  "DJW2012",
  "DMLA1988",
  "DPP",
  "DPPLT",
  "ESRA",
  "FR1999",
  "H1978",
  "HE1",
  "HE2",
  "HLJ",
  "HS",
  "I1981",
  "IA1983",
  "IB1983",
  "IC1983",
  "IP2002",
  "Josefsson",
  "JSP2011",
  "K1979",
  "K1980",
  "KC1980",
  "KSK1997",
  "L1970",
  "M1972",
  "M1976",
  "MAC",
  "METSTAT",
  "MI1985",
  "mKasten",
  "mP1982",
  "mP2000",
  "MRM4",
  "MRM5",
  "NAD1997",
  "PP1976",
  "PS2003",
  "PSA2000",
  "PSI",
  "PV1982",
  "REST2",
  "RS2000",
  "RSP1978",
  "S1976",
  "S1994",
  "SH1977",
  "SP1965",
  "sSOLIS",
  "W1978",
  "WL1967",
  "WMO1",
  "WMO2",
  "YH2001",
  "Z2006"
)

## Merge Sun geometry and Meteo objects
doMeteo <- function(sun, meteo){
        clMeteo <- class(meteo)
        switch(clMeteo,
               zoo = {
                   merge(sun, meteo)   
               },
               data.frame = {
                   if (nrow(meteo) != nrow(sun)){
                       if (nrow(meteo) == 1) {
                           nmsMeteo <- names(meteo)
                           meteo <- matrix(unlist(meteo),
                                           nrow = nrow(sun),
                                           ncol = ncol(meteo),
                                           byrow = TRUE)
                       } else {
                           stop('`meteo` and `sun` must have the same number of rows.')
                       }
                   }
                   sm <- cbind(sun, meteo)
                   names(sm) <- c(names(sun), nmsMeteo)
                   sm
               },
               matrix = {
                   if (nrow(meteo) != nrow(sun)){
                       if (nrow(meteo) == 1) {
                           nmsMeteo <- colnames(meteo)
                           meteo <- matrix(meteo,
                                           nrow = nrow(sun),
                                           ncol = ncol(meteo),
                                           byrow = TRUE)
                       } else {
                           stop('`meteo` and `sun` must have the same number of rows.')
                       }
                   }
                   sm <- cbind(sun, meteo)
                   if (is.null(nmsMeteo)) nmsMeteo <- seq_len(ncol(meteo))
                   names(sm) <- c(names(sun), nmsMeteo)
                   sm
               },
               numeric = {
                   nmsMeteo <- names(meteo)
                   meteo <- matrix(meteo,
                                   nrow = nrow(sun),
                                   ncol = length(meteo),
                                   byrow = TRUE)
                   sm <- cbind(sun, meteo)
                   if (is.null(nmsMeteo)) nmsMeteo <- seq_len(ncol(meteo))
                   names(sm) <- c(names(sun), nmsMeteo)
                   sm
               }
               )
    }

clearSky <- function(meteo, loc, model){
    stopifnot(model %in% csModels)
    if (missing(loc)) stop('`loc` is missing. You must provide a `list` with at least `lon` and `lat` components.')
    ## Sun geometry
    sun <- calcSol(lat = loc$lat, BTi = index(meteo))
    sun <- as.zooI(sun)
    ## Which sun geometry variables are already present in `meteo`?
    sunInMeteo <- names(sun) %in% names(meteo)
    ## Drop these variables from `sun` (`meteo` wins!) and merge with `meteo`
    data <- doMeteo(sun[, !sunInMeteo], meteo)
    ## Eval model with data
    z <- do.call(model,
                 list(data = data,
                      loc = loc)
                 )
    z
}

