##################################################################
## Clear Sky Models
##################################################################

## 1. ESRA 

# Rigollier, C., Bauer, O., Wald, L., 2000. On the clear sky model of the ESRA - European Solar Radiation Atlas-
# with respect to the heliosat method. Solar Energy 68(1), 33-48.

ESRA <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  y <- unique(year(index(data)))
  julianday <- as.numeric(julian(as.Date(index(data)), origin = as.Date(paste0(y,'-01-01')))) + 1
  
  TL <- data$TL
  
  elev <- loc$elev

  TLAM2 <- TL/0.8662
  # solar constant
  I0 <- 1367.0
  # day angle
  day_angle<- julianday*2.0*pi/365.2422
  # Correction accounting for eccentricity
  eo <- 1 + 0.03344*cos(day_angle - 2.80*pi/180.0)
  # relative air mass at sea level (Kasten's air mass), denoted m and ms1 at sea level
  # sun elev correction for refraction
  
  theta_c <- 0.061359*(0.1594 + 1.1230*(90 - theta) + 
                         0.065656*(90 - theta)^2)/(1 + 28.9344*(90 - theta) + 277.3971*(90 - theta)^2)
  theta_c <- theta + theta_c
    
  elev_c <- exp(-elev/8434.5)
  # air masses
  m <- elev_c/(cosThzS + 0.50572*((90 - theta) + 6.07995 )^-1.6364)
  m[(theta>=90)] <- 0
  ms1 <- 1/(cosThzS + 0.50572*((90 - theta) + 6.07995 )^-1.6364)
  # optical depth Rayleigh Page (1996)
  deltaR <- 1/(10.4 + 0.718*m)
  
  # pressure correction Page J., 2001 (modified rigollier at el. 2000 for pressure correction)  
  corr75 <- 1.248174 - 0.011997*ms1 + 0.00037*ms1^2
  corr50 <- 1.68219 - 0.03059*ms1 + 0.00089*ms1^2
  x <- (elev_c - 0.75)^2
  y <- (elev_c - 1.0)^2
  if(elev_c > 0.75){if(elev_c >= 0.99){press_corr <- rep(1, length(m))}else{press_corr <- (x + y*corr75)/(x + y)}}else{
    press_corr <- (x*corr50 + y*corr75)/(x + y)}
  ms1[is.na(ms1)] <- 0
  deltaR[(ms1<=20)] <- (1/(6.625928 + 1.92969*ms1[(ms1<=20)]  + -0.170073*(ms1[(ms1<=20)])^2 +
                        0.011517*(ms1[(ms1<=20)])^3 + -0.000285*(ms1[(ms1<=20)])^4))/press_corr[(ms1<=20)]
  
  a <- ((TLAM2 <= 1.0) | (m <= 0.0))
  Bn <- I0*eo*exp(-0.8662*TLAM2*m*deltaR)
  Bn[a] <- 0
  B0 <- Bn*cosThzS
  B0[(B0<0)] <- 0
  
  # correction of height
  TLAM2_c <- TLAM2*elev_c
  # Diffuse transmittance at zenith Trd
  T_d <- -1.5843e-2 + 3.0543e-2*TLAM2_c + 3.797e-4*TLAM2_c^2
  a <- matrix(c(2.6463e-1, -6.1581e-2, 3.1408e-3, 2.0402, 1.8945e-2,
                -1.1161e-2, -1.3025, 3.9231e-2, 8.5079e-3), nrow=3, ncol=3, byrow=T)
  A0 <- as.numeric(a[1,1] + a[1,2]*TLAM2_c + a[1,3]*TLAM2^2)
  A1 <- as.numeric(a[2,1] + a[2,2]*TLAM2_c + a[2,3]*TLAM2^2)
  A2 <- as.numeric(a[3,1] + a[3,2]*TLAM2_c + a[3,3]*TLAM2^2)
  
  A0[is.na(TLAM2_c)==TRUE] <- 2.0e-3/T_d

  Fd <- A0 + A1*cosThzS + A2*cosThzS^2
  D0 <- I0*eo*T_d*Fd 
  D0[D0<0] <- 0
  G0 <- D0 + B0 
  G0[G0<0] <- 0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}
  

## 2. simplified SOLIS. CHECKED with original

# Ineichen, P., 2008. A broadband simplified version of the Solis clear sky model. Solar Energy 82(8), 758–762.

sSOLIS <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  y <- unique(year(index(data)))
  julianday <- as.numeric(julian(as.Date(index(data)), origin = as.Date(paste0(y,'-01-01')))) + 1
  
  pwc <- data$pwc
  aod380 <- data$aod380
  aod500 <- data$aod500
  aod550 <- data$aod550
  aod700 <- data$aod700
  alpha <- data$alpha
  
  elev <- loc$elev
  
  # Bird, R.E., Huldstrom, R.L., 1980. Direct insolation models. 
  # Transactions ASME Journal of Solar Energy Engineering 103, 182-192.
  a <- ((is.na(aod700)==TRUE) & (is.na(aod380*aod500)==FALSE))
  aod700[a] <- 0.27583*aod380[a] + 0.35*aod500[a]
  b <- ((is.na(aod700)==TRUE) & (is.na(aod550*alpha)==FALSE))
  aod700[b] <- aod550[b]*(7/5.5)^-alpha[b]
  
  I0 <- 1367
  p2p0 <- exp(-elev/8434.5)
  day_angle <- julianday*2.0*pi/365.2422
  eo <- 1 + 0.03344*cos(day_angle - 2.80*pi/180)
  I00 <- 1.08*pwc^0.0051
  I01 <- 0.97*pwc^0.032
  I02 <- 0.12*pwc^0.56
  I0prima <- I0*(I02*aod700^2 + I01*aod700 + I00 + 0.071*log(p2p0))
  tb1 <- 1.82 + 0.056*log(pwc) + 0.0071*(log(pwc))^2
  tb0 <- 0.33 + 0.045*log(pwc) + 0.0096*(log(pwc))^2
  tbp <- 0.0089*pwc + 0.13
  taub <- tb1*aod700 + tb0 + tbp*log(p2p0)
  b1 <- 0.00925*aod700^2 + 0.0148*aod700 - 0.0172
  b0 <- -0.7565*aod700^2 + 0.5057*aod700 + 0.4557
  b <- b1*log(pwc) + b0
  Bn <- I0prima*exp(-taub/cosThzS^b)
  B0 <- Bn*cosThzS
  tg1 <- 1.24 + 0.047*log(pwc) + 0.0061*(log(pwc))^2
  tg0 <- 0.27 + 0.043*log(pwc) + 0.0090*(log(pwc))^2
  tgp <- 0.0079*pwc + 0.1
  taug <- tg1*aod700 + tg0 + tgp*log(p2p0)
  g <- -0.0147*log(pwc) - 0.3079*aod700^2 + 0.2846*aod700 + 0.3798
  G0 <- I0prima*exp(-taug/cosThzS^g)*cosThzS
  
  td4 <- -0.21*pwc + 11.6
  td3 <- 0.27*pwc - 20.7
  td2 <- -0.134*pwc + 15.5
  td1 <- 0.0554*pwc - 5.71
  td0 <- 0.0057*pwc + 2.94
  tdp <- -0.71*(1 + aod700)^(-15.0)
  
  a <- (aod700 < 0.05)
  td4[a] <- 86*pwc[a] - 13800
  td3[a] <- -3.11*pwc[a] + 79.4
  td2[a] <- -0.23*pwc[a] + 74.8
  td1[a] <- 0.092*pwc[a] - 8.86
  td0[a] <- 0.0042*pwc[a] + 3.12
  tdp[a] <- -0.83*(1 + aod700[a])^(-17.2)
  
  taud <- td4*aod700^4 + td3*aod700^3 + td2*aod700^2 + td1*aod700 + td0 + tdp*log(p2p0)
  dp <- 1/(18 + 152*aod700)
  d <- -0.337*aod700^2 + 0.63*aod700 + 0.116 + dp*log(p2p0)
  D0 <- I0prima*exp(-taud/cosThzS^d)
  
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 3. Modelo REST2 with some modifications


## Calculation of optical mass as per SMARTS.  character string
## indicating the process (between 'Ray' (Rayleigh dispersion), 'O3':
## absorption of ozone, 'NO2' absoption of NO2, 'Mix': absorption of
## gases uniformly mixed (02 and C02), Wva: absorption of water vapor,
## 'Aer' dispersion of aerosols, 'Kas': formula of Kaster for air
## optical mass

opticalMass <- function(theta, process){
  ## Matrix of coefficients
  At <- data.frame(
    Ray = c(4.5665e-1, 0.07, 96.4836, -1.697),
    O3 = c(2.6845e2, 0.5, 115.42, -3.2922),
    NO2 = c(6.023e2, 0.5, 117.96, -3.4536),
    Mix = c(4.5665e-1, 0.07, 96.4836, -1.697),
    Wva = c(0.10648, 0.11423, 93.781, -1.9203),
    Aer = c(0.16851, 0.18198, 95.318, -1.9542),
    Kas = c(0.50572, 0, 96.07995, -1.6364)
  )
  stopifnot(process %in% names(At))
  ## Choose the adequate variable
  A <- At[[process]]
  ## Sun zenith (theta) in degrees
  cosThzS <- cos(d2r(theta))
  m <- (cosThzS + A[1]*theta^A[2]*(A[3] - theta)^A[4])^-1
  return(m)
}

## Estimate pressure from elev and lat
press <- function(elev, lat){
  ## altitude in km
  zk <- elev/1000
  if (abs(abs(lat) - 45) >= 1e-3){
    PCOR1 <- 0.993 + 2.0783e-04*zk - 1.1589e-06*lat^2
    PCOR2 <- 8.855e-03 - 1.5236e-04*zk - 9.2907e-07*lat^2
    pcor <- PCOR1 + zk*PCOR2
  } else {
    pcor <- 1
  }
  ##  Pressure correction as function of altitude in km 
  if (zk > 10) {
    p <- exp((5.5797 - 0.14395*(zk - 10))/(1 - 0.0006335*(zk - 10)))
  } else {
    p <- 1013.25*pcor*exp(0.00177 - 0.11963*zk - 0.00136*zk^2)
  }
  p
}

# Gueymard, C.A., 2008. REST2: High-performance solar radiation model for cloudless-sky irradiance, illuminance,
# and photosynthetically active radiation – Validation with a benchmark dataset. Solar Energy 82, 272-285.

## Model parameters
# beta: Angstrom turbidity
# alpha1: alpha exponent in band 1, 290-700 nm
# alpha2: alpha exponent in band 2, 700-4000 nm
# pwc: precipitable water content (cm)
# elev: altitude from sea level (m)
# rhog: ground albedo 
# uo: ozone column (default 0.35 atm-cm)
# un: NO2 column    (default 0.0002 attm-cm)

## data is a `zoo` with intradaily values (produced by as.zooI and
## calcSol, or by fSolI). 

REST2 <- function(data, loc){
    ## Sun zenith angle
    cosThzS <- data$cosThzS
    ## in degrees
    theta <- r2d(acos(cosThzS))
  
    beta <- data$beta
    alpha1 <- data$alpha1
    alpha2 <- data$alpha2
    pwc <- data$pwc
    rhog <- data$rhog
    uo <- data$uo
    un <- data$un
    
    elev <- loc$elev
    lat <- loc$lat
    ## altitude in km
    zk <- elev/1000
    alati <- abs(lat)
    pcor <- 1
    if (abs(abs(lat) - 45) >= 1e-3){
        PCOR1 <- 0.993 + 2.0783e-04*alati - 1.1589e-06*lat^2
        PCOR2 <- 8.855e-03 - 1.5236e-04*alati - 9.2907e-07*lat^2
        pcor <- PCOR1 + zk*PCOR2}
    ## Pressure correction as function of altitude in km 
    if (zk > 10) {p <- exp((5.5797 - 0.14395*(zk - 10))/(1 - 0.0006335*(zk - 10)))}
    if (zk < 10) {p <- 1013.25*pcor*exp(0.00177 - 0.11963*zk - 0.00136*zk^2)}
    
    
    ## Retrieve or estimate pressure
    if ('p' %in% names(data)) {
        p <- data$p
    } else {
        p <- with(loc, press(elev, lat))
    }
    ##  relative air mass Rayleigh
    ## pres <- p/1013.25
    am <- p/1013.25 * opticalMass(theta, 'Ray')
    mo <- opticalMass(theta, 'O3')
    mw <- opticalMass(theta, 'Wva')
    ma <- opticalMass(theta, 'Aer')
    ##  extra-atmospheric irrradiances
    Bo0n1 <- 635.4
    Bo0n2 <- 709.7
    ##  Single scattering albedo
    omega1 <- 0.95
    omega2 <- 0.90
    
### Band 1 Calculation 290-700nm (Mostly Rayleigh and Mie
### scattering)
    
    ##Rayleigh
    T_r1 <- (1 + 1.8169*am - 0.033454*am^2)/(1 + 2.063*am + 0.31978*am^2)
    ##  Uniformly mixed gases
    T_g1 <- (1 + 0.95885*am + 0.012871*am^2)/(1 + 0.96321*am + 0.015455*am^2)
    ##   Ozone calculation
    f1 <- uo*(10.979 - 8.5421*uo)/(1 + 2.0115*uo + 40.189*uo^2)
    f2 <- uo*(-0.027589 - 0.005138*uo)/(1 - 2.4857*uo + 13.942*uo^2)
    f3 <- uo*(10.995 - 5.5001*uo)/(1 + 1.6784*uo + 42.406*uo^2)
    ##  Ozone transmittance
    T_o31 <- (1 + f1*mo + f2*mo^2)/(1 + f3*mo)
    uoLow <- (uo< 0.01)
    T_o31[uoLow] <- 1
    ##  NO2 Nitrogen dioxide
    g1 <- (0.17499 + 41.654*un - 2146.4*un^2)/(1 + 22295*un^2)
    g2 <- un*(-1.2134 + 59.324*un)/(1 + 8847.8*un^2)
    g3 <- (0.17499 + 61.658*un + 9196.4*un^2)/(1 + 74109*un^2)
    ##  nitrogen dioxide transmittance
    T_n1 <- (1 + g1*mw + g2*mw^2)/(1 + g3*mw)
    T_n1[T_n1>1] <- 1
                                        #    Tn1prima <- pmin(1, (1 + g1*1.66 + g2*1.66^2)/(1 + g3*1.66)) ##  todavia no
    ##  w vapor
    h1 <- pwc*(0.065445 + 0.00029901*pwc)/(1 + 1.2728*pwc)
    h2 <- pwc*(0.065687 + 0.0013218*pwc)/(1 + 1.2008*pwc)
    ##  w vapor transmittance
    T_w1 <- (1 + h1*mw)/(1 + h2*mw)
                                        # Tw1prima <- (1 + h1*1.66)/(1 + h2*1.66) ##  todavia no
    
  ## ##  Band 2 Calculation 700-4000 nm
    ##  Mostly gaseous absorption are relevant
    T_r2 <- (1 - 0.010394*am)/(1 - 0.00011042*am^2)
    T_g2 <- (1 + 0.27284*am - 0.00063699*am^2)/(1 + 0.30306*am)
    T_o32 <- 1
    T_n2 <- 1
    T_n2prima <- 1
    ##  Water vapor
    c1 <- pwc*(19.566 - 1.6506*pwc + 1.0672*pwc^2)/(1 + 5.4248*pwc + 1.6005*pwc^2)
    c2 <- pwc*(0.50158 - 0.14732*pwc + 0.047584*pwc^2)/(1 + 1.1811*pwc + 1.0699*pwc^2)
    c3 <- pwc*(21.286 - 0.39232*pwc + 1.2692*pwc^2)/(1 + 4.8318*pwc + 1.412*pwc^2)
    c4 <- pwc*(0.70992 - 0.23155*pwc + 0.096514*pwc^2)/(1 + 0.44907*pwc + 0.75425*pwc^2)
    ##  Water vapor transmittance
    T_w2 <- (1 + c1*mw + c2*mw^2)/(1 + c3*mw + c4*mw^2)
    T_w2prima <- (1 + c1*1.66 + c2*1.66^2)/(1 + c3*1.66 + c4*1.66^2) 
    pwcNeg <- (pwc < 0)
    T_w2prima[pwcNeg] <- 1
    T_w2[pwcNeg] < -1
    ## ##  Aerosols
    beta1 <- beta2 <- beta
    D0Alpha <- (alpha1 - alpha2)> 1e-3
    beta1[D0Alpha] <- beta1[D0Alpha] * 0.7^(alpha1 - alpha2)[D0Alpha]
    
    ua1 <- log(1 + ma*beta1)
    ua2 <- log(1 + ma*beta2)
    ##  Effective aerosol wavelength, band 1
    ##  as per fortran code REST2 8.3 version
    d0 <- 0.57664 - 0.024743*alpha1
    d1 <- (0.093942 - 0.2269*alpha1 + 0.12848*alpha1^2)/(1 + 0.6418*alpha1)
    d2 <- (-0.093819 + 0.36668*alpha1 - 0.12775*alpha1^2)/(1 - 0.11651*alpha1)
    d3 <- alpha1*(0.15232 - 0.087214*alpha1 + 0.012664*alpha1^2)/(1 - 0.90454*alpha1 + 0.26167*alpha1^2)
    
    a1Low <- (alpha1 <= 1.3)
    d0[a1Low] <- 0.544474
    d1[a1Low] <- 0.00877874
    d2[a1Low] <- 0.196771
    d3[a1Low] <- 0.294559
    
    lambda1 <- (d0 + d1*ua1 + d2*ua1^2)/(1 + d3*ua1^2)
    ##  Limits to lambda1
    lambda1[lambda1 < 0.3] <- 0.3
    lambda1[lambda1 > 0.65] <- 0.65
    
    taua1 <- beta1*lambda1^(-alpha1)
    ##  Aerosol transmittance
    T_a1 <- exp(-ma*taua1)
    T_as1 <- exp(-ma*omega1*taua1)
    ##  Effective aerosol wavelength, band 2
    ##  Cambiado segun codigo fuente fortran version REST2 8.3
    e0 <- (1.183 - 0.022989*alpha2 + 0.020829*alpha2^2)/(1 + 0.11133*alpha2)
    e1 <- (-0.50003 - 0.18329*alpha2 + 0.23835*alpha2^2)/(1 + 1.6756*alpha2)
    e2 <- (-0.50001 + 1.1414*alpha2 + 0.0083589*alpha2^2)/(1 + 11.168*alpha2)
    e3 <- (-0.70003 - 0.73587*alpha2 + 0.51509*alpha2^2)/(1 + 4.7665*alpha2)
    a2Low <- (alpha2 <= 1.3)
    e0[a2Low] <- 1.038076
    e1[a2Low] <- -0.105559
    e2[a2Low] <- 0.0643067
    e3[a2Low] <- -0.109243
    lambda2 <- (e0 + e1*ua2 + e2*ua2^2)/(1 + e3*ua2^2)
    ##  Limits to lambda2
    lambda2[lambda2 < 0.75] <- 0.75
    lambda2[lambda2 > 1.5] <- 1.5
    
    taua2 <- beta2*lambda2^(-alpha2)
    ##  Aerosol transmittance
    T_a2 <- exp(-ma*taua2)
    T_as2 <- exp(-ma*omega2*taua2)
    ##  Aerosol scattering correction factor , band 1
    g0 <- (3.7155 + 0.368*ma + 0.036294*ma^2)/(1 + 0.0009391*ma^2)
    g1 <- (-0.164 - 0.72567*ma + 0.20701*ma^2)/(1 + 0.0019012*ma^2)
    g2 <- (-0.052288 + 0.31902*ma + 0.17871*ma^2)/(1 + 0.0069592*ma^2)
    F1 <- (g0 + g1*taua1)/(1 + g2*taua1)
    ##  Aerosol scattering correction factor , band 2
    h0 <- (3.4352 + 0.65267*ma + 0.00034328*ma^2)/(1 + 0.034388*ma^1.5)
    h1 <- (1.231 - 1.63853*ma + 0.20667*ma^2)/(1 + 0.1451*ma^1.5)
    h2 <- (0.8889 - 0.55063*ma + 0.50152*ma^2)/(1 + 0.14865*ma^1.5)
    F2 <- (h0 + h1*taua2)/(1 + h2*taua2)
    ##  Sky albedo for backscattered radiation
    rhoa1 <- (0.13363 + 0.00077358*alpha1 + beta1*(0.37567 + 0.22946*alpha1)/
                  (1 - 0.10832*alpha1))/(1 + beta1*(0.84057 + 0.68683*alpha1)/(1 - 0.08158*alpha1))
    rhoa1[rhoa1 < 0.01] <- 0.01
    rhoa2 <- (0.010191 + 0.00085547*alpha2 + beta2*(0.14618 + 0.062758*alpha2)/
                  (1 - 0.19402*alpha2))/(1 + beta2*(0.58101 + 0.17426*alpha2)/(1 - 0.17586*alpha2))
    rhoa2[rhoa2 < 0.005] <- 0.005
    ##  Correction Sep 2010
    ## bet22 <- beta*beta
    AM6 <- am
    AM6[am > 6] <- 6
    R1a0 <- (0.02255 + 0.27375*beta - 0.41817*beta^2)/(1 + 4.5555*beta)
    R1a1 <- (.044773 + 0.15345*beta - 0.30157*beta^2)/(1 + 3.4223*beta)
    R1corr <- (1 + R1a0*(AM6 - 1))/(1 + R1a1*AM6)
    rhos1 <- rhoa1*R1corr
    ##  Reflected radiation is function of ground albedo rhog
    R1 <- rhoa1*rhog
    R2 <- rhoa2*rhog
###  Beam irrradiance
    Bn1 <- T_r1 * T_g1 * T_o31 * T_n1 * T_w1 * T_a1 * Bo0n1
    Bn2 <- T_r2 * T_g2 * T_o32 * T_n2 * T_w2 * T_a2 * Bo0n2
    Bn1[Bn1<0] <- 0
    Bn2[Bn2<0] <- 0
    Bn <- Bn1 + Bn2
    B01 <- Bn1 * cosThzS
    B02 <- Bn2 * cosThzS
    B0 <- B01 + B02
###  D0fuse
    BR1 <- 0.4625*(1.0024 - 0.0055808*am + 0.000051487*am^2)
    BR2 <- 0.588
    Ba <- 1 - exp(-0.6931 - 1.8326*cosThzS)
    Bo01 <- Bo0n1 * cosThzS
    Bo02 <- Bo0n2 * cosThzS
    
    betama <- ma*beta2
    Ratio <- 4.7488 - 0.84236*betama
    idxBM2 <- which(betama < 0.2)
    Ratio[idxBM2] <- 4.58033 + exp((0.606 - 2.9998*betama[idxBM2]^0.57)/
                                       (0.2 - betama[idxBM2]))
    EdpR1 <- BR1 * T_g1 * T_o31 * T_n1 * T_w1 * (1 - T_r1) * (T_a1^0.55) * Bo01
    EdpA1 <- Ba * T_n1 * T_w1 * (1 - T_as1^0.15) * T_g1 * T_o31 * (T_r1*Ratio) * Bo01
    EdpR1[EdpR1<0.2] <- 0.2 
    EdpA1[EdpA1<0.2] <- 0.2
    D01 <- EdpR1 + EdpA1
    EdpR2 <- BR2 * T_g2 * T_o32 * T_n2 * T_w2 * (1 - T_r2) * (T_a2^0.55) * Bo02
    EdpA2 <- Ba * T_n2 * T_w2 * (1 - T_as2^0.15) * T_g2 * T_o32 * (T_r2*Ratio)*Bo02
    EdpR2[EdpR2<0.1] <- 0.1
    EdpA2[EdpA2<0] <- 0    
    D02 <- EdpR2 + EdpA2
    ##  Adding the contribution of the reflected
    D01 <- D01 + R1 * (B01 + D01) / (1 - R1) 
    D02 <- D02 + R2 * (B02 + D02) / (1 - R2) 
    D0 <- D01 + D02
    ##  Global
    G01 <- B01 + D01
    G02 <- B02 + D02
    G0 <- G01 +  G02
    
    night <- (theta >= 90)
    B0[night] <- 0
    Bn[night] <- 0
    D0[night] <- 0
    G0[night] <- 0
    
    z <- cbind(G0, D0, B0, Bn)
    coredata(z)[is.na(z)] <- 0
    z
}

## 4. Modelo de Paulescu y Z. Schlett de 2003 CHECKED with original
## M. Paulescu and Z. Schlett. A simplified but accurate spectral solar irradiance model. 
## Theor. Appl. Climatol. 75, 203–212 (2003) DOI 10.1007/s00704-003-0731-y

PS2003 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  #p_sat <- data$p_sat
  RH <- data$RH
  pwc <- data$pwc
  p <- data$p
  uo <- data$uo
  beta <- data$beta
  k <- data$k # default a 0.5
  
  elev <- loc$elev
  
  # standard air mass using Kasten, F., Young, A.T., 1989. Revised optical air mass tables 
  # and approximation formula. Applied Optics 28(2), 4735-4738.
  m <- (1 - 1e-4*elev)/(cosThzS + 0.50572*(96.07995 - theta)^-1.6364)
  # water vapour column content (g/cm^2), Hann, J., 1901. Lehrbuch der Meteorologie, 225. Leipzig, Tauchnitz.
  #pwc <- 0.25*p_sat*RH
  T_r <- exp(-m*(p/1013)*(0.709 + 0.0013*m*(p/1013) - 0.5856*(m*(p/1013))^0.058))
  T_o3 <- exp(-m*uo*(0.0184 - 0.0004*m*uo + 0.022*((m*uo)^-0.66)))
  T_w <- exp(-m*pwc*(-0.002 + (1.67e-5)*m*pwc + 0.094*((m*pwc)^-0.693)))
  T_g <- exp(-m*(-5.4*10^-5 - (3.8e-6)*m + 0.0099*m^(-0.62)))
  T_a <- exp(-m*beta*(1.053 - 0.083*m*beta + 0.3345*((m*beta)^-0.668)))
  B0 <- Bo0 * T_r * T_a * T_w * T_o3 * T_g
  Bn <- B0/cosThzS
  D0 <- k * Bo0 * (1 - T_r*T_a) * T_w * T_o3 * T_g
  G0 <- B0 + D0

  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}

## 5. ASHRAE model Checked with Badescu

## ASHRAE, 1972. Handbook of Fundamentals. American Society of Heating, Refrigerating
## and Air-Conditioning Engineers, Inc., Atlanta, GA.


# obtained from Gueymard, C., 2012. Clear-sky irradiance predictions for solar resource mapping and
# large-scale applications: Improved validation methodology and
# detailed performance analysis of 18 broadband radiative models. Solar Energy 86, 2145–2169.

ASHRAE1972 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  month <- month(index(data))
  
  A <- c(1230, 1215, 1186, 1136, 1104, 1088, 1085, 1107, 1151, 1192, 1221, 1233)
  B <- c(0.142, 0.144, 0.156, 0.180, 0.196, 0.205, 0.207, 0.201, 0.177, 0.160, 0.149, 0.142)
  D <- c(0.058, 0.060, 0.071, 0.097, 0.121, 0.134, 0.136, 0.122, 0.092, 0.073, 0.063, 0.057)
  
  Bn <- A[month]*exp(-B[month]/cosThzS)
  D0 <- D[month]*Bn
  B0 <- Bn*cosThzS
  G0 <- B0 + D0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}



## 6. HLJ model CHECKED  with original

## HLJ results from the popular combination of the broadband Hottel model for direct transmittance (Hottel, 1976)
## and of Liu & Jordan’s Difuse transmittance expression (Hottel, 1976; Liu and Jordan, 1960).

# Hottel, H.C., 1976. A simple model for estimating the transmittance of direct solar radiation 
# through clear atmospheres. Solar Energy 18, 129–134.

# Liu, B.Y.H., Jordan, R.C., 1960. The interrelationship and characteristic distribution of direct, D0fuse and
# total solar radiation. Solar Energy 4, 1-19.

# obtained from Gueymard, C., 2012. Clear-sky irradiance predictions for solar resource mapping and
# large-scale applications: Improved validation methodology and
# detailed performance analysis of 18 broadband radiative models. Solar Energy 86, 2145–2169.


HLJ <- function(data, 
                loc = list(elev = NULL)){
  
  cosThzS <- data$cosThzS
  Bo0 <- data$Bo0
  
  elev <- loc$elev

  k0 <- 0.4327 - 0.00821*(((6000 - elev)/1000)^2)
  k1 <- 0.5055 + 0.00595*(((6500 - elev)/1000)^2)
  k2 <- 0.2711 + 0.01858*(((2500 - elev)/1000)^2)
  
  Bn <- (Bo0/cosThzS)*(k0 + k1*exp(-k2/cosThzS))
  B0 <- Bn*cosThzS
  D0 <- Bo0*(0.2710 - 0.2939*(k0 + k1*exp(-k2/cosThzS)))
  G0 <- B0 + D0

  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}


## 7. Kumar 1997 CHECKED with original

# Kumar, L, Skidmore AK and Knowles E 1997: Modelling topographic variation in solar radiation in
# a GIS environment. International Journal of Geographical Information Sciences 11(5), 475-497.

# obtained from Gueymard, C., 2012. Clear-sky irradiance predictions for solar resource mapping and
# large-scale applications: Improved validation methodology and
# detailed performance analysis of 18 broadband radiative models. Solar Energy 86, 2145–2169.

KSK1997 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  p <- data$p
  
  m <- ((1229 + (614*cosThzS)^2)^0.5 - 614*cosThzS)*p/1013.25
  T_b <- 0.56*(exp(-0.65*m) + exp(-0.095*m))
  Bn <- T_b*Bo0/cosThzS
  B0 <- Bn*cosThzS
  D0 <- Bo0*(0.271 - 0.294*T_b)
  G0 <- B0 + D0
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}

## 8. Fu and Rich model CHECKED with original and Gueymard

## Fu, P., Rich, P., 1999. Design and implementation of the Solar Analyst: an
## ArcView extension for modeling solar radiation at landscape scales. In:
## Proceedings of the 19th Annual ESRI User Conference.

# Tb bulk atmospheric transmittance, with a recommended value of 0.5
# P is the diffuse fraction of global normal irradiance, with a recommended value of 0.3

# obtained from Gueymard, C., 2012. Clear-sky irradiance predictions for solar resource mapping and
# large-scale applications: Improved validation methodology and
# detailed performance analysis of 18 broadband radiative models. Solar Energy 86, 2145–2169.

FR1999 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  Bo0 <- data$Bo0
  
  elev <- loc$elev

  k <- exp(-0.000118*elev - 1.638e-9*(elev^2))/cosThzS
  Bn <- Bo0*(0.5^k)
  B0 <- Bn*cosThzS
  D0 <- Bn*cosThzS*(0.3/(1 - 0.3))
  G0 <- B0 + D0
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}

## 9. Heliosat-1 model CHECKED with Gueymard

# obtained from Gueymard, C., 2012. Clear-sky irradiance predictions for solar resource mapping and
# large-scale applications: Improved validation methodology and
# detailed performance analysis of 18 broadband radiative models. Solar Energy 86, 2145–2169.

HE1 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  TL <- data$TL
  p <- data$p
  
  elev <- loc$elev

  TLAM2 <- TL/0.8662
  elev_c <- exp(-elev/8434.5)
  m <- elev_c/(cosThzS + 0.50572*((90 - theta) + 6.07995)^-1.6364)
  am <- m*p/1013
  foo <- 1/(6.62960 + 1.7513*am - 0.1202*am^2 + 0.0065*am^3 - 0.00013*am^4)
  Bn <- (Bo0/cosThzS)*exp(-am*foo*TL) 
  B0 <- Bn*cosThzS
  D0 <- (Bo0/cosThzS)*(0.0065 + (0.0646*TLAM2 - 0.045)*cosThzS - (0.0327*TLAM2 - 0.014)*((cosThzS)^2))
  G0 <- B0 + D0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}

## 10. HE- 2 model checked with Badescu et al.


# obtained from Gueymard, C., 2012. Clear-sky irradiance predictions for solar resource mapping and
# large-scale applications: Improved validation methodology and
# detailed performance analysis of 18 broadband radiative models. Solar Energy 86, 2145–2169.

# Christelle Rigollier, Mireille Lefevre, Sylvain Cros, Lucien Wald. Heliosat 2: an improved
# method for the mapping of the solar radiation from Meteosat imagery. 2002 EUMETSAT
# Meteorological Satellite Conference, Sep 2002, Dublin, Ireland. EUMETSAT, Darmstadt, Germany, 585-592

HE2 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  TL <- data$TL
  
  elev <- loc$elev

  TLAM2 <- TL/0.8662
  elev_c <- exp(-elev/8434.5)
  m <- elev_c/(cosThzS + 0.50572*((90 - theta) + 6.07995 )^-1.6364)
  m_c <- 1/(cosThzS + 0.50572*((90 - theta) + 6.07995 )^-1.6364)

  foo1 <- 1.248174 - 0.011997*m_c + 0.00037*m_c^2
  foo2 <- 1.68219 - 0.03059*m_c + 0.00089*m_c^2
  foo3 <- (elev_c - 0.75)^2
  foo4<- (elev_c - 1.0)^2
  if(elev_c > 0.75){if(elev_c >= 0.99){p_c <- rep(1, length(foo1))}else{p_c <- (foo3 + foo4*foo1)/(foo3 + foo4)}}else{
  p_c <- (foo3*foo2 + foo4*foo1)/(foo3 + foo4)}
  foo5 <- (1/(6.62960 + 1.7513*m - 0.1202*m^2 + 0.0065*m^3 - 0.00013*m^4))/p_c

  Bn <- (Bo0/cosThzS)*exp(-foo5*0.8662*TLAM2*m)
  B0 <- Bn*cosThzS

  TLAM2_c <- TLAM2*elev_c

  T_d <- -1.5843e-2 + 3.0543e-2*TLAM2_c + 3.797e-4*TLAM2_c^2
  a <- matrix(c(2.6463e-1, -6.1581e-2, 3.1408e-3, 2.0402, 1.8945e-2,
              -1.1161e-2, -1.3025, 3.9231e-2, 8.5079e-3), nrow=3, ncol=3, byrow=T)
  A0 <- as.numeric(a[1,1] + a[1,2]*TLAM2_c + a[1,3]*TLAM2^2)
  A1 <- as.numeric(a[2,1] + a[2,2]*TLAM2_c + a[2,3]*TLAM2^2)
  A2 <- as.numeric(a[3,1] + a[3,2]*TLAM2_c + a[3,3]*TLAM2^2)
  
  A0[((T_d*A0) < 2e-3)] <- 2e-3/T_d[((T_d*A0) < 2e-3)]
  
  D0 <- Bo0*T_d*(A0 + A1*cosThzS + A2*cosThzS^2)
  G0 <- B0 + D0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}

## 11. Adnot et al. 1979 CHECKED with Dazhi et al.
## direct: Adnot J, Bourges B, Campana D, Gicquel R (1979) Utilisation des courbes de frequence cumulees pour le calcul des
## installation solaires. In: Lestienne R (ed) Analise Statistique des Processus Meteorologiques Appliquee a l’Energie
## Solaire. Paris: CNRS, pp 9–40

## Difuse: Schulze, R.E., 1976. A physically based method of estimating solar radiation from suncards.
# Agricultural Meteorology 16(1), 85-101.

## obtained from:
# Dazhi, Y., Jirutitijaroen, P., Walsh, W.M., 2012. The Estimation of Clear Sky Global Horizontal Irradiance at
# the Equator. Energy Procedia 25, 141-148.


ABCGS <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  
  G0 <- 951.3*(cosThzS)^1.15 # Badescu et al uses 951.39 instead of 951.3
  D0 <- 94.23*(cosThzS^0.5)
  B0 <- G0 - D0
  B0[B0<0] <- 0
  Bn <- B0/cosThzS
  
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}
  

## 12. Bugler 1977 model CHECKED with original

## Bugler, J.W., 1977. The determination of hourly insolation on an inclined plane using a 
## diffuse irradiance model based on hourly measured global horizontal insolation. Solar Energy 19 (5), 477-491.

B1977 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))

  D0 <- 16*sqrt(90 - theta) - 0.4*(90 - theta)
  
  night <- (theta >= 90)
  D0[night] <- 0
  
  G0 <- B0 <- Bn <- zoo(rep(NA, length(D0)), index(D0))
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}
  

## 13. Samimi model
## Samimi, J. (1994) Estimation of Height-Dependent Solar Irradiation and Application to
## the Solar Climate of Iran. Solar Energy, 52, 401-409. 


S1994 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  elev <- loc$elev

  foo <- 1 - exp((-36/pi)*((pi/2) - d2r(theta)))
  B0 <- Bo0*((1 - 0.14*elev)*exp(-0.357/(cosThzS^0.678)) + 0.14*elev*foo)*cosThzS # REVISAR
  G0 <- 1.1*B0
  D0 <- G0 - B0
  Bn <- B0/cosThzS
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}

## 14. DPP model CHECKED with original 
## Daneshyar, M., 1978. Solar radiation statistics for Iran. Solar Energy 21(4), 345–349.
## G.W. Paltridge, D. Proctor. (1976) Monthly mean solar radiation statistics for Australia. Solar Energy 18:3, 235-243. 


DPP <- function(data, loc){
    
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
      
  B0 <- 950*(1 - exp(-0.075*(90 - theta)))*cosThzS
  D0 <- 14.29 + 21.04*((pi/2) - d2r(theta))
  G0 <- B0 + D0
  Bn <- B0/cosThzS
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}
  
 

## 15. Yang model CHECKED with original
## Yang, K., Huang, G.W., Tamai, N., 2001. A hybrid model for estimating global solar radiation.
## Solar Enery 70, 13-22.

YH2001 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  y <- unique(year(index(data)))
  julianday <- as.numeric(julian(as.Date(index(data)), origin = as.Date(paste0(y,'-01-01')))) + 1
  
  #p_sat <- data$p_sat
  RH <- data$RH
  p <- data$p
  m <- data$m
  pwc <- data$pwc
  uo <- data$uo
  beta <- data$beta
  
  elev <- loc$elev
  lat <- loc$lat
  
  #m[is.na(m)] <- (1 - 0.0001*elev)/(cosThzS + 0.15*(57.296*(90 - theta) + 3.885)^-1.253)
  # if beta not available
  #beta[is.na(beta)] <-  rep((0.025 + 0.1*cos(d2r(lat)))*exp(-0.7*elev/1000) + 0.04, length(theta)) # +/- 0.02-0.06
  # if uo not available
  #d <- julianday
  #d[julianday>300] <- julianday - 366
  #uo <- 0.44 - 0.16*(((lat - 80)/60)^2 + ((d -120)/(263 - lat))^2)^0.5
  # if w not available
  #pwc[is.na(pwc)] <- 0.25*p_sat*RH # using Hann, J., 1901. Lehrbuch der Meteorologie, 225. Leipzig, Tauchnitz.
  
  foo1 <- 0.547 + 0.014*(m*p/1013) - 0.00038*(m*p/1013)^2 + (4.6*10^-6)*(m*p/1013)^3
  foo2 <- 0.6777 + 0.1464*(m*beta) - 0.00626*(m*beta)^2
  T_r <- exp(-0.008735*m*((foo1)^-4.08)*p/1013)
  T_o3 <- exp(-uo*m*0.0365*((m*uo)^-0.2865))
  T_w <- exp(-(-log(0.909 - 0.036*log(m*pwc))))
  T_g <- exp(-0.0117*(m^0.3139))
  T_a <- exp(-beta*m*(foo2^-1.3))
  Bn <- (Bo0/cosThzS)*(T_o3*T_w*T_g*T_r*T_a - 0.013)
  B0 <- Bn*cosThzS
  D0 <- Bo0*(T_o3*T_g*T_w*(1 - T_a*T_r) + 0.013)
  G0 <- B0 + D0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}


## 16. MAC model CHECKED with Badescu et al
## J. A. Davies, D. C. McKay, G. Luciani, and M. Abdel-Wahab, Validation of models for estimating solar radiation
## on horizontal surfaces, IEA Task IX Final Report, Vol. 1, Atmospheric Environment Service, Downsview (1988).

# k default to 0.95 due to 
# Davies, J.A., Schertzer, W., Nunez, M., 1975. Estimating global solar radiation. Boundary-Layer Meteorology 9, 33-52

# model obtained from Gueymard, C., 1993. Critical analysis and performance assessment
# of clear sky solar irradiance models using theoretical and measured data. Solar Energy 51(2), 121-138.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.


MAC <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  p <- data$p
  pwc <- data$pwc
  k <- data$k # default 0.95
  m <- data$m
  rhog <- data$rhog
  
  m[is.na(m)] <- 35/(1 + 1224*(cosThzS^2)^0.5)
  am <- m*p/1013.25
  foo <- 8.688237*am^(0.0279286*log(am) - 0.806955)
  T_r <- foo/(1 + foo)
  T_o3 <- 1 - (((0.002118*3.5*m)/(1 + 0.0042*3.5*m + (3.23*10^-6)*(3.5*m)^2)) + 
                 ((0.1082*3.5*m)/((1 + 13.86*(3.5*m))^0.805)) + (0.00658*3.5*m/(1 + (10.36*3.5*m)^3)))
  T_a <- k^am 
  aw <- 0.29*pwc/((1 + 14.15*pwc)^0.635 + 0.5925*pwc) 
  Bn <- Bo0*(T_o3*T_r - aw)*T_a/cosThzS
  B0 <- Bn*cosThzS
  foo1 <- 0.93 - 0.21*log(m)
  foo2 <- 0.93 - 0.21*log(1.66)
  rhoa <- 0.0685 + (1 - T_a)*0.98*(1 - foo2)
  Dif1 <- Bo0*(0.5*T_o3*foo2*(1 - T_r) + 0.98*foo1*(T_o3*T_r - aw)*(1 - T_a))
  Dif2 <- rhog*rhoa*(Bn + Dif1)/(1 - rhoa*rhog)
  D0 <- Dif1 + Dif2
  G0 <- D0 + B0
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}

## 17. WMO 1 Checked with Rigollier et al.
## World Meteorological Organization Document 557 (1981), pp 120

# obtained from Rigollier, G., Wald, L., 2000. Selecting a clear-sky model to accurately map solar radiation
# from satellite images. Remote sensing in the 21st century. Economic and Environmental Applications, Casanova (Ed),
# Rotterdam, Nederlands.

WMO1 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  G0 <- 0.95*Bo0*cosThzS/(1 + 0.2/cosThzS)
  D0 <- B0 <- Bn <- zoo(rep(NA, length(G0)), index(G0))

  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}

## 18. Bourges model CHECKED with Rigollier

## Bourges, G., 1979. Reconstitution des courbes de fréquence cumulées de l'irradiation solaire globale horaire
## recue par une surface plane. Report CEE 295-77-ESF f Centre d'Enérgetique de l'Ecole Nationale Supérieure des Mines 
## de Paris, tome II. Paris, France.

# obtained from Rigollier, G., Wald, L., 2000. Selecting a clear-sky model to accurately map solar radiation
# from satellite images. Remote sensing in the 21st century. Economic and Environmental Applications, Casanova (Ed),
# Rotterdam, Nederlands.

B1979 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  day <- as.Date(index(data))
  
  G0 <- 0.7*Bo0*(cosThzS^1.15)
  B0 <- Bn <- D0 <- zoo(rep(NA, length(G0)), index(G0))

  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}


## 19. Perrin de Brichambaut and Vauge (1982) CHECKED with Rigollier
## Perrin de Brichambaut, C. Vauge, C., 1982. Le gisement solaire: Evaluation de la ressource energetique (Paris: 
## Technique et documentation).

# obtained from Rigollier, G., Wald, L., 2000. Selecting a clear-sky model to accurately map solar radiation
# from satellite images. Remote sensing in the 21st century. Economic and Environmental Applications, Casanova (Ed),
# Rotterdam, Nederlands.

PV1982 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0 
  
  G0 <- 0.81*Bo0*((cosThzS)^1.15)
  
  D0 <- B0 <- Bn <- zoo(rep(NA, length(G0)), index(G0))
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}

## 20. WMO2 CHECKED with Rigollier et al
## World Meteorological Organization Document 557 (1981), pp 124
# This model is suitable for elev higher than 20-30º.

# obtained from Rigollier, G., Wald, L., 2000. Selecting a clear-sky model to accurately map solar radiation
# from satellite images. Remote sensing in the 21st century. Economic and Environmental Applications, Casanova (Ed),
# Rotterdam, Nederlands.

WMO2 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  
  TL <- data$TL
  
  lat <- loc$lat
  
  eo <- fSolD(lat, as.Date(index(data)))$eo
  eo <- lapply(c(1:length(eo)), function(x){
    rep(as.numeric(eo[x]), length(which(as.Date(index(data)) == as.Date(index(eo)[x]))))}); eo <- do.call(c, eo)

  G0 <- eo*(1297 - 57*TL)*(cosThzS)^((36 + TL/0.8662)/33)
  D0 <- B0 <- Bn <- zoo(rep(NA, length(G0)), index(G0))
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}

## 21. DMLA model CHECKED with Gueymard, 1993
## J. A. Davies, D. C. McKay, G. Luciani, and M. Abdel-Wahab, Validation of models for estimating solar radiation
## on horizontal surfaces, IEA Task IX Final Report, Vol. 1, Atmospheric Environment Service, Downsview (1988).

# similar to MAC -- similar to Joseffson
# model obtained from Gueymard, C., 1993. Critical analysis and performance assessment
# of clear sky solar irradiance models using theoretical and measured data. Solar Energy 51(2), 121-138.

# k normally 0.95
  
DMLA1988 <- function(data, loc){
    
  cosThzS <- data$cosThzS
  Bo0 <- data$Bo0
  
  m <- data$m
  pwc <- data$pwc
  k <- data$k
  
  T_o3 <- 0.95545
  T_r <- 0.9768 - 0.0874*m + 0.010607552*m^2 - (8.46205e-4)*m^3 + (3.57246e-5)*m^4 - (6.0176e-7)*m^5
  T_a <- k^m  
  T_aa <- 1 - (1 - 0.75)*(1 - T_a)
  T_as <- 1 - 0.75*(1 - T_a)
  aw <- 0.0946*pwc^0.303
  Bn <- Bo0*(T_o3*T_r*T_aa*T_as - aw)/cosThzS
  Bn[Bn<0] <- 0
  B0 <- Bn*cosThzS
  # D0fuse for zero albedo
  D0 <- Bo0*(0.5*T_o3*T_aa*T_as*(1 - T_r) + 0.8832*(T_o3*T_r*T_aa - aw)*(1 - T_as))
  G0 <- D0 + B0
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}

## 22. Iqbal model A. CHECKED with Badescu et al.
# M. Iqbal, An introduction to solar radiation, Academic Press, Toronto (1983).

# similar to MAC 
# model obtained with some simplifications from Gueymard, C., 1993. Critical analysis and performance assessment
# of clear sky solar irradiance models using theoretical and measured data. Solar Energy 51(2), 121-138.


IA1983 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- data$theta
  Bo0 <- data$Bo0
  
  m <- data$m
  pwc <- data$pwc
  uo <- data$uo
  rhog <- data$rhog
  TempK <- data$Temp + 273
  p <- data$p
  k <- data$k
  alpha <- data$alpha
  beta <- data$beta
  am <- m*p/1013.25
  
  T_o3 <- 1 - (((0.002118*uo*m)/(1 + 0.0042*uo*m + (3.23e-6)*(uo*m)^2)) + 
                 ((0.1082*uo*m)/((1 + 13.86*(uo*m))^0.805)) + (0.00658*uo*m/(1 + (10.36*uo*m)^3)))
  
  w_p <- 10*pwc*m*p/1013.25
  w_p[is.na(TempK)==FALSE] <- 10*pwc*m*((p/1013.25)^0.75)*(273/TempK)^0.5
  aw <- (0.29*w_p)/((1 + 14.15*w_p)^0.635 + 0.5925*w_p)
  T_r <- 0.972 - 0.08262*am + 0.00933*am^2 - 0.00095*am^3 + 0.0000437*am^4 
  T_r[T_r>1]<-0
  T_a <- (0.12445*alpha - 0.0162) + (1.003 - 0.125*alpha)*exp(-beta*am*(1.089*alpha + 0.5123))
  T_a[is.na(k)==FALSE] <- k^m
  T_a[((is.na(k)[1]==TRUE) & (beta > 0.5))] <- 0.95^m
  T_ab <- 0.12445*alpha - 0.0162 + (1.003 - 0.125*alpha)*
    exp(-(1.66*p/1013.25)*beta*(1.089*alpha + 0.5123))
  
  rhoa <- 0.0685 + (1 - T_ab)*0.1581
  rhoa[is.na(alpha*beta)==TRUE] <- 0
  
  Bn <- (Bo0/cosThzS)*(T_o3*T_r - aw)*T_a
  Bn[Bn<0] <- 0
  B0 <- Bn*cosThzS
  D01 <- Bo0*(0.5*T_o3*T_ab*(1 - T_r) + 0.9*(0.93 - 0.21*log(m))*(T_r*T_o3 - aw)*(1 - T_a))
  D02 <- rhog*rhoa*(B0 + D01)/(1 - rhog*rhoa)
  D0 <- D01 + D02
  G0 <- B0 + D0
  G0[G0<0] <- 0
  D0[D0<0] <- 0
  B0[B0<0] <- 0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z

}


## 23. Iqbal model B. CHECKED
# M. Iqbal, An introduction to solar radiation, Academic Press, Toronto (1983).

# similar to MAC 
# model obtained with some simplifications from Gueymard, C., 1993. Critical analysis and performance assessment
# of clear sky solar irradiance models using theoretical and measured data. Solar Energy 51(2), 121-138.

IB1983 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  Bo0 <- data$Bo0

  m <- data$m
  pwc <- data$pwc
  uo <- data$uo
  rhog <- data$rhog
  TempK <- data$Temp + 273
  p <- data$p
  beta <- data$beta
  beta[beta>0.5] <- 0.5
  
  am <- m*p/1013.25
  w_p <- pwc*m*p/1013.25
  w_p[is.na(TempK)==FALSE] <- pwc*m*((p/1013.25)^0.75)*(273/TempK)^0.5
  T_r <- 0.615958 + 0.375566*exp(-0.221185*am)
  a_w <- 0.110*(w_p + 6.31e-4)^0.3 - 0.0121
  a_g <- 0.002351*(126*am + 0.0129)^0.26 - 7.5e-4 + 7.5e-3*am^0.875
  a_o3 <- 0.045*(uo*m + 8.34e-4)^0.38 - 3.1e-3
  T_as <- (-0.914 + 1.909267*exp(-0.667023*beta))^am
  a_a <- (1 - 0.95)*T_as
  Bn <- (Bo0/cosThzS)*(1 - (a_w + a_g + a_o3 + a_a))*T_r*T_as
  B0 <- Bn*cosThzS
  D01 <- Bo0*(1 - (a_w + a_g + a_o3 + a_a))*(0.5*(1 - T_r) + 0.75*(1 - T_as))
  a_o31 <- 0.045*(uo*1.66 + 8.34e-4)^0.38 - 3.1e-3
  T_as1 <- (-0.914 + 1.909267*exp(-0.667023*beta))^(1.66*p/1013.25)
  T_r1 <- 0.615958 + 0.375566*exp(-0.221185*(1.66*p/1013.25))
  Q <- B0 + D01
  a_g1 <- 0.00235*(126*(1.66*p/1013.25) + 0.0129)^0.26 - 7.5e-4 + 7.5e-3*(1.66*p/1013.25)^0.875
  T_as1 <- (-0.914 + 1.909267*exp(-0.667023*beta))^(1.66*p/1013.25)
  a_a1 <- (1 - 0.95)*T_as1
  D02 <- rhog*(B0 + D01)*(1 - (a_a1 + a_g1 + a_o31 + a_w))*(0.5*(1 - T_r1) + 0.25*(1 - T_as1))
  D0 <- D01 + D02 
  G0 <- D0 + B0
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}


## 24. Iqbal model C. CHECKED
# M. Iqbal, An introduction to solar radiation, Academic Press, Toronto (1983).

# similar to MAC 
# model obtained with some simplifications from Gueymard, C., 1993. Critical analysis and performance assessment
# of clear sky solar irradiance models using theoretical and measured data. Solar Energy 51(2), 121-138.


IC1983 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  Bo0 <- data$Bo0
  
  m <- data$m
  pwc <- data$pwc
  uo <- data$uo
  rhog <- data$rhog
  p <- data$p
  aod380 <- data$aod380
  aod500 <- data$aod500
  

  am <- m*p/1013.25
  T_r <- exp(-0.0903*(am)^0.84)*(1 + am - am^1.01)
  T_o3 <- 1 - (0.1611*uo*m)/(1 + 139.48*uo*m)^(0.3035) - 0.002715*uo*m/(1 + 0.044*uo*m + 0.0003*(uo*m)^2)
  T_g <- exp(-0.0127*(am^0.26))
  T_w <- 1 - 2.4959*pwc*m*(((1 + 79.034*pwc*m)^0.6828 + 6.385*pwc*m)^-1)
  ka <- 0.2758*aod380 + 0.35*aod500
  T_a <- exp(-(ka^0.873)*(1 + ka - ka^0.7088)*am^0.9108)
  T_aa <- 1 - (1 - 0.9)*(1 - am + am^1.06)*(1 - T_a)
  Bn <- 0.9751*(Bo0/cosThzS)*T_r*T_o3*T_g*T_w*T_a
  B0 <- Bn*cosThzS
  D01 <- 0.79*Bo0*T_o3*T_g*T_w*T_aa*0.5*(1 - T_r)/(1 - am + am^1.02)
  D02 <- 0.79*Bo0*T_o3*T_g*T_w*T_aa*0.84*(1 - T_a/T_aa)/(1 - am + am^1.02)
  rhoa <- 0.0685 + (1 - 0.84)*(1 - T_a/T_aa)
  D03 <- (B0 +  D01 + D02)*rhog*rhoa/(1 - rhog*rhoa)
  D0 <- D01 + D02 + D03
  G0 <- B0 + D0
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}




## 25. Atwater and Ball. CHECKED with original

# Atwater, M. A.; Ball, J. T., 1978. A Numerical Solar Radiation Model Based on
# Standard Meteorological Observations. Solar Energy 21, 163-170.

# obtained from Bird, R.E., Hulstrom, R.L., 1980. Direct insolation models, Solar Energy
# Research Institute (now NREL), Golden, CO, SERI/TR-335-344.
# http://www.nrel.gov/docs/legosti/old/344.pdf

AB1978 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  Bo0 <- data$Bo0
  
  pwc <- data$pwc
  p <- data$p
  alpha <- data$alpha
  
  M <- (1/35)*(((1224*(cosThzS)^2) + 1)^0.5)*p/1013
  T_d <- 1.041 - (0.15*((949e-5*p/10 + 0.051))^0.5)/M
  T_b <- 1.021 - 0.0824*((949e-5*p/10 + 0.051)/M)^0.5
  T_a <- exp(-alpha*M)
  Bn <- (Bo0/cosThzS)*(T_d*T_b - 0.077*(pwc/M)^0.3)*T_a
  B0 <- Bn*cosThzS
  B0[B0<0] <- 0
  Bn <- B0/cosThzS
  G0 <- D0 <- zoo(rep(NA, length(B0)), index(B0))

  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}

## 26. Hoyt model CHECKED with original
# Hoyt, D.V., 1978. A model for the calculation of solar global insolation. Solar Energy 21, 27–35.

# obtained from Bird, R.E., Hulstrom, R.L., 1980. Direct insolation models, Solar Energy
# Research Institute (now NREL), Golden, CO, SERI/TR-335-344.
# http://www.nrel.gov/docs/legosti/old/344.pdf

# uc pressure-corrected amount of carbon dioxide in the path [cm at standard temperature and pressure (STP) , 
# default 126 cm for air mass 1.0


H1978 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  pwc <- data$pwc
  beta <- data$beta
  m <- data$m
  p <- data$p
  uc <- data$uc
  uo <- data$uo
  am <- m*p/1013.25
 
  
  foo1 <- 0.11*(0.75*pwc + 6.31e-4)^0.3 - 0.0121
  foo2 <- 0.00235*(uc + 0.0129)^0.26 - 7.5e-4
  foo3 <- 0.045*(uo + 8.34e-4)^0.38 - 3.1e-3
  foo4 <- 7.5e-3*(m*p/1013)^0.875
  T_r <- 0.615958 + 0.375566*exp(-0.221185*am)
  s_d <- -((-0.914 + 1.909267*exp(-0.667023*beta))^m) - 1
  T_as <- (-0.914 + 1.909267*exp(-0.667023*beta))^(1.66*p/1013.25)
  foo5 <- 0.05*T_as
  Bn <- (Bo0/cosThzS)*(1 - (foo1 + foo2 + foo3 + foo4 + foo5))*T_r*T_as
  B0 <- Bn*cosThzS
  D0 <- Bo0*(1 - (foo1 + foo2 + foo3 + foo4 + foo5))*(0.5*(1 - T_r) + 0.75*(1 - T_as))
  G0 <- B0 + D0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}

## 27. Majumdar model CHECKED with original

# Majumdar, N. C. ; Mathur, B. L . ; Kaushik , S. B., 1972. Prediction of Direct
# Solar Radiation for Low Atmospheric Turbidity. Solar Energy 13, 383-394. 

# obtained from Bird, R.E., Hulstrom, R.L., 1980. Direct insolation models, Solar Energy
# Research Institute (now NREL), Golden, CO, SERI/TR-335-344.
# http://www.nrel.gov/docs/legosti/old/344.pdf

M1972 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  
  pwc <- data$pwc
  p <- data$p
  m <- data$m
  
  Bn <- 1331*((0.8644)^(m*p/1000))*(0.8507^((pwc*m)^0.25))
  B0 <- Bn*cosThzS
  
  G0 <- D0  <- zoo(rep(NA, length(B0)), index(B0))
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 28. Bird 1 Checked with original 

# Bird, R.E., Hulstrom, R.L., 1980. Direct insolation models, Solar Energy
# Research Institute (now NREL), Golden, CO, SERI/TR-335-344.
# http://www.nrel.gov/docs/legosti/old/344.pdf

# global from Badescu et al.


B1 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0

  pwc <- data$pwc
  p <- data$p
  uo <- data$uo
  aod380 <- data$aod380
  aod500 <- data$aod500
  rhog <- data$rhog

  m <- (cosThzS + 0.15*(93.885 - theta)^-1.25)^-1
  am <- m*p/1013.25
  tau_a <- 0.2758*aod380 + 0.35*aod500
  T_a <- exp((-tau_a^0.873)*(1 + tau_a - tau_a^0.7088)*m^0.9108)
  T_w <- 1 - 2.4959*pwc*m*((1 + 79.034*pwc*m)^0.6828 + 6.385*pwc*m)^-1
  T_r <- exp(-0.0903*((am)^0.84)*(1 + am - am^1.01))
  T_o3 <- 1 - 0.1611*uo*m/(1 + 139.48*uo*m)^0.3035 - 0.002715*uo*m*(1 + 0.044*uo*m + 0.0003*(uo*m)^2)^-1
  T_u <- exp(-0.0127*am^0.26)
  Bn <- (Bo0/cosThzS)*0.9662*T_r*T_o3*T_u*T_w*T_a
  B0 <- Bn*cosThzS

  T_aa <- 1 - 0.1*(1 - m + m^1.06)*(1 - T_a)
  T_as <- T_a/T_aa
  Tabs <- T_o3*T_u*T_aa*T_w
  Ts0 <- 0.79*Tabs*(0.5*(1 - T_r) + 0.84*(1 - T_as))/(1 - m + m^1.02)
  rhoa <- 0.0685 + 0.16*(1 - T_as)
  Dif0 <- Bo0*0.79*Tabs*(0.5*(1 - T_r) + 0.84*(1 - T_as))/(1 - m + m^1.02)
  G0 <- (B0 + Dif0)/(1 - rhog*rhoa)
  D0 <- G0 - B0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 29. Bird 2 CHECKED with original

# Bird, R.E., Hulstrom, R.L., 1980. Direct insolation models, Solar Energy
# Research Institute (now NREL), Golden, CO, SERI/TR-335-344.
# http://www.nrel.gov/docs/legosti/old/344.pdf

# global from Badescu


B2 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0

  pwc <- data$pwc
  p <- data$p
  uo <- data$uo
  m <- data$m
  aod380 <- data$aod380
  aod500 <- data$aod500
  rhog <- data$rhog
  
  m <- (cosThzS + 0.15*(93.885 - theta)^(-1.25))^-1
  am <- m*p/1013
  tau_a <- 0.2758*aod380 + 0.35*aod500
  T_a <- exp((-tau_a^0.873)*(1 + tau_a - tau_a^0.7088)*m^0.9108)
  aw <- 2.4959*pwc*m*((1 + 79.034*pwc*m)^0.6828 + 6.385*pwc*m)^-1
  T_r <- exp(-0.0903*(am^0.84)*(1 + am - am^1.01))
  T_m <- 1.041 - 0.15*(m*(9.368e-4*p + 0.051))^0.5
  Bn <- (Bo0/cosThzS)*0.9662*(T_m - aw)*T_a
  B0 <- Bn*cosThzS
  T_o3 <- 1 - 0.1611*uo*m/(1 + 139.48*uo*m)^0.3035 - 0.002715*uo*m*(1 + 0.044*uo*m + 0.0003*(uo*m)^2)^-1
  T_u <- exp(-0.0127*(m*p/1013)^0.26)
  T_aa <- 1 - 0.1*(1 - m + m^1.06)*(1 - T_a)
  T_as <- T_a/T_aa
  Tabs <- T_o3*T_u*T_aa*(1 - aw)
  Ts0 <- 0.79*Tabs*(0.5*(1 - T_r) + 0.84*(1 - T_as))/(1 - m + m^1.02)
  rhoa <- 0.0685 + 0.16*(1 - T_as)
  Dif0 <- Bo0*0.79*Tabs*(0.5*(1 - T_r) + 0.84*(1 - T_as))/(1 - m + m^1.02)
  G0 <- (B0 + Dif0)/(1 - rhog*rhoa)
  D0 <- G0 - B0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}

## 30. Kasten–Czeplak (KC) model (1980) CHECKED with original

# F. Kasten and G. Czeplak, "Solar and terrestrial radiation dependent on the amount and
# type of cloud," Solar Energy, vol. 24, pp. 177-189, 1980.

# obtained from Reno, M.J., Hansen, C. W., Stein, J.S., 2012. Global Horizontal Irradiance Clear Sky Models:
# Implementation and Analysis. SANDIA REPORT SAND2012-2389.

KC1980 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  
  G0 <- 910*cosThzS - 30
  D0 <- B0 <- Bn <- zoo(rep(NA, length(G0)), index(G0))
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}



## 31. Berger-Duffie model

# Badescu, 1998. Verification of some very simple clear and cloudy sky models to evaluate
# global solar irradiance. Solar Energy 61(4) 251-264.

# obtained from Reno, M.J., Hansen, C. W., Stein, J.S., 2012. Global Horizontal Irradiance Clear Sky Models:
# Implementation and Analysis. SANDIA REPORT SAND2012-2389.

BD <- function(data, loc){
  
  Bo0 <- data$Bo0
  
  G0 <- 0.7*Bo0
  D0 <- B0 <- Bn <- zoo(rep(NA, length(G0)), index(G0))
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}

## 32. Robledo-Soler model CHECKED with original

# L. Robledo and A. Soler, "Luminous efficacy of global solar radiation for clear skies,"
# Energy Conversion and Management, vol. 41, pp. 1769-1779, 2000.

# obtained from Reno, M.J., Hansen, C. W., Stein, J.S., 2012. Global Horizontal Irradiance Clear Sky Models:
# Implementation and Analysis. SANDIA REPORT SAND2012-2389.

RS2000 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  G0 <- 1159*(cosThzS^1.179)*exp(-0.0019*(90 - theta))
  D0 <- B0 <- Bn <- zoo(rep(NA, length(G0)), index(G0))
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}

## 33. Meinel 1976 model CHECKED with Meinel
# A. B. Meinel and M. P. Meinel, Applied solar energy. Reading, MA: Addison-Wesley Publishing Co., 1976.

# obtained from Reno, M.J., Hansen, C. W., Stein, J.S., 2012. Global Horizontal Irradiance Clear Sky Models:
# Implementation and Analysis. SANDIA REPORT SAND2012-2389.

M1976 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  Bn <- (Bo0/cosThzS)*0.7^((1/cosThzS)^0.678)
  B0 <- Bn*cosThzS
  G0 <- D0  <- zoo(rep(NA, length(Bn)), index(Bn))

  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 34. Laue 1970 model.
# E. G. Laue, "The measurement of solar spectral irradiance at different terrestrial
# elevations," Solar Energy, vol. 13, pp. 43-50, IN1-IN4, 51-57, 1970.

# obtained from Reno, M.J., Hansen, C. W., Stein, J.S., 2012. Global Horizontal Irradiance Clear Sky Models:
# Implementation and Analysis. SANDIA REPORT SAND2012-2389.

L1970 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  elev <- loc$elev/1000
  
  Bn <- (Bo0/cosThzS)*((1 - 0.14*elev)*0.7^((1/cosThzS)^0.678) + 0.14*elev)
  B0 <- Bn*cosThzS
  G0 <- D0 <- zoo(rep(NA, length(B0)), index(B0))
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}
  

## 35. Kasten 1980. CHECKED with Reno et al.
# F. Kasten, "A Simple Parameterization of the Pyrheliometric Formula for Determining
# the Linke Turbidity Factor," Meteorologische Rundschau, vol. 33, pp. 124-127, 1980.

# Bn has been taken from Ineichen and Perez, 2002

# obtained from Reno, M.J., Hansen, C. W., Stein, J.S., 2012. Global Horizontal Irradiance Clear Sky Models:
# Implementation and Analysis. SANDIA REPORT SAND2012-2389.

K1980 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  Bo0 <- data$Bo0
  
  TL <- data$TL
  m <- data$m
  
  elev <- loc$elev
  
  # correction for TL < 2 from Ineichen and Perez, 2002
  foo <- (TL < 2)
  TL[foo] <- TL[foo] - 0.25*(2 - TL[foo])^0.5
  G0 <- 0.84*Bo0*exp(-0.027*m*(exp(-elev/8000) + exp(-elev/1250)*(TL - 1)))
  Bn <- (0.664 + 0.163/(exp(-elev/8000)))*(Bo0/cosThzS)*exp(-0.09*m*(TL - 1))
  Bn <- pmin(Bn, G0*(1 - (0.1 - 0.2*exp(-TL))/(0.1 + 0.88/exp(-elev/8000)))/cosThzS)
  
  B0 <- Bn*cosThzS
  D0 <- G0 - B0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}

## 36. Ineichen and Perez 2002 CHECKED with original
# P. Ineichen and R. Perez, "A new airmass independent formulation for the Linke turbidity
# coefficient," Solar Energy, vol. 73, pp. 151-157, 2002.

# obtained from Reno, M.J., Hansen, C. W., Stein, J.S., 2012. Global Horizontal Irradiance Clear Sky Models:
# Implementation and Analysis. SANDIA REPORT SAND2012-2389.

IP2002 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  Bo0 <- data$Bo0
  
  TL <- data$TL
  m <- data$m
  p <- data$p
  
  am <- m*p/1013.25
  elev <- loc$elev
  
  
  foo <- (TL < 2)
  TL[foo] <- TL[foo] - 0.25*(2 - TL[foo])^0.5
  TLAM2 <- TL/0.8662
  Bn <- (0.664 + 0.163/(exp(-elev/8000)))*(Bo0/cosThzS)*exp(-0.09*am*(TLAM2 - 1))
  G0 <- (5.09e-5*elev + 0.868)*Bo0*exp(-(3.92e-5*elev + 0.0387)*
           am*(exp(-elev/8000) + exp(-elev/1250)*(TL - 1)))
  Bn <- pmin(Bn, G0*(1 - (0.1 - 0.2*exp(-TL))/(0.1 + 0.88/exp(-elev/8000)))/cosThzS)
  B0 <- Bn*cosThzS
  D0 <- G0 - B0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 37. ASHRAE 2005 Checked with Badescu

# ASHRAE's Handbook of Fundamentals, 2005 (chap. 31).

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

  
A2005  <- function(data, loc){
    
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  month <- month(index(data))
  
  A <- c(1202, 1187, 1164, 1130, 1106, 1092, 1093, 1107, 1136, 1166, 1190, 1204)
  B <- c(0.141, 0.142, 0.149, 0.164, 0.177, 0.185, 0.186, 0.182, 0.165, 0.152, 0.144, 0.141)
  D <- c(0.103, 0.104, 0.109, 0.12, 0.13, 0.137, 0.138, 0.134, 0.121, 0.111, 0.106, 0.103)
  Bn <- A[month]*exp(-B[month]/cosThzS)
  B0 <- Bn*cosThzS
  D0 <- D[month]*Bn
  G0 <- B0 + D0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}

## 38. Bashahu and Laplaze CHECKED with Badescu
# Bashahu, M., Laplaze, D., 1994. An atmospheric model for computing solar radiation.
# Renewable Energy 4(4), 455-458.


# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.


BL1994 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  Bo0 <- data$Bo0
  
  m <- data$m
  pwc <- data$pwc
  p <- data$p
  beta <- data$beta
    
  am <- m*p/1013.25
  amb <- (m*beta)
  T_s <- (5.228/am)*(exp(-0.0002254*am) - exp(-0.1409*am))*
    exp(-(m*beta)/(0.6489^1.3)) + 0.00165*am + 0.2022*(1 - (m*beta)/(2.875^1.3))
  T_a <- (0.1055 + 0.07053*log10(m*pwc + 0.07854))*exp(-(m*beta)/(1.519^1.3))
  Bn <- (Bo0/cosThzS)*(T_s - T_a)
  B0 <- Bn*cosThzS
  D0 <- 0.5*Bo0*(1 - T_s)^1.3333
  G0 <- B0 + D0

  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}

## 39. BCLS model CHECKED with original
# Barbaro, S., Coppolino, S., Leone, C., Sinagra, E., 1979. An atmospheric model for computing direct and 
# diffuse solar radiation. Solar Energy 22, 225-228.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.


BCLS1979 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  m <- data$m
  pwc <- data$pwc
  p <- data$p
  
  
  Bn <- (Bo0/cosThzS)*exp(-0.13491 + -0.00428*pwc - 
         3.68e-5*(- 200))*exp(-m*(0.13708 + 0.00261*pwc + 1.131e-4*(- 200)))
  B0 <- Bn*cosThzS
  foo <- Bo0*(0.938*exp(-0.154*pwc*m)) + 0.1*(41840/60)*(0.004*((10*pwc*m)^2.1) - 1.1086e-5*((10*pwc*m)^3) + 
         121.948*(1 + (10*pwc*m))/(1 + 10*((10*pwc*m)^2))) 
  # I use the formula by Badescu et al.
  D0 <- 0.5*(foo - Bn)*(cosThzS^1.33333) # Badescu et al. uses 1.33 instead of 1/3 of the original paper
  D0[D0<0] <- 0
  G0 <- (D0 + B0)/0.96 ## correction by Badescu et al.
  D0 <- G0 - B0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 40. Biga and Rosa CHECKED with Badescu
# Biga, A.J., Rosa., R., 1979. Contribution to the study of the solar radiation
# climate of Lisbon. Solar Energy 23:1, 61-67. 

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

BR1979 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  
  Bn <- 926*(cosThzS^0.29)
  B0 <- Bn*cosThzS
  D0 <- 131*(cosThzS^0.6)
  G0 <- B0 + D0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}

## 41. Chandra model CHECKED with Badescu

# Chandra, M., 1978. Dependence of solar radiation availability on atmospheric turbidity.
# In: Sun: Mankind's future source of energy; Proceedings of the International Solar Energy Congress, 
# New Delhi, India, January 16-21, 1978. Volume 1. (A79-17276 05-44) 
# Elmsford, N.Y., Pergamon Press, Inc., 1978, p. 430-434.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

C1978 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  TL <- data$TL
  
  # TL limited to 2-5 for convergence
  TL[TL>5] <- 5
  TL[TL<2] <- 2
  
  Bn <- 1.64*cosThzS - 1.10*(cosThzS)^2 + 2.64*(cosThzS)^3
  D0 <- 1.64*cosThzS - 5.13*(cosThzS)^2 + 7.31*(cosThzS)^3 - 3.41*(cosThzS)^4
  
  foo4 <- (TL < 4)
  Bn[foo4] <- 2.58*cosThzS[foo4] - 2.68*(cosThzS[foo4])^2 + 1.10*(cosThzS[foo4])^3
  D0[foo4] <- 0.816*cosThzS[foo4] - 1.30*(cosThzS[foo4])^2 + 1.45*(cosThzS[foo4])^3 - 0.657*(cosThzS[foo4])^4
  
  foo3 <- (TL < 3)
  Bn[foo3] <- 3.27*cosThzS[foo3] - 3.17*(cosThzS[foo3])^2 + 1.03*(cosThzS[foo3])^3
  D0[foo3] <- 0.745*cosThzS[foo3] - 1.56*(cosThzS[foo3])^2 + 1.85*(cosThzS[foo3])^3 - 0.767*(cosThzS[foo3])^4
  
  Bn <- Bn*1353*1.022/1.94
  B0 <- Bn*cosThzS
  D0 <- D0*1353*1.022/1.94
  G0 <- B0 + D0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 42. SH model checked with original and Badescu

# P. W. Suckling and J. E. Hay (1976) "Modelling direct, diffuse and total solar radiation
# for cloudless days." Atmosphere 14: 298-308.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.


SH1977 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  m <- data$m
  p <- data$p
  pwc <- data$pwc
  rhog <- data$rhog
  
  am <- m*p/1013.25
  T_r <- 0.972 - 0.08262*am + 0.00933*am^2 - 0.00095*am^3 + 4.37e-5*am^4
  T_r[T_r>1]<-0
  foo1 <- 1 - 0.0225*(pwc*am)
  foo2 <- 1 - 0.0225*(pwc*1.66)
  foo3 <- 1 - 0.077*((pwc*am)^0.3)
  foo4 <- 1 - 0.077*((pwc*1.66)^0.3)
  Bn <- (Bo0/cosThzS)*T_r*((0.975^am)^2)*foo3*foo1
  B0 <- Bn*cosThzS
  D01 <- 0.6*Bo0*foo3*(0.975^am)*(1 - foo1*(0.975^am)*T_r)
  D02 <- rhog*(Bn*cosThzS + D01)*0.4*foo4*(0.975^1.66)*(1 - foo2*(0.975^1.66)*0.85655)
  D0 <- D01 + D02
  G0 <- B0 + D01 + D02
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 43.  model. DPPLT model

# Goswami, T.K., & Klett, D. E., 1980. Comparison of several models for long term monthly average daily
# insolation on horizontal suraces and the estimation of horizontal
# surface insolation for 16 U.S. locations", ASME paper 80-WA/Sol-28.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

DPPLT <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))

  Bn <- 950*(1 - exp(-0.075*(90 - theta)))
  B0 <- Bn*cosThzS
  D0 <- 2.534 + 3.475*(90 - theta)
  G0 <- B0 + D0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}

## 44. model. Ideriah CHECKED with original

# Ideriah, F.J.K., 1981. A model for calculating direct and D0fuse solar radiation. Solar Energy 26(5), 447-452.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.


I1981 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  m <- data$m
  p <- data$p
  pwc <- data$pwc
  alpha <- data$alpha
  beta <- data$beta
  
  am <- m*p/1013.25
  amb <- m*beta
  T_s <- (5.228/am)*(exp(-0.0002254*am) - exp(-0.140875*am))*exp(-(m*beta)/(0.6489^alpha)) + 
    0.00165*am + 0.2022*(1 - (m*beta)/(2.875^alpha))
  T_abs <- (0.1055 + 0.07053*log10(m*pwc + 0.07854))*exp(-(m*beta)/(1.519^alpha))
  Bn <- (Bo0/cosThzS)*(T_s - T_abs)
  B0 <- Bn*cosThzS
  D0 <- 0.5*T_s*Bo0*(cosThzS^0.33333)
  G0 <- Bn*cosThzS + D0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0

  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 45. Josefsson model CHECKED with Badescu

## Davies, J.A., Mckay, D.C., Luciani, G., Abdel-Wahab, M., 1988. Validation of Models for Estimating Solar Radiation
## on Horizontal Surfaces. Atmospheric Environment Service, Downsview, Ont., IEA Task IX Final Report.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

Josefsson <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  m <- data$m
  p <- data$p
  pwc <- data$pwc
  rhog <- data$rhog
  
  am <- m*p/1013.25
  am66 <- 1.66*p/1013.25
  aw <- 0.29*(10*pwc)*m/((1 + 14.15*(10*pwc)*m)^0.635 + 0.5925*(10*pwc)*m)
  foo <- 0.5248 + 0.007981*(90 - theta)
  kk <- which(((90-theta) > 45) == TRUE)
  if (length(kk)>0){foo[kk] <- 0.856 + 0.000734*(90 - theta[kk])}
  
  T_r <- 0.9768 - 0.0874*am + 0.010607552*am^2 - 0.000846205*am^3 + 3.57246e-5*am^4 - 6.0176e-7*am^5
  T_o3 <- 0.95545
  T_a <- 0.95^m
  T_ab <- 0.95^am66
  T_aa <- 1 - (1 - 0.75)*(1 - T_a)
  T_as <- 1 - 0.75*(1 - T_a)
  foo1 <- 0.93 - 0.21*log(1.66)
  rhoa <- 0.0685 + (1 - T_ab)*0.75*(1 - foo1)
  
  Bn <- (Bo0/cosThzS)*(T_o3*T_aa*T_as*T_r - aw)
  Bn[Bn<0] <- 0
  B0 <- Bn*cosThzS
  D01 <- Bo0*(0.5*T_o3*T_aa*T_as*(1 - T_r) + foo*(T_r*T_o3*T_aa - aw)*(1 - T_as))
  D02 <- rhog*rhoa*(Bn*cosThzS + D01)/(1 - rhog*rhoa)
  D0 <- D01 + D02
  G0 <- B0 + D0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}

## 46. Modified Kasten model CHECKED with Badescu, 1997
# Badescu, V., 1997. Verification of some very simple clear and cloudy sky models to evaluate global solar irradiance
# Solar Energy 61(4), 251-264.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

mKasten <- function(data, loc){
  
  cosThzS <- data$cosThzS
  Bo0 <- data$Bo0
  
  m <- data$m
  p <- data$p
  pwc <- data$pwc
  
  am<- m*p/1013.25
  T_r <- 0.9768 - 0.0874*am + 0.010607552*am^2 - 0.000846205*am^3 + 3.57246e-5*am^4 - 6.0176e-7*am^5
  T_g <- 1 - (0.007413*m/(1 + 0.0147*m + 3.95675e-05*m^2)) - 
    (0.3787*m/((1 + 48.51*m)^0.805) + 0.02303*m/(1 + (36.26*m)^3))
  a_w <- 0.29*(10*pwc*am)/((1 + 14.15*(10*pwc*am))^0.635 + 0.5925*(10*pwc*am))
  T_a <- 0.9^am
  Bn <- (Bo0/cosThzS)*(T_r*T_g - a_w)*T_a
  Bn[Bn<0] <- 0
  B0 <- Bn*cosThzS
  D0 <- Bo0*(T_g*(1 - T_r)*0.5 + (T_g*T_r - a_w)*(1 - T_a)*0.75*(0.93 - 0.21*log(am)))
  G0 <- B0 + D0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 47. King model CHECK DNI with original

# King, R., Buckies, R.O., 1979. Direct solar tansmittance for a clear sky. Solar Energy 22, 297-301.

# Buckius R. O, King R. Diffuse solar radiation on a horizontal surface for a clear
# sky. Solar Energy 1978;21:503–9.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

K1979 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  Bo0 <- data$Bo0
  
  m <- data$m
  p <- data$p
  pwc <- data$pwc
  alpha <- data$alpha
  beta <- data$beta
  rhog <- data$rhog
  
  am <- m*p/1013.25
  amb <- m*beta
  a1 <- 0.5
  T_D <- (5.228/am)*(exp(-0.0002254*am) - exp(-0.1409*am))*exp(-(m*beta)/(0.6489^alpha)) + 
    0.00165*am + 0.2022*(1 - (m*beta)/(2.875^alpha))
  T_abs <- (0.1055 + 0.07053*log10(m*pwc + 0.07854))*exp(-(m*beta)/(1.519^alpha))
  Bn <- (Bo0/cosThzS)*(T_D - T_abs)
  B0 <- Bn*cosThzS
  T_l <- 0.8336 + 0.17257*beta - 0.64595*exp(-7.4296*(beta^1.5))
  T_ll <- exp(-am*T_l)
  foo1 <- 2*(1 + 2*cosThzS + (1 - 2*cosThzS)*T_ll)
  foo2 <- (rhog - 1)*(0.5*T_l - 4*T_l - 4) + 4*rhog
  D0 <- 0.634*(foo1/foo2 - T_ll)*Bo0
  G0 <- D0 + B0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 48. Machler Iqbal checked with Badescu

# Mächler, M.A., Iqbal, M., 1985. A modification of the ASHRAE clear sky irradiation model.
# ASHRAE Transactions 91, Pt 1.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

# vis visibility in km if NA is calculated

MI1985 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  alpha <- data$alpha
  beta <- data$beta
  vis <- data$vis
  p <- data$p
  pwc <- data$pwc
  
  ams <- (p/1013.25)/cosThzS
  # if visibility is NA is obtained with King & Buckies
  # King, R., Buckies, R.O., 1979. Direct solar tansmittance for a clear sky. Solar Energy 22, 297-301.
  a <- (is.na(vis) == TRUE)
  foo <- beta/(0.55^alpha)
  vis[a] <- (0.084987 - foo[a] + (foo[a]^2 - 0.169974*foo[a] + 0.0117554)^0.5)/0.000574492
  T_a <- 0
  a <- ((1 - 1.13/(vis^0.57)) > 0)
  T_a[a] <- (1 - 1.13/(vis[a]^0.57))^(ams[a]^0.85)
  T_w <- (1.0223 - 0.00149*pwc)^(ams^0.27)
  Bn <- (Bo0/cosThzS)*T_a*T_w*(0.775^(ams^0.5))
  B0 <- Bn*cosThzS
  D0 <- Bn*(0.1 + 3/vis)
  G0 <- D0 + B0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 49. METSTAT CHECKED with Badescu

# Maxwell, E. L. 1998. METSTAT--the solar radiation model used in the production 
# of the national solar radiation data base (nsrdb). Solar Energy, 62, 263-279.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

# taua - Unsworth-Monteith broadband aerosol optical depth

METSTAT <- function(data, loc){
  
  cosThzS <- data$cosThzS
  Bo0 <- data$Bo0
  
  m <- data$m
  p <- data$p
  pwc <- data$pwc
  uo <- data$uo
  rhog <- data$rhog
  taua <- data$taua/100
  
  am <- m*p/1013.25
  T_r <- exp(-0.0903*(am^0.84)*(1 + am - am^1.01))
  T_r[T_r>1] <- 0
  T_o3 <- 1 - 0.1611*m*uo/((1 + 139.48*m*uo)^0.3035) - 0.002715*m*uo/(1 + 0.044*m*uo + 0.0003*(m*uo)^2)
  T_g <- exp(-0.0127*am^0.26)
  T_w <- 1 - 1.668*m*pwc/((1 + 54.6*m*pwc)^0.637 + 4.042*m*pwc)
  T_a <- exp(-m*taua)
  Bn <- (Bo0/cosThzS)*0.9751*T_r*T_g*T_o3*T_w*T_a
  B0 <- Bn*cosThzS
  T_aa <- 1 - (1 - 0.9)*(1 - m + m^1.06)*(1 - T_a)
  T_as <- T_a/T_aa
  T_s <- (T_o3*T_g*T_aa*(1.34 - 0.5*T_r - 0.84*T_a))*(0.38 + 0.925*exp(-0.851*m))
  T_sg <- (0.9751*T_r*T_g*T_o3*T_w*T_a + T_s)*(rhog - 0.2)*(0.0685 + (1 - 0.84)*(1 - T_as))
  D0 <- (T_s + T_sg)*Bo0
  G0 <- (0.9751*T_r*T_g*T_o3*T_w*T_a + T_s + T_sg)*Bo0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 50. MRM4 model CHECKED with original

# Muneer, T., Gul, M., Kambezedis, H., 1998. Evaluation of an all-sky meteorological radiation
# model against long-term measured hourly data. Energy Conversion and Management 39(3–4), 303–317.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

MRM4 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  Bo0 <- data$Bo0
  
  m <- data$m
  p <- data$p
  pwc <- data$pwc
  uo <- data$uo
  
  k <- 0.387
  am <- m*(p/1013.25)
  T_r <- 0.972 - 0.08262*am + 0.00933*am^2 - 0.00095*am^3 + 0.000437*am^4
  T_w <- 1 - 2.4959*m*pwc/((1 + 79.034*m*pwc)^0.6828 + 6.385*m*pwc) 
  T_o3 <- 1 - (0.1611*uo*m/((1 + 139.48*uo*m)^0.3035) +
                 -0.002715*uo*m/(1 + 0.044*uo*m + 0.0003*(uo*m)^2))
  T_g <- exp(-0.0127*am^0.26)                                            
  T_a <- exp(-k^0.873*(1 + k - k^0.7088)*am^0.9108)
  T_aa <- 1 - 0.1*(1 - T_a)*(1 - m + m^1.06)
  T_as <- 10^(-0.045*am^0.7)
  Bn <- (Bo0/cosThzS)*T_a*T_r*T_o3*T_w*T_g
  B0 <- Bn*cosThzS
  G0 <- (Bo0*T_o3*T_w*T_g*T_aa/(1 - m + m^1.02))*(0.84*(1 - T_as) + 0.5*(1 - T_r))
  D0 <- G0 - B0
  a <- ((T_a*T_r*T_o3*T_w*T_g) >= 0.015)
  
  D0[a] <- (0.285/((T_a[a]*T_r[a]*T_o3[a]*T_w[a]*T_g[a])^1.006))*Bo0[a]*T_a[a]*T_r[a]*T_o3[a]*T_w[a]*T_g[a]
  G0[a] <- Bo0[a]*T_a[a]*T_r[a]*T_o3[a]*T_w[a]*T_g[a]*(1 + (0.285/((T_a[a]*T_r[a]*T_o3[a]*T_w[a]*T_g[a])^1.006)))
  B0[a] <- G0[a] - D0[a]
  Bn[a] <- B0[a]/cosThzS[a]
  
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}
  

## 51. MRM5 CHECKED with Badescu et al
# Badescu, V., 2008. Modeling Solar Radiation at the Earth's Surface. Springer Book.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

# atype 'rural' or 'urban' type of aerosol. default 'rural' 


MRM5 <- function(data, loc = 'rural'){
  
  cosThzS <- data$cosThzS
  Bo0 <- data$Bo0
  
  m <- data$m
  p <- data$p
  rhog <- data$rhog
  pwc <- data$pwc
  uo <- data$uo
  beta <- data$beta
  
  am <- m*p/1013.25
  m102 <- m^1.02
  rhog[(rhog<=0)] <- 0.2
  T_r <- exp(-(0.1128*(am^0.8346)*(0.9341 - (am^0.9868) + (0.9391*am))))
  T_w <- 1 - (3.014*m*pwc)/(((1 + (119.3*m*pwc))^0.644) + (5.814*m*pwc))
  T_o3 <- 1 - (0.2554*m*uo)/(((1 + (6107.26*m*uo))^0.204) + 0.471*m*uo)
  T_co2 <- 1 - (0.0721*am*350)/(((1 + (377.89*am*350))^0.5855) + 3.1709*am*350)
  T_co <- 1 - (0.0062*am*0.075)/(((1 + (243.67*am*0.075))^0.4246) + 1.7222*am*0.075)
  T_n2o <- 1 - (0.0326*am*0.28)/(((1 + (107.413*am*0.28))^0.5501) + 0.9093*am*0.28)
  T_ch4 <- 1 - (0.0192*am*1.6)/(((1 + (166.095*am*1.6))^0.4221)+ 0.7186*am*1.6)
  T_o2 <- 1 - (0.0003*am*2.095e05)/(((1 + (476.934*am*2.095e05))^0.4892) + 0.1261*am*2.095e05)
  
  # Aerosol scattering Yang, K., Huang, G.W., Tamai, N., 2001. A hybrid model for estimating global solar radiation.
  ## Solar Enery 70, 13-22.
  bm <- beta*am
  a3 <- 0.677 + 0.1464*bm - 0.00626*bm^2
  a3[(a3<=0)] <- 0.00001
  T_aer <- exp(-bm/(a3^1.3))

  #omega <- ifelse(loc == 'urban', 0.6, 0.9)
  omega <- 0.9
  T_ab <- 1 - ((1 - omega)*(1 - am + (am^1.06))*(1 - T_aer))
  bm166 <- beta*1.66
  a166 <- 0.677 + 0.1464*bm166 - 0.00626*bm166^2
  rsa <- rhog*(0.0685 + 0.16*(1 - exp(-bm166/(a166^1.3))))
  
  Bn <- Bo0*T_aer*T_r*T_w*T_o3*T_co2*T_co*T_n2o*T_ch4*T_o2/cosThzS
  B0 <- Bn*cosThzS
  D01 <- Bo0*T_w*T_o3*T_co2*T_co*T_n2o*T_ch4*T_o2*T_ab*(0.5*(1 - (T_aer/T_ab)*T_r))
  D02 <- (Bo0*T_aer*T_r*T_w*T_o3*T_co2*T_co*T_n2o*T_ch4*T_o2 + D01)*rsa/(1 - rsa)
  D0 <-  D01 + D02
  G0 <- Bo0*T_aer*T_r*T_w*T_o3*T_co2*T_co*T_n2o*T_ch4*T_o2 + D0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}
  

## 52. NAD model CHECKED with both

# Nijegorodov, N., Adedoyin, J.A., Devan, K.R.S., 1997. A new analytical-empirical model for the instantaneous
# Diffuse-radiation and experimental investigation of its validity. Renewable Energy 11, 341-350.

# with tranmsittances from Bird's model as obtained from Badescu et al.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

# Temp temperature in K

NAD1997 <- function(data, loc){

  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  m <- data$m
  p <- data$p
  pwc <- data$pwc
  uo <- data$uo
  alpha <- data$alpha
  beta <- data$beta
  TempK <- data$Temp + 273
  
  elev <- loc$elev
  
  RsH <- 6367/29.7
  foo <- (((1 + 2*RsH + (RsH^2)*(cosThzS^2))^0.5) - RsH*cosThzS)*exp(-0.034169*elev/TempK)
  am <- m*p/1013.25

  T_r <- exp(-0.0903*(am^0.84)*(1 + am - am^1.01))
  T_o3 <- 1 - 0.1611*(m*uo)/((1 + 139.48*(m*uo))^0.3035) - 0.002715*(m*uo)/(1 + 0.044*(m*uo) + 0.0003*(m*uo)^2)
  T_u <- exp(-0.0127*am^0.26)
  T_w <- 1 - 2.4959*(m*pwc)/((1 + 79.034*(m*pwc))^0.6828 + 6.385*(m*pwc))
  taua <- beta*(0.2758*(1/(0.38^alpha)) + 0.35*(1/(0.5^alpha)))
  T_a <- exp(-(taua^0.873)*(1 + taua - taua^0.7088)*(m*pwc)^0.9108)
  Bn <- (Bo0/cosThzS)*T_r*T_o3*T_u*T_w*T_a
  B0 <- Bn*cosThzS
  foo1 <- -log(T_r*T_o3*T_u*T_w*T_a)
  foo2 <- cosThzS/((2.634/(1 + 1.4196*foo)) - 0.0885)
  foo3 <- (0.5*(-log(T_r))/(0.7*(-log(T_a))))*(foo2 -1) + foo2
  D0 <- Bo0*(foo1*exp(-(foo - 1)*foo1))*(0.5*(-log(T_r)) + 0.7*(-log(T_a))*foo3)/(exp(foo1) - 1)
  D0[D0<0] <- 0
  G0 <- D0 + B0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}
  
## 53. BDG model CHECKED with original and complete model with Badescu et al.

# Belcher, B.N. and A. T. DeGaetano, 2007.
#A revised model to estimate solar radiation using automated surface 
# weather stations, Solar Energy, 81(3), 329-345.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

BDG2007 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  month <- month(index(data))
  
  p <- data$p
  pwc <- data$pwc
  vis <- data$vis
  RH <- data$RH
  
  k0 <- c(0.9603, 0.9682, 0.9402, 0.822)
  k1 <- c(-0.0024, -0.0028, -0.0019, 0)
  k2 <- c(0.0003, 0.0005, 0.00005, 0.0014)
  k3 <- c(0.000025, 0.000032, 0.00002, 0)
  k4 <- c(0.815, 0.9827, 1.0473, 0.67)
  k5 <- c(-0.0028, -0.006, -0.0058, 0)
  k6 <- c(0.0005, -0.0009, -0.0025, 0.0024)
  k7 <- c(0.000034, 0.000054, 0.000074, 0)
  
  m <- 35/(1224*cosThzS + 1)^0.5
  amp <- m*p/1013.25
  
  T_g <-  1.021 - 0.084*((m*(949*((p/1013.25)*101.325)*1e-5 + 0.051))^0.5)
  T_w <- 1 - 0.077*(m*pwc)^0.3
  
  foo <- rep(4, length(m))
  foo[((month >= 3) & (month <= 5))] <- 1
  foo[((month >= 6) & (month <= 8))] <- 2
  foo[((month >= 9) & (month <= 11))] <- 3
  
  T_a <- ((k0[foo] + k1[foo]*RH) + (k2[foo] + k3[foo]*RH)*theta)^m
  T_ap <- ((k4[foo] + k5[foo]*RH) + (k6[foo] + k7[foo]*RH)*theta)^m
  T_c <- 0.9999
  vism <- vis/1.6093
  T_s <- rep(1, length(vism))
  if (length(which(vism<10))>0){T_s[vism<10] <- 0.98}
  if (length(which(vism<1))>0){T_s[vism<1] <- 0.74}
  T_s[(vism<10)] <- 0.98
  T_s[(vism<1)] <- 0.74
  Bn <- (Bo0/cosThzS)*T_g*T_w*T_ap*T_c*T_s
  B0 <- Bn*cosThzS
  G0 <- Bo0*T_g*T_w*T_a*T_c*T_s
  D0 <- G0 - B0
  D0[D0<0] <- 0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}

## 54. Paltridge & Platt model CHECKED with Badescu

# Paltridge, G.W., Platt, C.M.R., 1976. Radiative processes in meteorology and climatology.
# Amsterdam-Oxford-New York, Elsevier Scientific Publishing Company.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

PP1976 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  
  Bn <- 1000*(1 - exp(-0.06*(90 - theta)))
  B0 <- Bn*cosThzS
  D0 <- 5 + 96*(1 - exp(-0.05*(90 - theta)))
  G0 <- D0 + B0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 55. modified Perrin CHECKED with Badescu

## Perrin de Brichambaut, C. Vauge, C., 1982. Le gisement solaire: Evaluation de la ressource energetique (Paris: 
## Technique et documentation).

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

mP1982 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  m <- data$m
  p <- data$p
  uo <- data$uo
  pwc <- data$pwc
  beta <- data$beta
  
  am <- m*p/1013.25
  akw <- m*pwc
  xw <- xwg <- 0
  xw[(akw>0)] <- log(m*pwc)
  xwg[((am*pwc) > 0)] <- log(am*pwc)

  Bn <- (Bo0/cosThzS)*exp(-0.031411 - 0.064331*am)*
    exp(-1.4327*m*beta)*(1 - (0.015 + 0.024*(m*uo)) - 
                           (0.1 + 0.03*xw + 0.002*xw^2) - (0.013 - 0.0015*xwg))
  B0 <- Bn*cosThzS
  foo <- -(log(Bn/(Bo0/cosThzS)))/((1/(9.4 + 0.9*am))*am)
  G0 <- (1270 - 56*foo)*cosThzS^((foo + 36)/33)
  D0 <- G0 - B0
  D0[(D0<0)] <- 0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 56. modified Psiloglou CHECKED with Badescu

# Psiloglou, B.E., Santamouris, M., Asimakopoulos, D.N., 2000. Atmospheric broadband model
# for computation of solar radiation at the earth's surface. Application to Mediterranean climate.  
# Pure and Applied Geophysics 157(5), 829-860.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

mP2000 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  m <- data$m
  p <- data$p
  uo <- data$uo
  pwc <- data$pwc
  beta <- data$beta
  rhog <- data$rhog
  
  am <- m*p/1013.25
  T_r <- exp(-0.1128*(am^0.8346)*(0.9341 - am^0.9868 + 0.9391*am))
  ako3 <- m*uo
  T_o3 <- 1 - 0.2554*ako3/((1 + 6107.26*ako3)^0.204 + 0.471*ako3)
  akco2 <- am*330
  T_co2 <-  1 - 0.0721*akco2/((1 + 377.89*akco2)^0.5855 + 3.1709*akco2)
  akco <- am*0.075
  T_co <- 1 - 0.0062*akco/((1 + 243.67*akco)^0.4246 + 1.7222*akco)
  akch4 <- am*1.6
  T_ch4 <- 1 - 0.0192*akch4/((1 + 166.095*akch4)^0.4221 + 0.7186*akch4)
  akn2o <- am*0.28
  T_n2o <- 1 - 0.0326*akn2o/((1 + 107.413*akn2o)^0.5501 + 0.9093*akn2o)
  ako2 <- am*2.095e+05
  T_o2 <- 1 - 0.0003*ako2/((1 + 476.934*ako2)^0.4892 + 0.1261*ako2)
  akw <- (m*pwc)
  T_w <- 1 - 3.014*akw/((1 + 119.3*akw)^0.6440 + 5.814*akw)
  # Aerosol transmittance from REST
  b1 <- (-0.013029 + 0.13126*beta)/(1 + 0.42003*beta)
  b2 <- (-0.0083581 + 0.40323*beta + 0.123*beta^2)/(1 + 0.42003*beta)
  T_a <- exp(-m*(beta*(1.6933 + m*b1)/(1 + b2*m)))
  T_aa <- 1.0 - 1.405e-3*m - 9.013e-5*m^2 + 2.2e-6*m^3
  T_as <- T_a/T_aa
  Bn <- (Bo0/cosThzS)*T_r*T_co2*T_co*T_ch4*T_n2o*T_o2*T_o3*T_w*T_a
  B0 <- Bn*cosThzS
  D01 <- Bo0*T_co2*T_co*T_ch4*T_n2o*T_o2*T_o3*T_w*T_aa*(1 - T_as*T_r)/2
  D02 <- (D01 + B0)*rhog*0.08503/(1 - rhog*0.08503)
  D0 <- D01 + D02
  G0 <- B0 + D0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 57. PSI CHECKED with Badescu

# Gueymard, C.A. 1993. Mathematically integrable parameterization of clear sky beam and global 
# irradiances and its use in daily irradiation applications. Solar Energy 50, 385–397.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

PSI <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  p <- data$p
  pwc <- data$pwc
  beta <- data$beta
  rhog <- data$rhog
  
  foo1 <- 1 + (0.1594 - 0.226*beta)*(1 - p/1013.25)
  foo2 <- 1 + (0.0752 - 0.107*beta)*(1 - p/1013.25)
  foo3 <- (1 - 0.2*(0.061 + 0.072*beta^0.5))/(1 - rhog*(0.061 + 0.072*beta^0.5))
  inc <- cosThzS
  inc[inc<0.0285] <- 0.0285
  b <- log(1 + 10*beta)
  
  k0 <- 0
  k1 <- 4.41547 - 5.09266*b + 1.471865*b^2
  k2 <- -18.45187 + 38.3584*b - 22.7449*b^2 + 4.3189*b^3
  k3 <- 31.2506 - 74.53835*b + 48.355*b^2 - 9.8657*b^3
  k4 <- -25.18762 + 64.3575*b - 43.6586*b^2 + 9.23151*b^3
  k5 <- 7.64179 - 20.41687*b + 14.20502*b^2 - 3.06053*b^3
  
  a <- (beta < 0.175)
  k0[a] <- exp(-1.62364 - 6.8727*b[a])
  k1[a] <- 2.94298 + 5.23504*b[a] - 18.23861*(b[a])^2 + 11.1652*(b[a])^3
  k2[a] <- -8.12158 - 15.8*b[a] + 69.2345*(b[a])^2 - 45.1637*(b[a])^3
  k3[a] <- 12.5571 + 25.44*b[a] - 123.3933*(b[a])^2 + 83.1014*(b[a])^3
  k4[a] <- -9.8044 - 20.3172*b[a] + 103.9906*(b[a])^2 - 71.3091*(b[a])^3
  k5[a] <- 3.00487 + 6.3176*b[a] - 33.3891*(b[a])^2 + 23.1547*(b[a])^3
  
  g0 <- 0.006
  g1 <- 0.387021 - 0.386246*b + 0.09234*b^2
  g2 <- 1.35369 + 1.533*b - 1.07736*b^2 + 0.23728*b^3
  g3 <- -1.59816 - 1.903774*b + 1.631134*b^2 - 0.3877*b^3
  g4 <- 0.66864 + 0.80172*b - 0.75795*b^2 + 0.18895*b^3
  
  Bn <- (Bo0/cosThzS)*(exp(0.155*(1 - pwc^0.24)))*foo1*(k0 + k1*inc + k2*inc^2 + k3*inc^3 + k4*inc^4 + k5*inc^5)
  B0 <- Bn*cosThzS
  G0 <- (Bo0/cosThzS)*(exp(0.155*(1 - pwc^0.24)))*foo2*(g0 + g1*inc + g2*inc^2 + g3*inc^3 + g4*inc^4)*foo3
  D0 <- G0 - B0
  D0[D0<0] <- 0

  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 58. Rodgers, Souster, Page model

# Rodgers, G.G., Souster, C.G., Page, J.K., 1978. The development of an interactive computer program SUN1 for 
# the calculation of solar irradiances and daily irradiations on horizontal surfaces on cloudless days for
# given conditions of sky clarity and atmospheric water content. In Report BS 28, Department of Building 
# Science, University of Sheffield, UK.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

# obtained Gueymard, C., 2003. Direct solar transmittance and irradiance predictions with
# broadband models. Part I: detailed theoretical performance assessment. Solar Energy 74, 355–379.

# taua - Unsworth-Monteith broadband turbidity coefficient

RSP1978 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  p <- data$p
  pwc <- data$pwc
  taua <- data$taua/100
  
  foo1 <- foo2 <- list(rep(0, length(p)))
  k0 <- c(-0.0642111, -0.0046883, 0.000844097)
  k1 <- c(-0.00801046, 0.00220414, -0.000191442)
  k2 <- c(0.000153069, -0.0000429818, 0.00000374176)
  k3 <- c(0, 47.382 , 29.671, -15.8621, 4.3463, -0.57764, 0.03472, -0.0007362)
  k4 <- c(297, 1.8313, -3.7082, 4.1233, -0.6409, 0.02855)
  
  m <- 1/cosThzS
  m[(theta>80)] <- exp(3.67985 - 24.4465*cosThzS[(theta>80)] + 154.017*(cosThzS[(theta>80)])^2 - 
                    742.181*(cosThzS[(theta>80)])^3 + 2263.36*(cosThzS[(theta>80)])^4 - 
                    3804.89*(cosThzS[(theta>80)])^5 + 2661.05*(cosThzS[(theta>80)])^6)
  
  am <- m*p/1013.25
  foo <- -0.129641 + 0.0412828*pwc - 0.0112096*pwc^2
  
  foo0 <- sapply(c(1:3), function(i){k0[i] + k1[i]*(10*pwc) + k2[i]*(10*pwc)^2})
  for(i in 2:7){foo1[[i]] <- foo1[[(i-1)]] + k4[(i-1)]*(((90 - theta)/10))^(i-2)}
  foo11 <- 1e-3*foo1[[7]]
  for(i in 2:8){foo2[[i]] <- foo2[[(i-1)]] + k3[i]*(((90 - theta)/10))^(i-1)}
  foo22 <- 2 + foo2[[8]]
  
  foo3 <- exp(foo + foo0[,1]*am + foo0[,2]*am^2 + foo0[,3]*am^3)
  Bn <- (Bo0/cosThzS)*foo3*exp(-am*taua)
  Bn[Bn>(Bo0/cosThzS)] <- Bo0/cosThzS
  B0 <- Bn*cosThzS
  D0 <- (Bo0/cosThzS/1347)*foo22 - foo11*Bn*cosThzS
  D0[(D0<0)] <- 0
  G0 <- B0 + D0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 59. Carroll model 1985 CHECKED with original
# Carroll, J.J., 1985. Global transmissivity and D0fuse fraction of solar radiation
# for clear and cloudy skies as measured and as predicted by bulk transmissivity models.
# Solar Energy 35(2), 105-118. 

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

C1985 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  m <- data$m
  p <- data$p
  pwc <- data$pwc
  rhog <- data$rhog
  beta <- data$beta
  
  am <- m*p/1013.25
  T_r <- 1/10^(0.054 - 0.0088*am + 0.00108*am^2 - 0.000051*am^3)
  T_w <- 1/(10^((0.04*pwc^0.1 + 0.01*pwc)*am))
  T_g <- 1/(10^(0.02*am))
  T_a <- 1/(10^(0.666*am*beta))
  Bn <- (Bo0/cosThzS)*T_g*T_w*T_r*T_a
  B0 <- Bn*cosThzS
  D01 <- (0.5 + 0.3*beta)*(cosThzS^0.33333)*(Bo0*T_g*T_w - B0)
  D02 <- (0.1*rhog + 0.36*beta*(rhog - 0.25) - 0.02)*(B0 + D01)
  D0 <- D01 + D02
  G0 <- B0 + D0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 60. Santamouris model

# Psiloglou, B.E., Santamouris, M., Asimakopoulos, D.N., 2000. Atmospheric broadband model
# to computation of solar radiation at the Earth's surface.  Application to Mediterranean 
# climates. Journal of Pure and Applied Geophysics, 157(5), 829-860.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.


PSA2000 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  m <- data$m
  p <- data$p
  pwc <- data$pwc
  uo <- data$uo
  rhog <- data$rhog
  
  am <- m*p/1013.25
  T_r <- exp(-0.1128*(am^0.8346)*(0.9341 - am^0.9868 + 0.9391*am))
  ako3 <- m*uo
  T_o3 <- 1 - 0.2554*ako3/((1 + 6107.26*ako3)^0.204 + 0.471*ako3)
  akco2 <- am*330
  T_co2 <- 1 - 0.0721*akco2/((1 + 377.89*akco2)^0.5855 + 3.1709*akco2)
  akco <- am*0.075
  T_co <- 1 - 0.0062*akco/((1 + 243.67*akco)^0.4246 + 1.7222*akco)
  akch4 <- am*1.6
  T_ch4 <- 1 - 0.0192*akch4/((1 + 166.095*akch4)^0.4221 + 0.7186*akch4)
  akn2o <- am*0.28
  T_n2o <- 1 - 0.0326*akn2o/((1 + 107.413*akn2o)^0.5501 + 0.9093*akn2o)
  ako2 <- am*2.095e+05
  T_o2 <- 1 - 0.0003*ako2/((1 + 476.934*ako2)^0.4892 + 0.1261*ako2)
  Trmg <- T_co2*T_co*T_ch4*T_n2o*T_o2
  T_w <- 1 - 3.014*m*pwc/((1 + 119.3*m*pwc)^0.6440 + 5.814*m*pwc)
  T_a <- 1 - 0.2579*m/((1/(1 + 0.04001*m)^2.8451) + 0.2748*m)
  T_aa <- 1 - 1.405e-3*m - 9.013e-5*m^2 + 2.2e-6*m^3
  Bn <- (Bo0/cosThzS)*T_r*T_co2*T_co*T_ch4*T_n2o*T_o2*T_o3*T_w*T_a
  B0 <- Bn*cosThzS
  D01 <- Bo0*T_co2*T_co*T_ch4*T_n2o*T_o2*T_o3*T_w*T_aa*(1 - (T_a/T_aa)*T_r)/2
  D02 <- (D01 + B0)*rhog*0.08503/(1 - rhog*0.08503)
  D0 <- D01 + D02
  G0 <- B0 + D0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 61. Schulze 1976 CHECKED with original

# Schulze, R.E., 1976. A physically based method of estimating solar radiation from suncards.
# Agricultural Meteorology 16(1), 85-101.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

S1976 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  Bn <- 1127*(0.888^(1/cosThzS)) # Badescu et al. uses 0.808 instead original 0.888
  B0 <- Bn*cosThzS
  D0 <- 94.23*(cosThzS^0.5)
  G0 <- B0 + D0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}
  
  
## 62. Sharma & Pal 1965 CHECKED with original

# Sharma, M.R., Pal, R.S., 1965. Total, Direct, and Diffuse Solar Radiation in the Tropics. 
# Solar Energy 9(4), 183-192.


SP1965 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  Bn <- 1.842*((Bo0/cosThzS)/2)*cosThzS/(0.3135 + cosThzS)
  B0 <- Bn*cosThzS
  G0 <- 4.5*((Bo0/cosThzS)/120) + 1.071*Bn*cosThzS
  D0 <- G0 - B0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}



## 63. Watt 1978

# Watt, A.D., 1978. On the nature and distribution of solar radiation. Report HCP/T2552-01, US DOE.

# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

W1978 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  m <- data$m
  p <- data$p
  pwc <- data$pwc
  uo <- data$uo
  rhog <- data$rhog
  alpha <- data$alpha
  beta <- data$beta
  
  K1 <- c(0, 0, 15, 20)
  K2 <- c(30, 1.5, 25, 40)
  am <- m*p/1013.25
  B <- beta*(2^alpha)*0.434294 - 0.02
  T_t <- (0.6*((beta*(2^alpha)*0.434294 - 0.02) - 0.01*pwc - 0.03) + 0.02)
  T_t[T_t<0] <- 0
  
  amm <- sapply(c(1:4), function(i){
    foo <- (((6.4e03/K2[i])*cosThzS)^2 + 2*(6.4e03/K2[i]) + 1)^0.5 - (6.4e03/K2[i])*cosThzS
    if(K1[i] > 0){
      foo1 <- (((6.4e03/K1[i])*cosThzS)^2 + 2*(6.4e03/K1[i]) + 1)^0.5 - (6.4e03/K1[i])*cosThzS
      foo <- (K2[i]*foo - K1[i]*foo1)/(K2[i] - K1[i])}
    foo})
  
  foo1 <- 0.045*(p/1013.25*amm[,1])^0.7
  foo2 <- 0.0071 + 0.001*10*uo*amm[,4]
  foo3 <- 0.0095*pwc*amm[,2]
  foo4 <- 0.6*(B - 0.01*pwc - 0.03)*(amm[,2]^0.7)
  foo5 <- 0.02*amm[,3]
  
  T_w <- 0.93 - 0.033*log10(pwc*amm[,2])
  T_n <- T_w*10^(- foo1 - foo2 - foo3 - foo4 - foo5)
  Bn <- (Bo0/cosThzS)*T_n
  B0 <- Bn*cosThzS
  acs <- (0.93 - 0.033*log10(pwc))*10^(-0.006*p/1013.25 - 0.4*T_t)
  as <- acs*(1 - 10^(-0.003*p/1013.25 - 0.01*pwc - 0.4*T_t))
  D0 <- (Bo0/cosThzS)*(0.8*as*(1 + rhog*as)*(1 + cosThzS)^0.5) + Bo0*0.5*as*(1 + acs*rhog)
  G0 <- D0 + B0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 64. Wesely & Lipschutz diffuse CHECKED with original 

# Wesely, M.L., Lipschutz, R.C., 1976. An experimental study of the effects of aerosols on Difuse and direct
# solar radiation received during the summer near Chicago. Atmospheric Environment 10(11), 981–987.

# Wesely, M.L., Lipschutz, R. C. 1976. A method for estimating hourly averages of diffuse and direct 
# solar radiation under a layer of scattered clouds. Solar Energy 18(5), 467–473.


# obtained from Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

# if vis NA then is calculated from angstrom alfa and beta: 
# King, R., Buckies, R.O., 1979. Direct solar tansmittance for a clear sky. Solar Energy 22, 297-301.

WL1967 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  m <- data$m
  p <- data$p
  vis <- data$vis
  alpha <- data$alpha
  beta <- data$beta
  
  am <- m*p/1013.25
  foo <- beta/0.55^alpha
  vis[is.na(vis)] <- (0.084987 - foo + (foo^2 - 0.169974*foo + 0.0117554)^0.5)/0.000574492
  
  Bn <- exp(-1.6/vis)*(1050/am)*exp(-0.1*am)
  B0 <- Bn*cosThzS
  D0 <- 0.9*Bn*(1 - exp(-0.1*am))
  G0 <- B0 + D0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 65. Zhang CHECKED with original

# Zhang, Q., 2006. Development of the typical meteorological database for Chinese locations.
# Energy and Buildings 38(11), 1320-1326.

# deltaT dry-bulb temperature T at time t, T(t), minus T(t-3 hours)

Z2006 <- function(data, loc){

  Bo0 <- data$Bo0
  cosThzS <- data$cosThzS

  #Ta <- data$Ta
  deltaT <- data$deltaT
  #deltaT <- diff(Ta, lag = 3, na.pad = TRUE)
  RH <- data$RH
  ## Global irradiance
  G0 <- (Bo0*(0.64223 + 0.025258*deltaT - 0.002675*RH) - 38.7179)/0.89794
  G0[G0 < 0] <- 0
  G0[G0 > 1360] <- 1360
  G0 <- zoo(G0, index(data))
  ## Direct Irradiance
  kt <- G0/Bo0
  foo1 <- -0.1556*cosThzS^2 + 0.1028*cosThzS + 1.3748
  foo2 <- 0.7973*cosThzS^2 + 0.1509*cosThzS + 3.035
  foo3 <- 5.4307*cosThzS + 7.2182
  foo4 <- 2.99 # corrected typo with Badescu et al., 2013
  xkn <- foo1 * foo2^(-foo3 * foo2^(-foo4*kt))
  xkn[xkn>1] <- 1
  B0 <- Bo0 * xkn
  Bn <- B0/cosThzS
  ## Diffuse Irradiance
  D0 <- G0 - B0
  D0[D0 < 0] <- 0
  z <- cbind(G0, D0, B0, Bn)
  coredata(z)[is.na(z)] <- 0
  z
}


## 66. Hourwitz and Schulze CHECKED with Badescu

# direct: Hourwitz, B., 1945. Insolation in relation to cloudiness and cloud density. 
# Journal of Meteorology 2, 154-156.
# Hourwitz, B., 1946. Insolation in relation to cloud type. Journal of Meteorology 3, 123-124.

# Diffuse: Schulze, R.E., 1976. A physically based method of estimating solar radiation from suncards.
# Agricultural Meteorology 16(1), 85-101.

# obtained from: Badescu, V., 1998. Verification of some very simple clear and cloudy sky model to evaluate
# global solar irradiance. Solar Energy 61(4), 251-264.
# also in: Badescu, V., Gueymard, C.A., Cheval, S., Oprea, C., Baciu, M., Dumistrescu, A.,
# Iacobescu, F., Milos, I., Rada, C., 2013. Accuracy analysis for fifty-four clear-sky solar 
# radiation models using routine hourly global irradiance measurements in Romania. Renewable Energy 55, 85–103.

HS <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  Bn <- 1098*exp(-0.057/cosThzS)
  B0 <- Bn*cosThzS
  D0 <- 94.23*(cosThzS^0.5)
  G0 <- B0 + D0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}
  

## 67. Janjai et al. 2011 CHECKED with original

# Janjai, S., Sricharoen, K., Pattarapanitchai, S., 2011. Semi-empirical models for the
# estimation of clear sky solar global and direct normal irradiances in the tropics. Applied Energy 88(12), 4749-4755.

JSP2011 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  m <- data$m
  p <- data$p
  uo <- data$uo
  pwc <- data$pwc
  alpha <- data$alpha
  beta <- data$beta
   
  am <- m*p/1013.25
  foo1 <- -0.102073 + 1.421569*beta + 0.080228*alpha + 0.011522*pwc + -0.208154*uo
  foo2 <- 0.70*(Bo0/cosThzS)*(cosThzS^0.25)
  Bn <- foo2*exp(-foo1*am)
  Bn[Bn>1000] <- 0
  B0 <- Bn*cosThzS
  D0 <- (Bo0/cosThzS)*(0.18738 + 0.375057*beta + 0.011316*alpha + 
                         0.006227*pwc - 0.51979*uo)*(cosThzS^0.75)
  D0[D0<0] <- 0
  G0 <- B0 + D0
  
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}


## 68. Dazhi et al. CHECKED with original

# Dazhi, Y., Jirutitijaroen, P., Walsh, W.M., 2012. The Estimation of Clear Sky Global Horizontal Irradiance at
# the Equator. Energy Procedia 25, 141-148.

DJW2012 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  y <- unique(year(index(data)))
  julianday <- as.numeric(julian(as.Date(index(data)), origin = as.Date(paste0(y,'-01-01')))) + 1
  # day angle
  day_angle<- julianday*2.0*pi/365.2422
  # Correction accounting for eccentricity
  eo <- 1 + 0.03344*cos(day_angle - 2.80*pi/180.0)
  
  G0 <- 0.8277*(1366.1*eo)*(cosThzS^1.3644)*exp(-0.0013*(90 - theta))
  
  D0 <- B0 <- Bn <- zoo(rep(NA, length(G0)), index(G0))
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}

## 69. Dai and Fang CHECKED with original

# Dai, Q., Fang, X., 2014. A simple model to predict solar radiation under clear sky conditions.
# Advances in Space Research 53, 1239–1245.

DF2014 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  Bo0 <- data$Bo0
  
  p <- data$p
  pwc <- data$pwc
  aod500 <- data$aod500
  
  m <- (cosThzS +  0.15*(3.885 + (90 -theta))^-1.253)^-1
  am <- m*p/1013
  tau <- 0.744*aod500
  Bn <- (Bo0/cosThzS)*exp(-0.103*am^0.571 - 0.081*(pwc*m)^0.213 - tau^0.91*m^0.87)
  B0 <- Bn*cosThzS
  D0 <- (0.143 + 0.113*cosThzS - 0.0485*pwc + tau)*((Bo0/cosThzS) - Bn)*cosThzS
  G0 <- B0 + D0
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}

## 70. Badescu & Dumitrescu 2014

# Badescu, V., Dumitrescu, A., 2014. Simple models to compute solar global irradiance 
# from the CMSAF product Cloud Fractional Coverage. Renewable Energy 66, 118-131.


BD2014 <- function(data, loc){
  
  cosThzS <- data$cosThzS
  theta <- r2d(acos(cosThzS))
  
  
  G0 <- 1007.31*cosThzS - 74.09
  G0[G0<0] <- 0
  D0 <- B0 <- Bn <- zoo(rep(NA, length(G0)), index(G0))
  
  night <- (theta >= 90)
  B0[night] <- 0
  Bn[night] <- 0
  D0[night] <- 0
  G0[night] <- 0
  
  
  z <- cbind(G0, D0, B0, Bn) 
  coredata(z)[is.na(z)] <- 0
  z
}

