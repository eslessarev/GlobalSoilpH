# Soil pH Buffer Calculation Script
# Accompanies Manuscript: "Water Balance Creates a Threshold in Soil pH at the Global Scale"
# Eric Slessarev 11/20/2016


#============================================================================
# INPUT DATA DESCRIPTION
# the aluminum buffer model is empirical, and so requires data
# use the spatially resampled NCSS data, "ssn"

# required fields:
# pH: pH in 1:1 water/soil slurry
# ECEC: effective cation exchange capacity (cmolc/kg)
# EXAL: exchangeable aluminum (cmolc/kg)


# example:

# cellindex   DS      PID       cclon    cclat    pH    AWB_PM  ...  ECEC  EXAL
#   2906     NCSS  NCSS24576    -163.5    64.5    6.7    3.29   ...   NA    NA
#   2906     NCSS  NCSS24577    -163.5    64.5    6.6    3.29   ...   NA    NA
#   3082     NCSS  NCSS33677    -162.5    68.5    7.1    2.89   ...  20.3  12.9
#    ...

#load("C:/Users/.../soildata.RData")

#============================================================================


#============================================================================
#load rootSolve
require(rootSolve)
#============================================================================


#============================================================================
# CALCITE BUFFER

# parameters
K1 = 10^(-6.35) # first Ka for bicarbonate (mol L-1)
K2 = 10^(-10.33) # second Ka for bicarbonate (mol L-1)
Ks = 10^(-8.48) # solubility product for calcite (mol2 L-2)
kH = 0.033 # Henry's constant (mol L-1 atm-1)
Kw = 10^-14 # dissociation of water (mol2 L-2)
pCO2_lab = 0.000345 # pCO2 in 1985, median collection date for soil profile data (atm)


# quartic function
calcpH_calcite = function(H,pCO2,parameters){
  
  # unpack parameters
  for(p in 1:length(parameters)){
    assign(names(parameters)[p],parameters[p])
  }
  
  # terms 1-4
  t1 = (2*Ks/(K1*K2*kH*pCO2))*H^4 
  t2 = H^3
  t3 = (Kw+K1*kH*pCO2)*H 
  t4 = 2*(K1*K2*kH*pCO2)
  
  # function for H+
  H = t1  + t2 - t3 - t4
  return(H)
  
}

# bundle parameters
Ca_parameters = c(K1=K1,K2=K2,Ks=Ks,kH=kH,Kw=Kw)

# Get pH of open calcite system under lab conditions
HCa = uniroot.all(calcpH_calcite,interval=c(10^-8,10^-8.5),pCO2=pCO2_lab,parameters=Ca_parameters)

# convert to pH
pHCa = -1*log10(HCa)

# round
pHCa = round(pHCa,1)
#============================================================================


#============================================================================
# Al EXCHANGE BUFFER

# use spatially resampled subsoil horizons (0.5 m depth) from the NCSS dataset

# NOTE:
# Several restrictions apply before using any of the models below:
# values of exchangeable aluminum and calcium must be nonzero, or the models are undefined
# thus, values for EXAL and ECEC must not equal zero, nor may they equal each other.
# missing data must be excluded.

exData = ssn[ssn$EXAL > 0 & 
               ssn$ECEC > 0 &
               !is.na(ssn$EXAL) &
               !is.na(ssn$ECEC) &
               !ssn$EXAL >= ssn$ECEC,]

# extract values for exchangeable Al and Ca, plus pH
exCa = exData$ECEC - exData$EXAL # exchangeable Ca values
exAl = exData$EXAL # exchangeable Al values
pHobs = exData$pH # pH values

# the Gaines-Thomas exchange model (as adapted by Ruess et al. 1990)
GT = log10((exCa^3)/(exAl^2))
GT_Model = summary(lm(pHobs~GT))

# the Gapon exchange model (as adapted by Ruess et al. 1990)
Gapon = log10((exCa)/(exAl))
Gapon_Model = summary(lm(pHobs~Gapon))

# the Gapon model fits better (contrary to what Ruess et al. found)
# get the fitted Gapon coefficients 
b1 = Gapon_Model$coef[2,1]
b0 = Gapon_Model$coef[1,1]  

# it is simplest to just take the mean pH of soils with exchangeable aluminum
pHAl = mean(pHobs)

# round
pHAl = round(pHAl,1)
#============================================================================
