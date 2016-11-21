# PET calculation script
# Accompanies Manuscript: "Water Balance Creates a Threshold in Soil pH at the Global Scale"
# Eric Slessarev 11/20/2016

# This script calculates PET using a modified Penman-Monteith-Leuning model
# It also obtains Priestley-Taylor PET using two different radiation datasets

#============================================================================
# INPUT DATA DESCRIPTION

# Input data are stored as global arrays of monthly climatological means
# The geographic extent of the data is global: -90 90 -180 180
# the resolution of the data is 1 decimal degree.
# each array is 360 (lon) X 180 (lat) X 12 (months)
# cell centers are -179.5:179.5 (lon) -89.5:89.5 (lat) 1/15:12/15 (months)
# Note that the maps are rotated 90 degrees clockwise
# i.e. the south pole comprises the first column of the data
# this is for compatability with the image() command in R
# The input data must be spatially and temporally aggregated upstream of this script.

# monthly_tmp: CRU Ts3.13 mean air temperature (degrees C)
# monthly_vap: CRU Ts3.13 mean vapor pressure (kPa)
# monthly_Rn_CERES: CERES net surface radiation (MJ m-2 d-1)
# monthly_SW_CERES: CERES surface shortwave downwelling radiation (MJ m-2 d-1)
# monthly_Rn_GEWEX: GEWEX net surface radiation (MJ m-2 d-1)
# monthly_LAI: GLASS LAI dataset (m2 m-2)
# monthly_COV: MOD12 land cover classes, repeated for each month*
# DEM: ETOPO1 digital elevation model, repeated each month* (m)

# * these static data are stored in 260X180X12 arrays for compatability.  
#============================================================================


#============================================================================
# PARAMETERS

# general parameters, from Allen et al. 1998
l = 2.45 # latent heat of vaporization MJ kg-1
cp = 0.001013 #specific heat of air Mj kg-1
eps = 0.622 # ratio MW water vapor/dry air 
R = 0.287058 # specific gas constant kJ kg-1 K-1
a = 1.26 # Priestley-Taylor constant
std_atm = 101.3 # standard atmospheric pressure (kPa)
abs_zero = -273.16 # absolute zero (degrees C)

# parameters for PML model
Gs_max = 0.006 # maximum stomatal conductance, Kelliher et al. 1995 (m s-1)
ka = 0.6 # extinction coefficient for radiation, Zhang et al. 2008 (m-1)
kQ = 0.6 # extinction coefficient for light, Zhang et al. 2008 (m-1)
Q50 = 2.6 # half saturation constant for light, Zhang et al. 2008 (MJ m-2 d-1)
D50 = 0.8 # half saturation constant for vapor pressure deficit, Zhang et al. 2008 (kPa)

# Biome lookup table for Ga:
#____________________________
# biome       code  value
#____________________________
# Water        0    NA
# ENF          1    0.033
# EBF          2    0.033
# DNF          3    0.033
# DBF          4    0.033
# MF           5    0.033
# CSH          6    0.0125
# OSH          7    0.0125
# WL           8    0.033
# SV           9    0.033
# Grass        10   0.01
# Wetland      11   0.01
# Crop         12   0.01
# Urban        13   NA
# Crop_mosaic  14   0.01
# Ice          15   NA
# Barren       16   0.01
# Unclass      254  NA
# Fill         255  NA
#____________________________

biome = c("Water","ENF","EBF","DNF","DBF","MF","CSH","OSH","WL","SV","Grass",
          "Wetland","Crop","Urban","Crop_mosaic","Ice","Barren","Unclass","Fill")
code = c(0,1,2,3,4,5,6,7,8,9,10,11,12,
         13,14,15,16,254,255)
value = c(NA,0.033,0.033,0.033,0.033,0.033,0.0125,0.0125,0.033,0.033,0.01,
          0.01,0.01,NA,0.01,NA,0.01,NA,NA)
Ga_BLUT = data.frame(biome,code,value)
#============================================================================


#============================================================================
# FUNCTIONS

# vt() obtains saturation vapor pressure (kPa)
# Allen et al. 1998 equation 14
vt = function(Ta){
  vt = 0.6108*exp((17.27*Ta)/(Ta+237.3))
  return(vt)
}

# Delta() obtains the slope of the saturation vapor pressure curve (kPa C-1)
# Allen et al. 1998 equation 20
Delta = function(Ta){
  vt = 0.6108*exp((17.27*Ta)/(Ta+237.3))
  D = (4098*vt)/(Ta+237.3)^2
  return(D)
}

# atmospheric pressure, from elevation (kPa)
# Allen et al. 1998 equation 7
Patm = function(elevation,std_atm){
  P = std_atm*((293-0.0065*elevation)/293)^5.26
  return(P)
}

# function for assigning parameter values from a biome lookup table
assign_param =  function(BLUT,monthly_COV){
  param = array(NA,dim=c(360,180,12))
  for(row in 1:nrow(BLUT)){
    param[monthly_COV == BLUT$code[row]] = BLUT$value[row] 
  }
  return(param)
}
#==============================================================================



#==============================================================================
# DERIVED DATA

# these data are obtained using the parameters and input data above
# they are stored in 360X180X12 arrays just like the input data

# air temperature in Kelvin (degrees K)
monthly_tmpK = monthly_tmp - abs_zero

# saturation vapor pressure (kPa)
monthly_vst = vt(monthly_tmp)

# vapor pressure deficit (kPa)
monthly_vpd = monthly_vst - monthly_vap

# atmospheric pressure (kPa)
P = Patm(DEM,std_atm)

# density of dry air (kg m-3)
pa = P/(monthly_tmpK*R) 

# slope of vapor pressure curve (kPa C-1)
D = Delta(monthly_tmp) 

# psychrometric constant
y = (cp*P)/(eps*l)

# aerodynamic conductance, Zhang et al. 2008 (m s-1)
Ga = assign_param(Ga_BLUT,monthly_COV)
#==============================================================================


#==============================================================================
# PET CALCULATIONS

#global Priestley-Taylor PET, CERES data
monthly_PT = (1/l)*a*(D*monthly_Rn_CERES)*(1/(D+y)) # (mm d-1)
monthly_PT = monthly_PT*30.41667 #scale days to months (mm month-1)
annual_PT = apply(monthly_PT,c(1,2),sum) # get annual sum as 360X180X1 array (mm year-1)

#global Priestley-Taylor PET, GEWEX data
monthly_PT_gw = (1/l)*a*(D*monthly_Rn_GEWEX)*(1/(D+y)) # (mm d-1)
monthly_PT_gw = monthly_PT_gw*30.41667 # scale days to months (mm month-1)
annual_PT_gw = apply(monthly_PT_gw,c(1,2),sum) # get annual sum as 360X180X1 array (mm year-1)

# Modified Penman-Monteith-Leuning PET
# Soil evaporation is accounted for by Priestly-Taylor model.

# Soil available energy, Leuning et al. 2008 (MJ m-2 d-1)
As = monthly_Rn_CERES*exp(-ka*monthly_LAI)

# Canopy available energy, Leuning et al. 2008 (MJ m-2 d-1)
Ac = monthly_Rn_CERES - As 

# Photosynthetically active radiation, Leuning et al. 2008 (MJ m-2 d-1)
Qh = (1/2)*monthly_SW_CERES*(1 - exp(-ka*monthly_LAI)) 

# Canopy conductance, from Leuning et al. 2008 (m s-1)
Gc = (Gs_max/ka)*log((Qh+Q50)/(Qh*exp(-ka*monthly_LAI) + Q50))*(1/(1+monthly_vpd/D50))

# convert conductances from m s-1 to m d-1
Gc = Gc*60*60*24
Ga = Ga*60*60*24

# canopy evaporation, Penman Monteith equation (mm d-1)
Ec = (1/l)*(Ac*D + cp*pa*monthly_vpd*Ga)*(1/(D + y*(1 + Ga/Gc)))

# soil evaporation, Priestley-Taylor method (mm d-1)
Es = (1/l)*a*(As)*(D/(D+y)) 

monthly_PML = Ec + Es # sum canopy and soil evaporation (mm d-1)
monthly_PML = monthly_PML*30.41667 # scale daily values to monthly values (mm month-1)
annual_PML = apply(monthly_PML,c(1,2),sum) # get annual sum as 360X180X1 array (mm year-1)
#==============================================================================


#==============================================================================
# REFERENCES

# Kelliher, F.M., Leuning, R., Raupach, M.R., Schulze, E.D. 
#    Maximum conductances for evaporation from global vegetation types. 
#    Agri. For. Meteorol. 73, 1-16 (1995).

# Leuning, R., Zhang, Y.Q., Rajaud, A., Cleugh, H., Tu, K. 
#    A simple surface conductance model to estimate regional evaporation 
#    using MODIS leaf area index and the Penman Monteith equation. 
#    Water Resour. Res. 44, 1-17 (2008).

# Zhang, Y.Q., Chiew, H.S., Zhang, L., Leuning, R., Cleugh, H. 
#    Estimating catchment evaporation and runoff using MODIS leaf area index 
#    and the Penman Monteith equation. Water Resour. Res. 44, (2008).

#==============================================================================
