# cellarea()
# Function to get area of a grid cell given latitude and resolution 
# Accompanies Manuscript: "Water Balance Creates a Threshold in Soil pH at the Global Scale"
# Eric Slessarev 11/20/2016

cellarea = function(lat,resolution){

#lat is the latitude at the cell center

#lat1 is the latitude of the lower bound of the cell
lat1 = lat-resolution/2

#lat2 is the latitude of the upper bound of the cell
lat2 = lat+resolution/2

re = 6378.1 # radius of the earth

# vertical distance from lat1 to the pole
h1 = re - sin(lat1*pi/180)*re

# vertical distance from lat2 to the pole
h2 = re - sin(lat2*pi/180)*re

#area of a spherical cap with base at lat1
A1 = 2*pi*re*h1

#area of a spherical cap with base at lat2
A2 = 2*pi*re*h2

#area of the row of grid cells at latitude lat
Ar = A1 - A2

#area of one grid cell at latitude lat
Ac = Ar*(resolution/360)

return(Ac)}