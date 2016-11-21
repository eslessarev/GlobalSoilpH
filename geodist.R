# geodist()
# Home-made Haverinse Formula geographic distance function
# Accompanies Manuscript: "Water Balance Creates a Threshold in Soil pH at the Global Scale"
# Eric Slessarev 11/20/2016

# This function accepts a single target point and a set of surrounding points
# retunrs all of the distances between the target point and the surrounding points.
# Latitude and longitude are assumed to be in decimal degree format

geodist = function(targetlon,targetlat,searchlon,searchlat){
  
  # targetlon and targetlat are the coordinates of the target point
  # searchlon and searchlat are multiple coordinates to be searched  

  #radius of the earth
  re = 6378.1
  
  #cell centers in radians
  targetlat = targetlat * (2*pi/360)
  targetlon = targetlon * (2*pi/360)
  searchlat = searchlat * (2*pi/360)
  searchlon = searchlon * (2*pi/360)
  
  # haversine forumla
  t1 = sin((targetlat-searchlat)/2)^2
  t2 = cos(targetlat)*cos(searchlat)*sin((targetlon-searchlon)/2)^2
  d = 2*re*asin(sqrt(t1+t2))
  
  return(d)
}