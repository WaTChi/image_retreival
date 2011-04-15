function [lats,lons] = moveLocation(lat,lon,bear,dist)

% [lat,lon] = moveLocation(lat,lon,bear,dist)
% 
% DESCRIPTION
%   This function computes a new location in latitude and longitude from an
%   old location with a bearing and distance (in meters) to move. The
%   bearing is in degrees, with 0 indicating north and 90 indicating east.
%   Matrix inputs for bear and dist are accepted.

radius = 6371000;
lat = pi/180 * lat;
lon = pi/180 * lon;
bear = pi/180 * bear;

lats = asin( sin(lat)*cos(dist/radius) + ...
    cos(lat)*sin(dist/radius).*cos(bear) );
lons = lon + atan2( sin(bear).*sin(dist/radius)*cos(lat) , ...
    cos(dist/radius) - sin(lat)*sin(lats) );

lats = 180/pi * lats;
lons = 180/pi * lons;