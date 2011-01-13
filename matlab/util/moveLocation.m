function [lat,lon] = moveLocation(lat,lon,bear,dist)

% [lat,lon] = moveLocation(lat,lon,bear,dist)
% 
% DESCRIPTION
%   This function computes a new location in latitude and longitude from an
%   old location with a bearing and distance (in meters) to move. The
%   bearing is in degrees, with 0 indicating north and 90 indicating east.

radius = 6371000;
lat = pi/180 * lat;
lon = pi/180 * lon;
bear = pi/180 * bear;

lat = asin( sin(lat)*cos(dist/radius) + ...
    cos(lat)*sin(dist/radius)*cos(bear) );
lon = lon + atan2( sin(bear)*sin(dist/radius)*cos(lat) , ...
    cos(dist/radius) - sin(lat)*sin(lat) );

lat = 180/pi * lat;
lon = 180/pi * lon;