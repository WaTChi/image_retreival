function [d] = latlonDistance(lat1,lon1,lat2,lon2)

% [d] = latlonDistance(lat1,lon1,lat2,lon2)
% 
% DESCRIPTION
%   This function returns the distance between latitude and longitude
%   points. Matrix inputs are okay. The output d is in meters.

radius = 6371000;
lat1 = pi/180 * lat1;
lat2 = pi/180 * lat2;
lon1 = pi/180 * lon1;
lon2 = pi/180 * lon2;

d = radius * acos( sin(lat1).*sin(lat2) + ...
    cos(lat1).*cos(lat2).*cos(lon1-lon2) );