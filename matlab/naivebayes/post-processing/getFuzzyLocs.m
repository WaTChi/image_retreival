function [lats,lons] = getFuzzyLocs(lat,lon,rad,spacing)

% [lats,longs] = getFuzzyLocs(lat,long,fuzzy_rad)
% 
% DESCRIPTION
%   This function generates a set of latitudes and longitudes on a grid
%   (where grid spacing is 1 meter) using an original location and a
%   maximum radius.

grid = -floor(rad):spacing:floor(rad);
maxN = length(grid)^2;
lats = zeros(maxN,1);
lons = zeros(maxN,1);
last = 0;
for x = grid
    for y = grid
        d = sqrt(x^2+y^2);
        if d < rad
            last = last+1;
            b = 180/pi*atan2(y,x);
            [lats(last),lons(last)] = moveLocation(lat,lon,b,d);
        end
    end
end
lats = lats(1:last);
lons = lons(1:last);