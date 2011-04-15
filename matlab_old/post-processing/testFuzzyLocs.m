lat=37;
lon=-122;

[lats,lons] = getFuzzyLocs(37,-122,200,4);

N = length(lats);
d = latlonDistance( repmat(lat,[N,1]),repmat(lon,[N,1]),lats,lons );

lat = pi/180 * lat;
lon = pi/180 * lon;
lats = pi/180 * lats;
lons = pi/180 * lons;
dLon = lons-lon;

y = sin(dLon) .* cos(lats);
x = cos(lat)*sin(lats) - sin(lat)*cos(lats).*cos(dLon);
b = atan2(y,x);

x = d.*cos(b);
y = d.*sin(b);

plot(x,y,'.')