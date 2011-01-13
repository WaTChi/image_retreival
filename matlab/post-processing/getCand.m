function [candIdx] = getCand(lat,lon,img_lat,img_lon,vote,filt_dist,exN)

% [candIdx] = getCand(lat,lon,img_lat,img_lon,vote,filt_dist,exN)
% 
% DESCRIPTION
%   This function gets up to 10 candidate images (eliminating ties) by
%   looking at combined vote totals and applying a distance filter.
% 
% INPUTS
%   lat:        Latitude of reported location
%   lon:        Longitude of reported location
%   img_lat:    List of image latitude coordinates (descending votes)
%   img_lon:    List of image longitude coordinates (descending votes)
%   vote:       List of vote totals for the images (descending order)
%   filtDist:   Distance in meters of removal filter
%   exN:        Boolean that decides between exactly N and up to 10
% 
% OUTPUTS
%   candIdx:    Indices in [image,vote] which are the candidates

N = 10; % maximum number of candidate images
candIdx = zeros(N+1,1);
n = 0;
k = 1;
while n < N+1
    img_dist = latlonDistance(lat,lon,img_lat(k),img_lon(k));
    if img_dist < filt_dist
    	n = n+1;
        candIdx(n) = k;
    end
    k=k+1;
end

if exN
    candIdx = candIdx(1:N);
else
    if vote(candIdx(N)) == vote(candIdx(N+1))
        n = find( vote(candIdx) == vote(candIdx(N)) , 1 , 'first' ) - 1;
        candIdx = candIdx(1:n);
    else
        candIdx = candIdx(1:N);
    end
end