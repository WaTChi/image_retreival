function [cand,cand_vote] = getCand(lat,lon,img_lat,img_lon,image,vote,filter,exN)

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
%   image:      List of images (descending vote order)
%   vote:       List of vote totals for the images (descending order)
%   filter:   	Type of filter: 'none' or 'cutoff' or 'exponential'
%   exN:        Boolean that decides between exactly N and up to 10
% 
% OUTPUTS
%   candIdx:    Indices in [image,vote] which are the candidates

N = 10; % maximum number of candidate images
cand = cell(11,1);
cand_vote = zeros(N+1,1);
filt_vote = zeros(N+1,1);

nimg = length(vote);
k=1;
minf_vote = 0;
minf_idx = 1;
while k <= nimg && vote(k) > minf_vote
    img_dist = latlonDistance(lat,lon,img_lat(k),img_lon(k));
    if strcmp(filter,'cutoff')
        w = (img_dist<100);
    elseif strcmp(filter,'exponential')
        w = exp(-img_dist/50);
    else
        w = 1;
    end
    if w*vote(k) > minf_vote
        filt_vote(minf_idx) = w*vote(k);
        cand_vote(minf_idx) = vote(k);
        cand(minf_idx) = image(k);
        [minf_vote,minf_idx] = min(filt_vote);
    end
    k = k+1;
end

[cand_vote,perm] = sort(cand_vote,'descend');
cand = cand(perm);

if exN
    cand = cand(1:N);
    cand_vote = cand_vote(1:N);
else
    if cand_vote(N) == cand_vote(N+1)
        n = find( cand_vote == cand_vote(N) , 1 , 'first' ) - 1;
        candIdx = candIdx(1:n);
    else
        candIdx = candIdx(1:N);
    end
end