function [image,vote,img_lat,img_lon,nfeats] = ...
    cellCombine(cell_files,vote_dir,ncand)

% [image,vote,img_lat,img_lon,nfeats] = ...
%     cellCombine(cell_files,vote_dir,ncand)
% 
% DESCRIPTION
%   This function combines vote totals across the cell results listed in
%   cell_files and compiles a list of the image locations.

im = [];
vo = [];

% Gather cell vote tallies
ncells = length(cell_files);
for k=1:ncells
    cell = cell_files{k};
    [v,i] = parseVotes(cell,vote_dir);
    im = [im;i];
    vo = [vo;v];
end

% Combine image vote totals
image = [];
vote = [];
while ~isempty(im)
    image = [image;im(1)];
    idx = strcmp(im,im(1));
    try
        vote = [vote;sum(vo(idx))];
    catch e
        
        find(idx)
        throw(e)
    end
    im(idx) = [];
    vo(idx) = [];
end

% Sort in descending order and normalize
[vote,sort_idx] = sort(vote,'descend');
image = image(sort_idx);
nfeats = sum(vote) / ncells; % number of features in query image
vote = vote / nfeats;

% Get lat/lon of the top N images
if nargin < 3
    N = length(vote);
else
    N = min(length(vote),ncand);
end
vote = vote(1:N);
image = image(1:N);
img_lat = zeros(N,1);
img_lon = zeros(N,1);
for k=1:N
    idx = strfind(image{k},',');
    img_lat(k) = str2double(image{k}(1:idx-1));
    img_lon(k) = str2double(image{k}(idx+1:end-13));
end