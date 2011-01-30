function [images,features] = getFeatures(query_num,vote_dir)

% [features] = getFeatures(query_num,query_set)
% 
% First coded 12 Dec 2010 by Aaron Hallquist.
% Latest revision 6 Jan 2011 by Aaron Hallquist.
% 
% DESCRIPTION:
%   This function reads vote totals for potential database matches. For
%   every photo which receives votes in at least one cell search, this
%   gathers the three largest vote totals of that photo (zeros if not 3).
%       
% 
% INPUT:
%   query_num:  The number of the query to analyze
%   vote_dir:   Directory with vote results
% 
% OUTPUT:
%   images:     cell array of image filenames
%   features:   matrix containing the top 3 image vote totals

addpath E:\Research\app\code\matlab\util\

query_files = dir(vote_dir);
query_files = strvcat(query_files.name);
file_idx = strmatch(num2str(query_num),query_files(:,5:8));
query_files = cellstr(query_files(file_idx,:));
ncells = length(query_files);

photos = cell(0,1);
votes = zeros(0,1);
for k=1:ncells
    [v,p] = parseVotes(query_files{k},vote_dir);
    photos = [photos;p];
    votes = [votes;v];
end

images = cell(0,1);
features = zeros(0,3);
while ~isempty(photos)
    images = [images;photos(1)];
    idx = strcmp(photos,photos(1));
    v = sort(votes(idx),'descend')';
    if length(v)<3
        v(end+1:3) = 0;
    end
    features = [features;v(1:3)];
    photos(idx) = [];
    votes(idx) = [];
end