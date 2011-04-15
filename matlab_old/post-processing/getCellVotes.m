function [image,vote] = getCellVotes(query,lat,long,vote_dir,max_dist)

% [image,vote] = getCellVotes(query_num,lat,long,vote_dir,max_dist)
% 
% First coded 8 Jane 2011 by Aaron Hallquist.
% Latest revision 8 Jan 2011 by Aaron Hallquist.
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