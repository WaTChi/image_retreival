function [qVote,qPhoto] = parseVotes(filename,vote_dir)

% [qVote,qPhoto] = parseVotes(filename,vote_dir)
% 
% First coded 9 Dec 2010 by Aaron Hallquist
% Latest revision 9 Dec 2010 by Aaron Hallquist
% 
% K = number of queries
% 
% DESCRIPTION:
%   This function parses the vote results from flann for the given query
%   and returns the database files and associated votes.
% 
% INPUTS
%   filename:   Filename for query's vote data
%   vote_dir:   Directory where vote data is stored
% 
% OUTPUTS
%   qVote:      Array of vote values
%   qPhoto:     Cell array of database files

% Read and parse vote data
vote_data = textread([vote_dir,filename],'%s');
vote_data = reshape(vote_data,[2,length(vote_data)/2])';
qVote = str2num(strvcat(vote_data(:,1)));
qPhoto = vote_data(:,2);

end