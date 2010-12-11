function [num,name,vote,cand] = getCand(filter,comb_method,query_num)

% [num,name,vote,cand] = getCand(filter,comb_method,query_num)
% 
% First coded 9 Dec 2010 by Aaron Hallquist
% Latest revision 9 Dec 2010 by Aaron Hallquist
% 
% K = number of queries
% 
% DESCRIPTION
%   This function reads the vote results for all queries and generates the
%   candidate list for post-processing, using the query set number,
%   combination method, and candidate filter as inputs.
% 
% INPUTS
%   filter:         String listing the filter method that determines the
%                   candidate list from the flann vote results
%   comb_method:    String listing the combination method that determines
%                   the vote results from flann.
%   query_num:      Integer that determines the query set from which to
%                   generate results
% 
% OUTPUTS
%   num:        Kx1 array of query numbers
%   name:     	Kx1 cell array of query filenames
%   vote:     	Structure containing the flann vote results...
%       .method:    String containing the combination method for votes
%       .files:     Kx1 cell array of cell arrays of strings containing the
%                   database image filenames associated with votes
%       .values:    Kx1 cell array of numerical arrays containing the vote
%                   results for each database image and each query
%   cand:       Structure containing the candidate images and scores...
%       .filter:    String containing the filter method
%       .files:     Kx1 cell array of cell arrays of strings containing the
%                   candidate database image filenames
%       .votes:     Kx1 cell array of numerical arrays containing the flann
%                   vote results for each candidate image
%       .scores:    Kx1 cell array of numerical arrays containing the score
%                   for each candidate image from post-processing
%                   (these are arrays of zeros coming out of this function)

% Get vote files
vDir = ['.\',comb_method,'\query',num2str(query_num),'\'];
vote_files = dir(vDir);
vote_files = vote_files(3:end);
nq = length(vote_files); % number of queries

% Parse vote files and set num, name, vote, cand
num = zeros(nq,1);
name = cell(nq,1);
vote.method = comb_method;
vote.files = cell(nq,1);
vote.values = cell(nq,1);
cand.filter = filter;
cand.files = cell(nq,1);
cand.votes = cell(nq,1);
cand.scores = cell(nq,1);
for k=1:nq
    % Set query number and names
    fn = vote_files(k).name;
    num(k) = str2num(fn(5:8));
    name{k} = fn(1:end-12);
    % Get the vote information for query k
    [kVotes,kFiles] = parseVotes(fn,vDir);
    vote.files{k} = kFiles;
    vote.values{k} = kVotes;
    % Apply candidate filter to query k vote results
    if strcmp(filter,'upto10')
        minVote = kVotes(10);
        minVoteIdx = find(kVotes==minVote,1,'last');
        if minVoteIdx > 10
            topN = find(kVotes==minVote,1,'first') - 1;
        else
            topN = 10;
        end
        cand.files{k} = kFiles(1:topN);
        cand.votes{k} = kVotes(1:topN);
    elseif strcmp(filter,'above1')
        idx = find(kVotes>1);
        cand.files{k} = kFiles(idx);
        cand.votes{k} = kVotes(idx);
    end
end

end