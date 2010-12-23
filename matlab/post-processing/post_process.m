function [results] = post_process(filter,comb_method,query_num)

% [results] = post_process(filter,comb_method,query_num)
% 
% First coded 9 Dec 2010 by Aaron Hallquist
% Latest revision 9 Dec 2010 by Aaron Hallquist
% 
% DESCRIPTION:
%   This file takes in a combination method and query set number and runs
%   post-processing on the flann vote results.
% 
% INPUTS
%   filter:         String choosing filter for obtaining candidate images
%   comb_method:    String choosing the combination method, one of...
%           'no-comb':      Nearest cell only
%           'dup-comb':     Finds duplicate images across multiple cells
%           'vote-comb':    Sums vote results across multiple cells
%   query_num:      Integer choosing the query set to use
% 
% OUTPUTS
%   results:  	Structure containing the details of post-processing
%       - type 'help ppStruct' for details about its contents

% Fixed directories
qDir = 'E:\Research\collected_images\query\';
gtDir = 'E:\Research\app\code\matlab\ground-truth\';
dbDir = 'E:\Research\collected_images\earthmine-new,culled\37.871955,-122.270829\';

% Adjust directories based on inputs
qDir = [qDir,'query',num2str(query_num),'\'];
gtDir = [gtDir,'query',num2str(query_num),'\'];

% Post-processing code initialization
run 'C:\vlfeat-0.9.9\toolbox\vl_setup'

% Code dependencies
addpath('.\..\util\')

% Initialize results structure
results.directory = qDir;
[ results.num , results.name , results.vote , results.cand ] = ...
    getCand(filter,comb_method,query_num);
nq = length(results.num); % number of queries
results.numqueries = nq;

% Get and store the ground truth files
results.gtruth = cell(nq,1);
if query_num==3
    [gtIdx,gtFile] = parseGT('GYRBO',gtDir);
else
    [gtIdx,gtFile] = parseGT('A',gtDir);
end
for k=1:nq
    idx = find(gtIdx==results.num(k));
    results.gtruth{k} = gtFile( gtIdx(idx,2) : gtIdx(idx,3) );
end

% compute the post-processing scores for the candidate images
for k=1:nq
    
    disp(['Recombination on query ',num2str(k)])
    
    querySift = [results.name{k},'sift.txt'];
    [~,da] = vl_ubcread(strcat(qDir,querySift));
    
    candFiles = results.cand.files{k};
    nc = length(candFiles); % number of candidates
    scores = zeros(nc,1);
    for j=1:nc
        candSift = [candFiles{j},'sift.txt'];
        [~,db] = vl_ubcread(strcat(dbDir,candSift));
        scores(j) = length(vl_ubcmatch(da,db));
    end
    
    results.cand.scores{k} = scores;
    
end

outFile = ['.\',comb_method,'\query',num2str(query_num),'_results.mat'];
save(outFile,'results');

end