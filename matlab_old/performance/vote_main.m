% Script for running analysis

vote_dir = 'E:\cellsg=50,r=100,d=86.6,query3kdtree4matches,threshold=70k\';
gt_dir = 'E:\Research\app\code\matlab\';
match_file = 'E:\Research\g=100,r=d=236.6cells,query3matches.txt';

exclude = [];

mFilter = [ 10 , 0.3 , 8 , 0.5 ];
% mFilter(1) = maximum number of images retained for matching
% mFilter(2) = ratio threshold for matching
% mFilter(3) = minimum number of absolute votes for matching
% mFilter(4) = minimum number of votes as a percentage of total votes

qFilter = [ 1.5 ];
% qFilter(1) = minimum percent of votes the maximum vote attains for query

disp(' ')
disp('Running ground truth analysis...')
[goodMask,noMatch,tmp] = analyzeGT(100,exclude,mFilter,vote_dir,gt_dir);
exclude = [exclude,noMatch];
disp(' ')
disp('Running vote analysis...')
filtQ = analyzeVotes(200,exclude,goodMask,qFilter,vote_dir);
disp(' ')
disp('Running filtered ground truth analysis...')
[~,~,tmp2] = analyzeGT(300,filtQ,mFilter,vote_dir,gt_dir);
disp(' ')

elim = [];
for k = 1:length(tmp2)
    elim = [elim,find(tmp==tmp2(k))];
end
tmp(elim) = [];

tmp