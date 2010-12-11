function [outputs] = post_process(query_num,comb_method)

% [outputs] = post_process(comb_method)
% 
% First coded 9 Dec 2010 by Aaron Hallquist
% Latest revision 9 Dec 2010 by Aaron Hallquist
% 
% DESCRIPTION:
%   This file takes in a
% 
% INPUTS
%   divs:       String containing the ground truth divisions we want to
%               inspect. One character for each division out of:
%                   - G : Green division
%                   - Y : Yellow division
%                   - R : Red division
%                   - B : Blue division
%                   - O : Orange division
%               An example input is divs = 'GY' to include only the green 
%               and yellow divisions.
%   root_dir:   This is the directory where the ground truth files reside.
% 
% OUTPUTS
%   gtIdx:      Ground truth index. This is a Kx3 matrix containing...
%                   gtIdx(:,1): vector of query numbers
%                   gtIdx(:,2): pointer to the start of database matches
%                   gtIdx(:,3): pointer to the end of database matches
%   gtFile:     Ground truth database match files. This is a cell array of
%               strings of query matches. It is set up so that
%                   gtFile( gtIdx(k,2) : gtIdx(k,3) )
%               refers to the database matches for query #gtIdx(k,1).
%               Examples: 
%                   gtIdx(k,:) = [ 8847 13 12 ]
%                       Query #8847 has no matches in the divisions
%                       specified [ gtFile( 13:12 ) is empty ]
%                   gtIdx(j,:) = [ 8860 65 68 ]
%                       Query #8860 has 4 matches in the divisions
%                       specified [ the entries of gtFile( 65:68 ) ]

% Combination technique to benchmark
comb = 'no-comb';

% Directories
gt_dir = 'E:\Research\app\code\matlab\ground-truth\';
dbDir = 'E:\Research\collected_images\earthmine-new,culled\37.871955,-122.270829\';
qDir = 'E:\Research\collected_images\query\query3\';

% Post-processing code initialization
run 'C:\vlfeat-0.9.9\toolbox\vl_setup'

% Code dependencies
addpath('.\..\util\')

[gtIdx,gtFile] = parseGT('GYRBO',gt_dir);
[candIdx,candFile] = parseCand(cand_dir);
nq = length(candIdx);
nc = length(candFile);

tic
% 10 is hard coded
scores = zeros(nq,2);
rIdx = [];
yIdx = [];
gIdx = [];
disp(' ')
for k=1:nq
    
    disp(['Recombination on query ',num2str(k)])
    
    querySift = getQuerySift(candIdx(k,1),qDir);
    [~,da] = vl_ubcread(strcat(qDir,querySift));
    
    topN = 1;
    cand = candFile( candIdx(k,2) : candIdx(k,3) );
    for j=1:length(cand)
        candSift = [cand{j},'sift.txt'];
        [~,db] = vl_ubcread(strcat(dbDir,candSift));
        scores(k,j) = length(vl_ubcmatch(da,db));
    end
    [~,perm] = sort(scores(k,:),2,'descend');
    winners = cand( perm(1:topN) );
    
    idx = find(candIdx(k,1)==gtIdx(:,1));
    gt = gtFile( gtIdx(idx,2) : gtIdx(idx,3) );
    if textMatch(winners,gt)
        gIdx = [gIdx,candIdx(k)];
    elseif textMatch(cand,gt)
        yIdx = [yIdx,candIdx(k)];
    else
        rIdx = [rIdx,candIdx(k)];
    end
    
end
disp(' ')
ng = length(gIdx);
ny = length(yIdx);
nr = length(rIdx);
disp(['Successful recombination for ',num2str(ng),' queries out of ',...
    num2str(nq),' (',num2str(dRound(100*ng/nq,0)),'%).'])
disp(['Filter lost ', num2str(ny),' queries out of ',...
    num2str(ng+ny),' (',num2str(dRound(100*ny/(ny+ng),0)),'%).'])
elapsed_time = toc;
disp(['Average match time of ',dRound(num2str(elapsed_time/nc),-1),' seconds.'])