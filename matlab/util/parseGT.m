function [gtIdx,gtFile] = parseGT(divs,root_dir)

% [gtIdx,gtFile] = parseGT(divs,root_dir)
% 
% First coded 9 Dec 2010 by Aaron Hallquist
% Latest revision 9 Dec 2010 by Aaron Hallquist
% 
% K = number of queries
% 
% DESCRIPTION:
%   This file parses the ground truth text files for use in Matlab. For
%   each specified query, this function retrieves the query number and a 
%   list of database image file names which appear in the specified ground
%   truth categories.
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
%               refers to the database matches for query k.
%               Examples: 
%                   gtIdx(k,:) = [ 8847 13 12 ]
%                       Query #8847 has no matches in the divisions
%                       specified [ gtFile( 13:12 ) is empty ]
%                   gtIdx(j,:) = [ 8860 65 68 ]
%                       Query #8860 has 4 matches in the divisions
%                       specified [ the entries of gtFile( 65:68 ) ]

% Default directory
if nargin<2
    root_dir = 'E:\Research\app\code\matlab\ground-truth\';
end

% Ground truth filename format
fnPrefix = 'groundtruth';
fnSuffix = '.txt';

% Number of queries and ground truth divisions
nq = length(importdata([root_dir,fnPrefix,divs(1),fnSuffix]))
nd = length(divs);

% Fetch ground truth data
gt = cell(nq,nd);
for j=1:nd
    gt(:,j) = importdata([root_dir,fnPrefix,divs(j),fnSuffix]);
end

% Parse the ground truth data into the output format
gtIdx = zeros(nq,3);
gtFile = cell(10000,1);
idx = 1;
for k=1:nq
    gtIdx(k,1) = str2double(gt{k,1}(5:8));
    gtIdx(k,2) = idx;
    for j=1:nd
        tmp = strfind(gt{k,j},'''');
        tmp = reshape(tmp,[2,length(tmp)/2]);
        for i=tmp
            gtFile{idx} = gt{k,j}(i(1)+1:i(2)-1);
            idx=idx+1;
        end
    end
    gtIdx(k,3) = idx-1;
end
gtFile = gtFile(1:idx-1);
gtIdx = sortrows(gtIdx);

end