function [] = binaryFilter3(dataset, varargin)

%   binaryFilter3(dataset, varargin) 
%
%   Input:
%       - dataset: string that specifies the directory and the file names.
%       Typically of the form "date_dataset"
%       - A .txt file that contains results from running fabmap
%       Typically of the form "result'dataset'.txt"
%
%   Output:
%       - A .txt file that indicates correct image pairs by image names and
%       lists them in descending counts order (unsorted refers to names)
%

% number of candidate pairs to be looked at; default is 1000 (which is 
% essentially all pairs)
if size(varargin, 2)
    top = varargin{1};
else
    top = 1000;
end

dir = ['.\', dataset, '\']; % must end in file separator (e.g. '\'); this 
            % is an input directory that must consist of *.mat files 
            % specifying SIFT features according to a format that Matt 
            % Carlberg was using

siftFeatures = ls(fullfile(dir, '*.feat')); 

threshold = 0.09;

fileName = strcat('result', dataset);
fileName = strcat(fileName, '.txt');

fid = fopen(fileName, 'r');
textscan(fid,'%*[^\n]', 1);
match = textscan(fid, '%d %d %*[^\n]', 'delimiter', '-');
i = match{1};
j = match{2};
top = min(top, size(i, 1));

candList = []; % stores image pairs
queryList = [];
fprintf('Examining the top %d candidates for loop closure...\n', top);
for k = 1:top
    a = i(k);
    b = j(k);
    [score1, Q1] = nnRatio(i(k), j(k), dataset);
    [score2, Q2] = nnRatio(j(k), i(k), dataset);
    if score1 > threshold*Q1 || score2 > threshold*Q2
        fprintf('Images %d and %d make a genuine loop closure.\n', a, b);
        candList = [candList, a];
        queryList = [queryList, b];
    else
        fprintf('Images %d and %d are not guaranteed to make a genuine loop closure.\n', a, b);
    end
end
fclose(fid);

reportName = strcat('unsortedMatches', dataset);
reportName = strcat(reportName, '.txt');
report = fopen(reportName, 'wt');
C = size(candList, 2);
for x = 1:C
    u = candList(x);
    v = queryList(x);
    fprintf(report,'%s-%s', siftFeatures(u,:), siftFeatures(v,:));
    fprintf(report,'\n');
end
fclose(report);
fprintf('Complete.\n');