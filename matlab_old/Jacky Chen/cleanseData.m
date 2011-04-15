function [] = cleanseData(fileName)

%   cleanseData(fileName) 
%
%   Input:
%       - fileName: name of .mat file that stores observations from the
%       output of convertImagesToObservations
%
%   Output:
%       - Nothing is returned. The .mat file is changed and resides in the
%       same folder.
%

if ~exist(fileName)         
    error('File not found');
else
    load(fileName);
    N = size(observations, 1);
end

%We remove features that appear for every image, no images, and only one
%image, because they are not particularly helpful in cooccurrences.
remove = [N 0 1];

for i = 1:size(remove, 2)
    s = sum(observations);
    ind = find(s == remove(i));
    [m, n] = size(ind);
    
    while(n ~= 0)
        if mod(n, 100) == 1
            fprintf('%d columns to remove\n', n-1);
        end
        observations(:, ind(n)) = [];
        n = n-1;
    end
end

save(fileName, 'observations');