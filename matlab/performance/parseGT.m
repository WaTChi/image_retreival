function [gtIdx,gtFile] = parseGT(divs,root_dir)

if nargin<2
    root_dir = 'E:\Research\app\code\';
end
fnPrefix = 'groundtruth';
fnSuffix = 'Matlab.txt';
nq = length(importdata([root_dir,fnPrefix,divs(1),fnSuffix]));
nd = length(divs);
gt = cell(nq,nd);
for j=1:nd
    gt(:,j) = importdata([root_dir,fnPrefix,divs(j),fnSuffix]);
end

gtIdx = zeros(nq,3);
gtFile = cell(10000,1);
idx = 1;
for k=1:nq
    gtIdx(k,1) = str2num(gt{k,1}(5:8));
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