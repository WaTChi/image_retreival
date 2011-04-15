function [gt] = getGT(query,query_set,gt_dir,f139)

% TMP F139 ADDITION FOR BUG FIX

% [gt] = getGT(query,query_set,gt_dir)
% 
% DESCRIPTION
%   This function assembles a list of the ground truth matches for the
%   given query and query set. It adds 'sift.txt' to the ends to match up
%   with the cell query format.

if query_set == 3
    [gtIdx,gtFile] = parseGT('GYRBO',[gt_dir,'query3\']);
else
    [gtIdx,gtFile] = parseGT('A',[gt_dir,'query',num2str(query_set),'\']);
end

if query == 139 && f139
    idx = find(gtIdx==query,2,'first');
    idx = idx(2);
else
    idx = find(gtIdx==query,1,'first');
end
gt = gtFile(gtIdx(idx,2):gtIdx(idx,3));
for k=1:length(gt)
    gt{k} = [gt{k},'sift.txt'];
end