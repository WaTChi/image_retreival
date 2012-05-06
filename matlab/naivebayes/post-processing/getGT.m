function [gt] = getGT(query,query_set,gt_dir)
% [gt] = getGT(query,query_set,gt_dir)
% 
% DESCRIPTION
%   This function assembles a list of the ground truth matches for the
%   given query and query set. It adds 'sift.txt' to the ends to match up
%   with the cell query format.

% Remove 'sift.txt' from query name
query = query(1:end-8);

% Ground truth filename format
fn = [gt_dir,'gt_',query_set,'.txt'];
textdata = textread(fn,'%s');
n = length(textdata);
idx = find(~cellfun('isempty',strfind(textdata,query)),1);
if isempty(idx)
    error('Could not find query in ground truth file.')
end

gt = cell(0,1);
stop = false;
dbFormat = {'-0002';'-0003';'-0004';'-0008';'-0009';'-0010'};
while ~stop
    idx = idx+1;
    if idx>n
        stop = true;
    else
        file = textdata{idx};
        suffix = file(end-4:end);
        if any(strcmp(suffix,dbFormat))
            gt{end+1,1} = [file,'sift.txt'];
        else
            stop = true;
        end
    end
end