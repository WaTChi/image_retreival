function [scores,img_scores] = getRatioScores(...
    cand,img_scores,qSift,qDir,dbDir)

% [scores] = getRatioScores(cand,query,query_set)
% 
% DESCRIPTION
%   This function obtains the ratio scores for the input candidate images.
%   If the scores already exist, they are drawn. If not, then they are 
%   generated and stored.

ncand = length(cand);
scores = zeros(ncand,1);
for k=1:ncand
    mask = strcmp(img_scores(:,1),cand{k});
    if sum(mask)==0 % generate score
        fprintf('|')
        [~,da] = vl_ubcread(strcat(qDir,qSift));
        [~,db] = vl_ubcread(strcat(dbDir,cand{k}));
        scores(k) = length(vl_ubcmatch(da,db));
        img_scores{end+1,1} = cand{k};
        img_scores{end,2} = scores(k);
    else
        scores(k) = img_scores{mask,2};
    end
end