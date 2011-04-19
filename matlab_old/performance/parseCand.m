function [candIdx,candFile] = parseCand(cand_dir)

qCand = dir(cand_dir);
idx = [];
for k=1:length(qCand)
    if ~isempty(strfind(qCand(k).name,'.top'))
        idx = [idx,k];
    end
end
qCand = qCand(idx);
nq = length(qCand);

candIdx = zeros(nq,3);
candFile = cell(10000,1);
idx = 1;
for k=1:nq
    candIdx(k,1) = str2num(qCand(k).name(5:8));
    candIdx(k,2) = idx;
    cand = textread([cand_dir,qCand(k).name],'%s');
    nc = length(cand);
    candFile(idx:idx+nc-1) = cand;
    idx = idx+nc;
    candIdx(k,3) = idx-1;
end
candFile = candFile(1:idx-1);
candIdx = sortrows(candIdx);

end