function [qNum,qName] = parseQuery(vote_dir)

query = dir(vote_dir);
idx = [];
for k=1:length(query)
    if ~isempty(strfind(query(k).name,'.res'))
        idx = [idx,k];
    end
end
query = query(idx);
qName = strvcat(query.name);
qNum  = str2num(qName(:,5:8));