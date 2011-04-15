dump = textread('dump.txt','%s');

nq=112;
rp=30;
match = dump(1:rp:end);
match = cellfun('isempty',strfind(match,'False'));
query = strvcat(dump(2:rp:end));
query = str2double(cellstr(query(:,5:8)));
[query,idx] = sort(query,'ascend');
match = match(idx);

sum(match)/length(match)