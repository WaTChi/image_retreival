function filename = getQuerySift(qNum,qDir)

query = dir(qDir);
idx = 0;
for k=1:length(query)
    if ~isempty(strfind(query(k).name,num2str(qNum))) && ...
        ~isempty(strfind(query(k).name,'sift.txt')) && ...
        isempty(strfind(query(k).name,'.hdf5'))
        idx = k;
    end
end
filename = query(idx).name;