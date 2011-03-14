qDir = 'Z:\query4-matlab\';

% Get a list of queries
query = struct2cell(dir(qDir));
query = query(1,:)';
sift_idx = ~cellfun('isempty',strfind(query,'sift.txt'));
hdf5_idx = ~cellfun('isempty',strfind(query,'.hdf5'));
query_name = query( sift_idx & ~hdf5_idx );
query = strvcat(query_name); % query numbers
query = str2double(cellstr(query(:,5:8)));
nq = length(query); % number of queries

nfeat = zeros(nq,1);
for k=1:nq
    featid = fopen([qDir,query_name{k}]);
    nfeat(k) = fscanf(featid,'%d',1);
    fclose(featid);
end