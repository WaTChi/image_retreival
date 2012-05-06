function query = getQuery(qDir)

% Grabs the queries from the query directory folder

query = struct2cell(dir(qDir));
query = query(1,:)';
sift_idx = ~cellfun('isempty',strfind(query,'sift.txt'));
hdf5_idx = ~cellfun('isempty',strfind(query,'.hdf5'));
query = query( sift_idx & ~hdf5_idx );