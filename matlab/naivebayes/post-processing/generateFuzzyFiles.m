% Generates fuzzy locations for each query and saves them out to a file

addpath('.\..\util\')
addpath('.\..\bayes\')

% Directories
set = 4;
qDir = 'Z:\query4\';
fDir = '.\tmp\';

% Get a list of queries
query = struct2cell(dir(qDir));
query = query(1,:)';
sift_idx = ~cellfun('isempty',strfind(query,'sift.txt'));
hdf5_idx = ~cellfun('isempty',strfind(query,'.hdf5'));
query = query( sift_idx & ~hdf5_idx );
nq = length(query); % number of queries

for k=1:nq

    % Get query location
    [qLat,qLon] = getQueryLocation(query{k},set);

    % Get fuzzy point locations
    [lats,lons] = getFuzzyLocs(qLat,qLon,200,8);
    
    % Write to file
    fid = fopen([fDir,query{k}(1:end-8),'.fuz'],'w');
    for j=1:length(lats)
        fprintf(fid,'%1.15f\t%1.14f\n',lats(j),lons(j));
    end
    fclose(fid);
    
end