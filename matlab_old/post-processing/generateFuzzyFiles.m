% Generates fuzzy locations for each query and saves them out to a file

addpath('.\..\util\')
addpath('.\..\bayes\')

% Directories
qDir = 'Z:\query4-matlab\';
fDir = '.\tmp\';

% Get a list of queries
query = struct2cell(dir(qDir));
query = query(1,:)';
sift_idx = ~cellfun('isempty',strfind(query,'sift.txt'));
hdf5_idx = ~cellfun('isempty',strfind(query,'.hdf5'));
dsc_idx = ~cellfun('isempty',strfind(query,'DSC_'));
query_name = query( sift_idx & ~hdf5_idx & dsc_idx );
query = strvcat(query_name); % query numbers
query = str2double(cellstr(query(:,5:8)));
nq = length(query); % number of queries

for k=1:nq

    % Get query location
    idx = strfind(query_name{k},',');
    qLat = str2double(query_name{k}(idx(1)+1:idx(2)-1));
    qLon = str2double(query_name{k}(idx(2)+1:end-8));

    % Get fuzzy point locations
    [lats,lons] = getFuzzyLocs(qLat,qLon,200,8);
    
    % Write to file
    fid = fopen([fDir,query_name{k}(1:end-8),'.fuz'],'w');
    for j=1:length(lats)
        fprintf(fid,'%1.15f\t%1.14f\n',lats(j),lons(j));
    end
    fclose(fid);
    
end