function [classifier] = trainQueryClassifier(method,newtrain)

% [classifier] = trainQueryClassifier(method,ndist,nvote)
% 
% First coded 31 Jan 2011 by Aaron Hallquist
% Latest revision 10 Feb 2011 by Aaron Hallquist
% 
% DESCRIPTION
%   This function trains the distributions of an image search query with 
%   the method specified in order to determine the best match.
% 
% INPUTS
%   method:     Structure detailing what query method to use; Uses the same
%               format as in the function post_process
%   nsamps:     Number of training samples per query
%   newtrain:   Optional flag (1 indicates classifier reset; default = 0)
% 
% OUTPUTS
%   classifier: Trained classifier (saved out as well)

% Set reset flag
if nargin < 2
    newtrain = 0;
end

% Fixed directory
gtDir = 'C:\matlab_local\ground-truth\';

% Adjusted directories based on inputs
qDir = ['Z:\query',num2str(method.set),'\'];
vDir = ['C:\matlab_local\results\query',num2str(method.set),'\'];

% Get list of cells and cell locations
[~,cLat,cLon] = getCells;

% Code dependencies
addpath('.\..\util\')
addpath('.\..\bayes\')

% Initialize size parameters
C = 2; % number of classes
M = 2; % number of features

% Parse distribution parameter
idx = strfind(method.distribution,'-');
distr = method.distribution(1:idx-1);
distr_prm = str2double(method.distribution(idx+1:end));

% Load the current classifier 
bayes_file = ['.\bayes\classifier\',distr, ...
              num2str(dRound(method.cell_dist,0)),'_bayes.mat'];
if newtrain
    classifier.spacing = [25,.01];
else
    try
        load(bayes_file)
    catch
        classifier.spacing = [25,.01];
    end
end

% Get a list of queries
query = struct2cell(dir(qDir));
query = query(1,:)';
sift_idx = ~cellfun('isempty',strfind(query,'sift.txt'));
hdf5_idx = ~cellfun('isempty',strfind(query,'.hdf5'));
query_name = query( sift_idx & ~hdf5_idx );
query = strvcat(query_name); % query numbers
query = str2double(cellstr(query(:,5:8)));
nq = length(query); % number of queries

% Initialize variables
if ~isfield(method,'cell_dist')
    cell_dist = 336.6;
else
    cell_dist = method.cell_dist;
end
div = 25; % used to determine sample spacing; nsamps ~ pi * div^2

fprintf('\nRunning classifier training...\n')

for k=1:nq
    
    fprintf(['\nProcessing query ',num2str(k),'... '])
    
    % Get the number of features
    featid = fopen([qDir,query_name{k}]);
    nfeat = fscanf(featid,'%d',1);
    fclose(featid);
    
    % Get ground truth matches
    gt = getGT(query(k),method.set,gtDir);
    
    % Get query location
    idx = strfind(query_name{k},',');
    qLat = str2double(query_name{k}(idx(1)+1:idx(2)-1));
    qLon = str2double(query_name{k}(idx(2)+1:end-8));
    
    % Get fuzzy point locations
    if strcmp(distr,'exact')
        lats = qLat; lons = qLon;
    elseif strcmp(distr,'unif')
        rad = distr_prm;
        [lats,lons] = getFuzzyLocs(qLat,qLon,rad,rad/25);
    else % if strcmp(distr,'expo')
        rad = 4*distr_prm;
        [lats,lons] = getFuzzyLocs(qLat,qLon,rad,rad/25);
    end
    
    % Create cell groupings from fuzzy points and iterate through them
    cellgroups = groupFuzzy(lats,lons,cLat,cLon,cell_dist);
    
    for cg=cellgroups
        
        % Get cell combination results, files, and candidate locations
        [cand,cand_vote,cand_lat,cand_lon] = getCand(cg.idx,query(k),vDir);
        cand_vote = cand_vote / nfeat;
        ncand = length(cand);

        % Iterate through each fuzzy point and gather features
        for j=1:cg.npts
            
            % Get the candidate distances
            cand_dist = latlonDistance(cg.lats(j),cg.lons(j),cand_lat,cand_lon);
            
            % Get the classes for these candidates
            bayes_classes = 2*ones(ncand,1); % Class 2 = no match
            for c = 1:ncand
                if ~isempty(find(strcmp(cand{c},gt),1))
                    bayes_classes(c) = 1; % Class 1 = match
                end
            end   
            
            % Update classifier
            weights = ones(ncand,1);
            if strcmp(distr,'expo')
                d = latlonDistance(qLat,qLon,cg.lats(j),cg.lons(j));
                weights(:) = 1 / (2*pi*(div/4)^2) * exp(-d/distr_prm);
            elseif strcmp(distr,'unif')
                weights(:) = 1 / (pi*div^2);
            end
            bayes_features = [cand_dist,cand_vote];
            classifier = trainbayes(bayes_features,bayes_classes,classifier,weights);
            
        end

    end
    
end

% store classifier
save(bayes_file,'classifier')

fprintf('\nRun complete.\n\n')