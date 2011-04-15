function [results] = post_process(method,reset)

% [results] = post_process(method)
% 
% First coded 8 Jan 2011 by Aaron Hallquist
% Latest revision 3 Apr 2011 by Aaron Hallquist
% 
% DESCRIPTION
%   This function computes the post-processing results on the input query
%   set with the given method parameters. The output is a results structure
%   which is saved on the hard drive.
% 
% INPUTS
%   method:     Structure detailing what method to use. Contents depend on
%               the values of certain required fields, listed below:
%       .set:           String determing which set to use
%       .cell_dist:     Maximum distance between location and cell center
%                       to search that cell. If equal to zero, we use no
%                       combination and only the nearest cell is searched.
%                       This is now optional with default = 336.6
%       .decision:      Decision method. Current modes supported...
%           'vote-':        Uses the vote total only.
%           'bayes-xyz':    Uses 2D Bayes classifier : must be trained
%                           - If 'bayes' is chosen, the 'xyz' parameter
%                             refers to which features will be used in the
%                             bayes decision.
%                               - 'd' for distance
%                               - 'v' for vote
%                           - e.g. 'dv' decides based on distance and vote
%                             while 'v' is based on vote only
%       .canddiv:       Query number of candidate divisions for classifier.
%                       - e.g. canddiv = [10 100] splits queries up into
%                         those with fewer than 10 candidates, those with
%                         more than 100, and those between 10 and 100.
%       .distribution:  Noisy location distribution. Currently supports:
%           'exact-':   Using the camera GPS coordinate only
%           'unif-x':   Generates noisy locations in a uniform distribution
%                       - Distribution parameter .param must be specified,
%                         indicating the maximum radius from GPS location
%           'expo-x':   Generates noisy locations, exponential distribution
%                       - Distribution parameter .param must be specified,
%                         indicating the mean distance from GPS location
%                       - Samples cap at 4 times mean from GPS location
%                           (e.g. caps at 200m if mean is 50m)
%   reset:      Optional flag (1 indicates results reset; default = 0)
% 
% OUTPUTS
%   results:    Structure containing the details of post-processing
%               - type 'help ppstruct' for information on the format of the
%                 results structure (depends on the contents of method)

% Code dependencies
addpath('.\..\util\')
addpath('.\..\bayes\')

% Set reset flag
if nargin < 2
    reset = 0;
end

% Fixed parameters
ntop = 10; % maximum top n for results
div = 25; % used to determine sample spacing; nsamps ~ pi * div^2

% Fixed directory
gtDir = 'C:\matlab_local\ground-truth\';

% Adjusted directories based on inputs
qDir = ['Z:\',method.set,'\'];
vDir = ['C:\matlab_local\results\',method.set,'\'];

% Get list of cells and cell locations
[~,cLat,cLon] = getCells;

% Parse the decision parameter
idx = strfind(method.decision,'-');
decision = method.decision(1:idx-1);
decis_prm = method.decision(idx+1:end);

% Parse distribution parameter
idx = strfind(method.distribution,'-');
distr = method.distribution(1:idx-1);
distr_prm = str2double(method.distribution(idx+1:end));

% Flag to record top n in files
recordFlag = strcmp(distr,'exact') && strcmp(decis_prm,'dv');
if recordFlag
    time = datestr(now);
    idx = strfind(time,':');
    time = time([1:idx(1)-1,idx(1)+1:idx(2)-1]);
    rdir = ['.\exact-runs\',method.set,'\',time,'\'];
    mkdir(rdir);
end

% Get a list of queries
query = getQuery(qDir);
nq = length(query); % number of queries

% Load results structure
results_file = ['.\',decision,'\',method.set,'_',distr, ...
                num2str(dRound(method.cell_dist,0)),'results.mat'];
if reset
    results.run = cell(1,0);
    results.match = zeros(ntop,0);
    results.total = zeros(1,0);
    results.match_pct = zeros(ntop,0);
    results.query_pct = zeros(nq,0);
else
    try
        load(results_file)
    catch
        results.run = cell(1,0);
        results.match = zeros(ntop,0);
        results.total = zeros(1,0);
        results.match_pct = zeros(ntop,0);
        results.query_pct = zeros(nq,0);
    end
end

% Load the classifier if necessary
if strcmp(decision,'bayes')
    bayes_file = ['.\',decision,'\classifier\',distr, ...
        num2str(dRound(method.cell_dist,0)),'_',decision,'.mat'];
    try
        load(bayes_file)
    catch
        error('No Bayes classifier trained.')
    end
end


% ---------------------------------------------
% Post-processing | Main part
% ---------------------------------------------
    
% Initialize variables
if ~isfield(method,'cell_dist')
    cell_dist = 336.6;
else
    cell_dist = method.cell_dist;
end
total = 0;
match = zeros(ntop,1);
query_pct = zeros(nq,1);

fprintf('\nRunning post processing...\n')

% Iterate through each query
for k=1:nq
    
    fprintf(['\nProcessing query ',num2str(k),'... '])
    
    % Get the number of features
    featid = fopen([qDir,query{k}]);
    nfeat = fscanf(featid,'%d',1);
    fclose(featid);

    % Get ground truth matches
    gt = getGT(query{k},method.set,gtDir);

    % Get query location
    [qLat,qLon] = getQueryLocation(query{k},method.set);

    % Get fuzzy point locations
    if strcmp(distr,'exact')
        lats = qLat; lons = qLon;
    elseif strcmp(distr,'unif')
        rad = distr_prm;
        [lats,lons] = getFuzzyLocs(qLat,qLon,rad,rad/div);
    else % if strcmp(dist,'expo')
        rad = 4*distr_prm;
        [lats,lons] = getFuzzyLocs(qLat,qLon,rad,rad/div);
    end

    % Create cell groupings from fuzzy points and iterate through them
    cellgroups = groupFuzzy(lats,lons,cLat,cLon,cell_dist);

    ma = 0;
    to = 0;
    for cg=cellgroups
        
        % Get cell combination results, files, and candidate locations
        [cand,cand_vote,cand_lat,cand_lon] = getCand(cg.idx,query{k},vDir);
        ncand = length(cand);
        cand_vote = cand_vote / nfeat;

        % Iterate through each fuzzy point and post process
        for j=1:cg.npts
            
            % Get the candidate distances
            cand_dist = latlonDistance(cg.lats(j),cg.lons(j),cand_lat,cand_lon);

            % Post process the features parameters to get a final ranking
            param = [cand_dist,cand_vote];
            [cand_nres,p,s] = rankNcand(cand,param,method,ncand,ntop);
            
            % Check for matches
            m = zeros(ntop,1);
            nres = min(ntop,length(cand_nres));
            cand_idx = 1;
            while cand_idx<=nres && ~textMatch(cand_nres(cand_idx),gt)
                cand_idx = cand_idx+1;
            end
            if cand_idx <= nres
                m(cand_idx:end)=1;
            end
            
            % weigh result if exponential distribution
            if strcmp(distr,'expo')
                d = latlonDistance(qLat,qLon,cg.lats(j),cg.lons(j));
                weight = 1 / (2*pi*(div/4)^2) * exp(-d/distr_prm);
            elseif strcmp(distr,'unif')
                weight = 1 / (pi*div^2);
            else % strcmp(distr,'exact')
                weight = 1;
            end
            
            % Update match information
            match = match + weight*m;
            ma = ma + weight*m(1);
            total = total + weight;
            to = to + weight;
            
            % Write to file if conditions are correct
            if recordFlag
                fn = [rdir,query{k}(1:end-8),'.bay'];
                fid = fopen(fn,'w');
                for res_i = 1:nres
                    res_p = s(res_i);
                    res_v = p(res_i,2);
                    res_d = p(res_i,1);
                    res_name = cand_nres{res_i}(1:end-8);
                    fprintf(fid,'%.4f\t%d\t%d\t%s\n', ...
                        res_p,round(nfeat*res_v),round(res_d),res_name);
                end
                fclose(fid);
            end

        end

    end

    query_pct(k) = ma / to;

end

% store results
results.run{1,end+1} = [decis_prm,'-',mat2str(method.canddiv)];
results.total(1,end+1) = total;
results.match(:,end+1) = match;
results.match_pct(:,end+1) = match / total;
results.query_pct(:,end+1) = query_pct;

save(results_file,'results')

fprintf('\nRun complete.\n\n')