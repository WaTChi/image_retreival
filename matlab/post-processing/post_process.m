function [results] = post_process(method,reset)

% [results] = post_process(method)
% 
% First coded 8 Jan 2011 by Aaron Hallquist
% Latest revision 28 Jan 2011 by Aaron Hallquist
% 
% DESCRIPTION
%   This function computes the post-processing results on the input query
%   set with the given method parameters. The output is a results structure
%   which is saved on the hard drive.
% 
% INPUTS
%   method:     Structure detailing what method to use. Contents depend on
%               the values of certain required fields, listed below:
%       .set:           Integer indicating which query set to use
%       .cell_dist:     Maximum distance between location and cell center
%                       to search that cell. If equal to zero, we use no
%                       combination and only the nearest cell is searched
%       .decision:      Decision method. Current modes supported...
%           'linear-a.b':   Linear combination of votes and scores.
%                           - If 'linear' is chosen, the a.b parameter
%                             refers to the ratio score weight. The vote
%                             score weight will be 1, while the ratio score
%                             weight will be str2double('a.b')
%           'bayes-xyz':    Uses Naive Bayes classifier : must be trained
%                           - If 'bayes' is chosen, the 'xyz' parameter
%                             refers to which features will be used in the
%                             bayes decision.
%                               - 'd' for distance
%                               - 'v' for vote
%                               - 'r' for ratio score
%                           - e.g. 'dv' decides based on distance and vote
%                             while 'vs' is based on vote and ratio score
%       .distribution:  Noisy location distribution. Currently supports:
%           'exact':    Using the camera GPS coordinate only
%           'unif':     Generates noisy locations in a uniform distribution
%                       - Distribution parameter .param must be specified,
%                         indicating the maximum radius from GPS location
%           'expo':     Generates noisy locations, exponential distribution
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

% Fixed parameters
ncand = 100; % number of candidate images from vote only
nfilt = 25; % number of images to filter down to from candidates
ntop = 10; % maximum top n for results

% Set reset flag
if nargin < 2
    reset = 0;
end

% Fixed directory
dbDir = 'E:\Research\collected_images\earthmine-fa10.1,culled\37.871955,-122.270829\';
gtDir = 'E:\Research\app\code\matlab\ground-truth\';

% Adjusted directories based on inputs
qDir = ['E:\query',num2str(method.set),'\'];
vDir = ['E:\Research\results\query',num2str(method.set), ...
    '\matchescells(g=100,r=d=236.6),query',num2str(method.set), ...
    ',kdtree4,threshold=70k,searchparam=1024\'];

% Ratio score code initialization
run 'C:\vlfeat-0.9.9\toolbox\vl_setup'

% Code dependencies
addpath('.\..\util\')
addpath('.\..\bayes\')

% Read the decision parameter
idx = strfind(method.decision,'-');
decision = method.decision(1:idx-1);
decis_prm = method.decision(idx+1:end);

% Get a list of queries
query = dir(qDir);
query = strvcat(query(3:end).name);
query = unique(str2double(cellstr(query(:,5:8)))); % query numbers
query(isnan(query)) = [];
nq = length(query); % number of queries

% Load results structure
results_file = ['.\',decision,'\query',num2str(method.set), ...
                method.distribution, ...
                num2str(dRound(method.cell_dist,0)),'_results.mat'];
if reset
    results.run = cell(1,0);
    results.match = zeros(ntop,0);
    results.total = zeros(ntop,0);
    results.match_pct = zeros(ntop,0);
    results.query_pct = zeros(nq,0);
else
    try
        load(results_file)
    catch
        results.run = cell(1,0);
        results.match = zeros(ntop,0);
        results.total = zeros(ntop,0);
        results.match_pct = zeros(ntop,0);
        results.query_pct = zeros(nq,0);
    end
end

% Load query ratio scores
scores_file = ['.\scores\query',num2str(method.set),'_scores.mat'];
try
    load(scores_file)
catch
    query_scores.scores = cell(nq,1);
    query_scores.number = query;
end

% Load the classifier if necessary
if strcmp(decision,'bayes')
    bayes_file = ['.\bayes\classifier\',method.distribution, ...
                       num2str(dRound(method.cell_dist,0)),'_bayes.mat'];
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
cell_dist = method.cell_dist;
run = decis_prm;
match = zeros(ntop,1);
total = zeros(ntop,1);
match_pct = zeros(ntop,1);
query_pct = zeros(nq,1);

% Set run

fprintf('\nRunning post processing...\n')

% Iterate through each query
for k=[13,16,21,24,28,68]%1:nq

    fprintf(['\nProcessing query ',num2str(k),'... '])

    % Get ground truth matches
    gt = getGT(query(k),method.set,gtDir);

    % Get list of cells and cell locations
    [qCell,cLat,cLon] = getCells(query(k),vDir);

    % Get query sift file name and query location
    idx = strfind(qCell{1},',');
    query_name = qCell{1}(1:idx(3)-1);
    qLat = str2double(query_name(idx(1)+1:idx(2)-1));
    qLon = str2double(query_name(idx(2)+1:end-8));

    % Get fuzzy point locations
    if strcmp(method.distribution,'exact')
        lats = qLat; lons = qLon;
    elseif strcmp(method.distribution,'unif')
        if ~isfield(method,'param')
            rad = 75;
        else
            rad = method.param;
        end
        [lats,lons] = getFuzzyLocs(qLat,qLon,rad,rad/25);
    else % if strcmp(method.distribution,'expo')
        if ~isfield(method,'param')
            rad = 200;
        else
            rad = 4*method.param;
        end
        [lats,lons] = getFuzzyLocs(qLat,qLon,rad,rad/25);
    end

    % Create cell groupings from fuzzy points and iterate through them
    cellgroups = groupFuzzy(lats,lons,cLat,cLon,cell_dist);

    ma = 0;
    to = 0;
    for cg=cellgroups

        % Combine cells in this grouping
        comb_cells = qCell(cg.idx);
        [cand,cand_vote,cand_lat,cand_lon,nfeat] = cellCombine(comb_cells,vDir,ncand);

        % Iterate through each fuzzy point and post process
        for j=1:cg.npts
            
            % Get the candidate distances and filter down to nfilt
            if strcmp(decision,'linear')
                cand = cand(1:nfilt);
                cand_vote = cand_vote(1:nfilt);
            else
                cand_dist = latlonDistance(cg.lats(j),cg.lons(j),cand_lat,cand_lon);
                param = [cand_dist,cand_vote,nan(ncand,1)];
                [cand,param,~] = rankNcand(cand,param,method,nfilt);
                cand_dist = param(1,:);
                cand_vote = param(2,:);
            end
            
            % Get the candidate ratio scores | read from file if possible
            idx = find(query_scores.number==query(k),1,'first');
            img_scores = query_scores.scores{idx};
            if ~iscell(img_scores)
                img_scores = cell(0,2);
            end
            [cand_score,img_scores] = getRatioScores(...
                cand,img_scores,query_name,qDir,dbDir);
            query_scores.scores{idx} = img_scores;
            cand_score = cand_score / nfeat;

            % Post process the features parameters to get a final ranking
            if strcmp(decision,'linear')
                param = [cand_vote,cand_score];
            else % if strcmp(decision,'bayes')
                param = [cand_dist,cand_vote,cand_score];
            end
            cand = rankNcand(cand,param,method,ntop);
            
            % Check for matches
            m = zeros(ntop,1);
            cand_idx = 1;
            while cand_idx<=ntop && ~textMatch(cand(cand_idx),gt)
                cand_idx = cand_idx+1;
            end
            m(cand_idx:end)=1;
            
            % weigh result if exponential distribution
            if strcmp(method.distribution,'expo')
                d = latlonDistance(qLat,qLon,cg.lats(j),cg.lons(j));
                if ~isfield(method,'param')
                    param = 50;
                else
                    param = method.param;
                end
                weight = exp(-d/param);
            else % exact or uniform
                weight = 1;
            end
                
            match = match + weight*m;
            total = total + weight;

            ma = ma + weight*m(1);
            to = to + weight;

        end

    end

    query_pct(k) = ma / to;
    
    % store query ratio scores
    save(scores_file,'query_scores')

end

% store results
results.run(1,end+1) = decis_prm;
results.match(:,end+1) = match;
results.total(:,end+1) = total;
results.match_pct(:,end+1) = match./total;
results.query_pct(:,end+1) = query_pct;

save(results_file,'results')

fprintf('\nRun complete.\n\n')