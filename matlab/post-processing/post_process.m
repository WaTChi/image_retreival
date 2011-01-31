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
%           'linear':   Linear combination of votes and scores.
%                       - If 'linear' is chosen, the field .ws must be
%                         specified, indicating the weight on the ratio
%                         score; the vote weight is set to 1
%           'bayes':    Uses Naive Bayes classifier : must be trained
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

% Set reset flag
if nargin < 2
    reset = 0;
end

% Fixed directory
dbDir = 'E:\Research\collected_images\earthmine-new,culled\37.871955,-122.270829\';
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

% Load results structure
results_file = ['.\',method.decision,'\query',num2str(method.set), ...
                method.distribution, ...
                num2str(dRound(method.cell_dist,0)),'_results.mat'];
if reset
    results = struct;
else
    try
        load(results_file)
    catch
        results = struct;
    end
end

% Get a list of queries
query = dir(qDir);
query = strvcat(query(3:end).name);
query = unique(str2double(cellstr(query(:,5:8)))); % query numbers
query(isnan(query)) = [];
nq = length(query); % number of queries

% Load query ratio scores
scores_file = ['.\scores\query',num2str(method.set),'_scores.mat'];
try
    load(scores_file)
catch
    query_scores.scores = cell(nq,1);
    query_scores.number = query;
end

% Load the classifier if necessary
if strcmp(method.decision,'bayes')
    bayes_file = ['.\bayes\classifier\',method.distribution, ...
                       num2str(dRound(method.cell_dist,0)),'_bayes.mat'];
    try
        load(bayes_file)
    catch
        train_method = method;
        train_method.set = 2;
        fprintf('\nMust train classifier. Doing this now...\n')
        bayes = trainQueryClassifier(train_method,100);
    end
end

               
               
% ---------------------------------------------
% Post-processing | Main part
% ---------------------------------------------
    
% Initialize variables
ncand = 50; % number of candidate images from vote only
nfilt = 20; % number of images to filter down to from candidates
ntop = 10; % maximum top n for results
cell_dist = method.cell_dist;
if strcmp(method.decision,'linear')
    ws = method.ws;
    nruns = length(ws);
    wv = ones(nruns,1);
else % if strcmp(method.decision,'bayes')
    nruns = 1;
end
match = zeros(ntop,nruns);
total = zeros(ntop,nruns);
match_pct = zeros(ntop,nruns);
query_pct = zeros(nq,nruns);

fprintf('\nRunning post processing...\n')

% Iterate through each query
for k=1:nq

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
        rad = method.param;
        [lats,lons] = getFuzzyLocs(qLat,qLon,rad,rad/50);
    else % if strcmp(method.distribution,'expo')
        rad = 4*method.param;
        [lats,lons] = getFuzzyLocs(qLat,qLon,rad,rad/50);
    end

    % Create cell groupings from fuzzy points and iterate through them
    cellgroups = groupFuzzy(lats,lons,cLat,cLon,cell_dist);

    ma = 0;
    to = 0;
    for cg=cellgroups

        % Combine cells in this grouping
        comb_cells = qCell(cg.idx);
        [cand,cand_vote,cand_lat,cand_lon,nfeat] = cellCombine(comb_cells,ncand,vDir);

        % Iterate through each fuzzy point and post process
        for j=1:cg.npts
            
            % Get the candidate distances
            cand_dist = latlonDistance(cg.lats(j),cg.lons(j),cand_lat,cand_lon);
            
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

            % Evaluate the post processing
            if strcmp(method.decision,'linear')
                [m,p] = evaluateScores(cand_vote,cand_score,cand, ...
                    gt,wv,ws,method.rerank);
                % weigh result if exponential distribution
                if strcmp(method.distribution,'exponential')
                    d = latlonDistance(qLat,qLon,cg.lats(j),cg.lons(j));
                    weight = exp(-d/50);
                else % exact or uniform
                    weight = 1;
                end
            else % if strcmp(method.decision,'bayes')
                bayes_features = [cand_dist,cand_vote,cand_score];
                [~,condP] = classify(bayes_features,bayes);
                [condP,prob_idx] = sort(condP(:,1),'descend');
                cand = cand(prob_idx);
            end
                
            match = match + weight*m;
            poss = poss + weight*p;
            total = total + weight;

            ma = ma + weight*m(1);
            to = to + weight;

        end

    end

    mismatch_pct = ma / to;
    fprintf(fid,[query_name(1:end-8),'\t%1.4f\n'],mismatch_pct);

end

% store query ratio scores
save(scores_file,'query_scores')
% fclose(fid);

% store results
results.match(end-nruns+1:end,:) = match;
results.poss(end-nruns+1:end,:) = poss;
results.total(end-nruns+1:end,:) = total;
save(results_file,'results')

fprintf('\nRun complete.\n\n')