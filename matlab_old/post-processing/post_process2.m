function [results] = post_process2(query_set,method)

% [results] = post_process(query_set,method)
% 
% First coded 8 Jan 2011 by Aaron Hallquist
% Latest revision 18 Jan 2011 by Aaron Hallquist
% 
% DESCRIPTION
%   This function computes the post-processing results on the input query
%   set with the given method parameters. The output is a results structure
%   which is saved on the hard drive.
% 
% INPUTS
%   query_set:  Integer determining which query set to post process.
%   method:     Structure detailing what method to use. Contents depend on
%               the values of certain required fields, listed below:
%       .fuzzy:         Boolean to use fuzzy reported locations
%       .cell_dist:     Maximum distance between location and cell center
%                       to search that cell. If equal to zero, we use no
%                       combination and only the nearest cell is searched
%       .decision:      Decision method. Current modes supported...
%           'linear':   Linear combination of votes and scores.
%       .rerank:        Boolean to rerank 10 and analyze
%       .filter:        Boolean to filter nearby db images or not.
%       .distribution:  Fuzzy point distribution: 'uniform' or 'exponential'
%       - type 'help ppstruct' for additional fields used based on the
%         values of the above fields
% 
% OUTPUTS
%   results:    Structure containing the details of post-processing
%               - type 'help ppstruct' for information on the format of the
%                 results structure (depends on the contents of method)

% Fixed directory
dbDir = 'E:\Research\collected_images\earthmine-new,culled\37.871955,-122.270829\';
gtDir = 'E:\Research\app\code\matlab\ground-truth\';

% Adjusted directories based on inputs
qDir = ['E:\query',num2str(query_set),'\'];
vDir = ['E:\Research\results(query',num2str(query_set), ...
    ')\matchescells(g=100,r=d=236.6),query',num2str(query_set), ...
    ',kdtree4,threshold=70k,searchparam=1024\'];

% Post-processing code initialization
run 'C:\vlfeat-0.9.9\toolbox\vl_setup'

% Code dependencies
addpath('.\..\util\')

% Load results structure
if strcmp(method.distribution,'exponential') && method.filter
    results_file = ['.\',method.decision,'\query',num2str(query_set), ...
        'expo',num2str(dRound(method.cell_dist,0)),'filt_results.mat'];
elseif strcmp(method.distribution,'exponential') && ~method.filter
    results_file = ['.\',method.decision,'\query',num2str(query_set), ...
        'expo',num2str(dRound(method.cell_dist,0)),'nofilt_results.mat'];
elseif method.fuzzy && ~method.rerank && method.filter
    results_file = ['.\',method.decision,'\query',num2str(query_set), ...
        'fuzzy',num2str(dRound(method.cell_dist,0)),'_results.mat'];
elseif method.fuzzy && method.rerank && method.filter
    results_file = ['.\',method.decision,'\rerank\query',num2str(query_set), ...
        'fuzzy',num2str(dRound(method.cell_dist,0)),'_results.mat'];
elseif method.fuzzy && ~method.rerank && ~method.filter
    results_file = ['.\',method.decision,'\query',num2str(query_set), ...
        'fuzzy',num2str(dRound(method.cell_dist,0)),'nofilt_results.mat'];
elseif method.fuzzy && method.rerank && ~method.filter
    results_file = ['.\',method.decision,'\rerank\query',num2str(query_set), ...
        'fuzzy',num2str(dRound(method.cell_dist,0)),'nofilt_results.mat'];
else
    results_file = ['.\',method.decision,'\query',num2str(query_set), ...
        'exact',num2str(dRound(method.cell_dist,0)),'_results.mat'];
end
try
    load(results_file)
catch
    results = struct;
end

% Get a list of queries
query = dir(qDir);
query = strvcat(query(3:end).name);
query = unique(str2double(cellstr(query(:,5:8)))); % query numbers
query(isnan(query)) = [];
nq = length(query); % number of queries

% Load query ratio scores
scores_file = ['.\scores\query',num2str(query_set),'_scores.mat'];
try
    load(scores_file)
catch
    query_scores.scores = cell(nq,1);
    query_scores.number = query;
end

% Post-processing | Main part
if method.fuzzy && strcmp(method.decision,'linear')
    
    % Initialize variables
    wv = method.wv;
    ws = method.ws;
    cell_dist = method.cell_dist;
    if strcmp(method.distribution,'exponential')
        fuzzy_rad = 200;
        spacing = 4;
    else
        fuzzy_rad = 75;
        spacing = 1;
    end
    N = 10; % maximum number of candidate images
    nruns = length(wv);
    
    % Set up results structure and fields
    if ~isfield(results,'wv')
        results.wv = wv;
        results.ws = ws;
        results.match = zeros(nruns,N);
        results.poss = zeros(nruns,N);
        results.total = zeros(nruns,N);
    else % check for results already done
%         alreadydone = [];
%         for k=1:length(results.wv)
%             idxv = find(wv==results.wv(k));
%             idxs = find(ws==results.ws(k));
%             idx = intersect(idxv,idxs);
%             alreadydone = [alreadydone,idx];
%         end
%         if ~isempty(alreadydone)
%             disp('Some data points have already been computed.')
%             wv(alreadydone) = [];
%             ws(alreadydone) = [];
%             nruns = length(wv);
%         end
        results.wv(1,end+1:end+nruns) = wv;
        results.ws(1,end+1:end+nruns) = ws;
        results.match(end+1:end+nruns,:) = zeros(nruns,N);
        results.poss(end+1:end+nruns,:) = zeros(nruns,N);
        results.total(end+1:end+nruns,:) = zeros(nruns,N);
    end
    match = zeros(nruns,N);
    poss = zeros(nruns,N);
    total = zeros(nruns,N);
    
    if cell_dist == 0
        fprintf(['\nRunning fuzzy post processing with linear ', ...
            'decision and no combination...\n'])
    else
        fprintf(['\nRunning fuzzy post processing with linear ', ...
            'decision and vote combination...\n'])
    end
    
    mismatch_file = 'E:\Research\app\code\matlab\post-processing\linear\rerank\mismatch.txt';
    fid = fopen(mismatch_file,'w');
    
    % Iterate through each query
    for k=1:nq
        
        fprintf(['\nProcessing query ',num2str(k),'... '])
        
        % Get ground truth matches
        gt = getGT(query(k),query_set,gtDir);
        
        % Get list of cells and cell locations
        [qCell,cLat,cLon] = getCells(query(k),vDir);
        
        % Get query sift file name and query location
        idx = strfind(qCell{1},',');
        query_name = qCell{1}(1:idx(3)-1);
        qLat = str2double(query_name(idx(1)+1:idx(2)-1));
        qLon = str2double(query_name(idx(2)+1:end-8));
        
        % Get fuzzy point locations
        [lats,lons] = getFuzzyLocs(qLat,qLon,fuzzy_rad,spacing);
        
        % Create cell groupings from fuzzy points and iterate through them
        cellgroups = groupFuzzy(lats,lons,cLat,cLon,cell_dist);
        
        ma = 0;
        to = 0;
        for cg=cellgroups
            
            % Combine cells in this grouping
            comb_cells = qCell(cg.idx);
            [image,vote,img_lat,img_lon] = cellCombine(comb_cells,vDir);
            
            % Iterate through each fuzzy point and post process
            for j=1:cg.npts
                
                % Get the candidate images from the combined votes
                if ~method.filter
                    [cand,cand_vote] = getCand( cg.lats(j),cg.lons(j), ...
                        img_lat,img_lon,image,vote,'none',method.rerank );
                elseif strcmp(method.distribution,'exponential')
                    [cand,cand_vote] = getCand( cg.lats(j),cg.lons(j), ...
                        img_lat,img_lon,image,vote,'exponential',method.rerank );
                else
                    [cand,cand_vote] = getCand( cg.lats(j),cg.lons(j), ...
                        img_lat,img_lon,image,vote,'cutoff',method.rerank );
                end
                
                % Get the candidate ratio scores
                idx = find(query_scores.number==query(k),1,'first');
                img_scores = query_scores.scores{idx};
                if ~iscell(img_scores)
                    img_scores = cell(0,2);
                end
                [cand_score,img_scores] = getRatioScores(...
                    cand,img_scores,query_name,qDir,dbDir);
                query_scores.scores{idx} = img_scores;
                
                % Evaluate the post processing
                [m,p] = evaluateScores(cand_vote,cand_score,cand, ...
                    gt,wv,ws,method.rerank);
                % weigh result if exponential distribution
                if strcmp(method.distribution,'exponential')
                    d = latlonDistance(qLat,qLon,cg.lats(j),cg.lons(j));
                    weight = exp(-d/50);
                else
                    weight = 1;
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
    
end

% store query ratio scores
save(scores_file,'query_scores')
fclose(fid);

% store results
% results.match(end-nruns+1:end,:) = match;
% results.poss(end-nruns+1:end,:) = poss;
% results.total(end-nruns+1:end,:) = total;
% save(results_file,'results')

fprintf('\nRun complete.\n\n')