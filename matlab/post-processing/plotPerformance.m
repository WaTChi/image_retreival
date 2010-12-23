function results = plotPerformance(filter,comb_method,query_num,fignum)

% results = plotPerformance(filter,comb_method,query_num,fignum)
% 
% First coded 9 Dec 2010 by Aaron Hallquist
% Latest revision 9 Dec 2010 by Aaron Hallquist
% 
% DESCRIPTION
%   This function plots the performance of post-processing.
% 
% INPUT
%   filter:         String choosing filter for obtaining candidate images
%   comb_method:    String choosing the combination method, one of...
%           'no-comb':      Nearest cell only
%           'dup-comb':     Finds duplicate images across multiple cells
%           'vote-comb':    Sums vote results across multiple cells
%   query_num:      Integer choosing the query set to use
%   figum:          Integer choosing where to begin the figures
% 
% OUTPUT
%   results:        Structure containing post-processing data

addpath('.\..\util\')

% Default values
if nargin<4
    fignum=0;
end

% Run post-processing if necessary
file = ['.\',comb_method,'\query',num2str(query_num),'_results.mat'];
try
    load(file)
catch e
    disp('Post-processing not complete. Running this now...')
    results = post_process(filter,comb_method,query_num);
end

getResults = true;
% Evaluate performance if necessary
if getResults || ~isfield(results,'performance')
    disp('Performance not yet evaluated. Running this now...')
    results = getPerformance(filter,comb_method,query_num);
end

% Plot or display performance
if strcmp(filter,'upto10')
    nq = results.numqueries;
    gNum = results.performance.gNum;
    yNum = results.performance.yNum;
    figure(fignum+2)
    hold off
    bar(1:10,gNum+yNum,'FaceColor','y')
    hold on
    bar(1:10,gNum,'FaceColor','g')
    axis([0.5 10.5 0 nq])
    title('Post-processing performance on top N')
    xlabel('Top N')
    ylabel('Number of successful matches; Y=possible, G=actual')
    [~,bestIdx] = max(gNum); bestIdx = min(bestIdx);
elseif strcmp(filter,'above1')
    nq = results.numqueries;
    gNum = results.performance.gNum;
    disp(' ')
    disp(['Post-processing succeeds on ',num2str(gNum),'/',num2str(nq),...
          ' queries (',num2str(dRound(100*gNum/nq,0)),'%).'])
      disp(' ')
    bestIdx = 1;
end

% Plot min/max scores of good and bad queries
gIdx = results.performance.gQuery{bestIdx,2};
yIdx = results.performance.yQuery{bestIdx,2};
rIdx = results.performance.rQuery{bestIdx,2};
minScore = results.performance.minScore(:,bestIdx);
maxScore = results.performance.maxScore(:,bestIdx);
figure(fignum+1)
hold off
bar(gIdx,maxScore(gIdx),'FaceColor','g')
hold on
bar(gIdx,minScore(gIdx),'FaceColor','w','EdgeColor','w')
bar(yIdx,maxScore(yIdx),'FaceColor','y')
bar(yIdx,minScore(yIdx),'FaceColor','w','EdgeColor','w')
bar(rIdx,maxScore(rIdx),'FaceColor','r')
bar(rIdx,minScore(rIdx),'FaceColor','w','EdgeColor','w')
title('Post-processing score range for queries')
xlabel('Queries; G=success, Y=possible, R=impossible')
ylabel('Post-processing score range')