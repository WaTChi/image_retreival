function [results] = getPerformance(filter,comb_method,query_num)

% [results] = getPerformance(filter,comb_method,query_num)
% 
% First coded 9 Dec 2010 by Aaron Hallquist
% Latest revision 9 Dec 2010 by Aaron Hallquist
% 
% DESCRIPTION
%   This function analyzes the post-processing results and computes
%   performance statistics. The analysis is different for each run
%   depending on the combination method and candidate filter used.
% 
% INPUT
%   filter:         String choosing filter for obtaining candidate images
%   comb_method:    String choosing the combination method, one of...
%           'no-comb':      Nearest cell only
%           'dup-comb':     Finds duplicate images across multiple cells
%           'vote-comb':    Sums vote results across multiple cells
%   query_num:      Integer choosing the query set to use
% 
% OUTPUT
%   results:  	Structure containing the details of post-processing
%       - type 'help ppStruct' for details about its contents

file = ['.\',comb_method,'\query',num2str(query_num),'_results.mat'];
try
    load(file)
catch e
    disp('Post-processing not complete. Running this now...')
    results = post_process(filter,comb_method,query_num);
end

if strcmp(filter,'upto10')
    
    results
    
    % Iterate through each query to get its performance
    nq = results.numqueries; % number of queries
    winner = cell(nq,10);
    minScore = zeros(nq,10);
    maxScore = zeros(nq,10);
    gQuery = cell(10,2);
    yQuery = cell(10,2);
    rQuery = cell(10,2);
    for k=1:nq
        
        % Get candidate information
        candFiles = results.cand.files{k};
        candVotes = results.cand.votes{k};
        candScores = results.cand.scores{k};
        
        % Get ground truth information
        gtFiles = results.gtruth{k};
        
        % Iterate through each n in top N analysis
        for N=(10:-1:1)
            
            % Reduce the size of the candidates based on N
            if N<length(candVotes) && candVotes(N)>candVotes(end)
                candFiles = candFiles(1:N);
                candVotes = candVotes(1:N);
                candScores = candScores(1:N);
            end
            
            minScore(k,N) = min(candScores);
            [maxScore(k,N),winnerIdx] = max(candScores);
            winnerIdx = min(winnerIdx);
            winner(k,N) = candFiles(winnerIdx);
            
            if textMatch(winner(k,N),gtFiles)
                gQuery{N,1} = [gQuery{N,1}, results.num(k)];
                gQuery{N,2} = [gQuery{N,2}, k];
            elseif textMatch(candFiles,gtFiles)
                yQuery{N,1} = [yQuery{N,1}, results.num(k)];
                yQuery{N,2} = [yQuery{N,2}, k];
            else
                rQuery{N,1} = [rQuery{N,1}, results.num(k)];
                rQuery{N,2} = [rQuery{N,2}, k];
            end
            
        end
        
    end
    
    % Set the results performance structure
    gNum = zeros(10,1); gRatio = zeros(10,1);
    yNum = zeros(10,1); yRatio = zeros(10,1);
    rNum = zeros(10,1);
    for N=1:10
        gNum(N) = length(gQuery{N,1});
        yNum(N) = length(yQuery{N,1});
        rNum(N) = length(rQuery{N,1});
        gRatio(N) = gNum(N) / nq;
        yRatio(N) = ( gNum(N)+yNum(N) ) / nq;
    end
    results.performance.gRatio = gRatio;
    results.performance.yRatio = yRatio;
    results.performance.gQuery = gQuery;
    results.performance.yQuery = yQuery;
    results.performance.rQuery = rQuery;
    results.performance.winner = winner;
    results.performance.gNum = gNum;
    results.performance.yNum = yNum;
    results.performance.rNum = rNum;
    results.performance.minScore = minScore;
    results.performance.maxScore = maxScore;
    
elseif strcmp(filter,'above1')
    
    % Iterate through each query to get its performance
    nq = results.numqueries; % number of queries
    winner = cell(nq,1);
    minScore = zeros(nq,1);
    maxScore = zeros(nq,1);
    gQuery = cell(1,2);
    yQuery = cell(1,2);
    rQuery = cell(1,2);
    for k=1:nq
        
        % Get candidate information
        candFiles = results.cand.files{k};
        candVotes = results.cand.votes{k};
        candScores = results.cand.scores{k};
        
        % Get ground truth information
        gtFiles = results.gtruth{k};

        minScore(k) = min(candScores);
        [maxScore(k),winnerIdx] = max(candScores);
        winnerIdx = min(winnerIdx);
        winner(k) = candFiles(winnerIdx);

        if textMatch(winner(k),gtFiles)
            gQuery{1} = [gQuery{1}, results.num(k)];
            gQuery{2} = [gQuery{2}, k];
        elseif textMatch(candFiles,gtFiles)
            yQuery{1} = [yQuery{1}, results.num(k)];
            yQuery{2} = [yQuery{2}, k];
        else
            rQuery{1} = [rQuery{1}, results.num(k)];
            rQuery{2} = [rQuery{2}, k];
        end
        
    end
    
    % Set the results performance structure
    gNum = length(gQuery{1});
    yNum = length(yQuery{1});
    rNum = length(rQuery{1});
    gRatio = gNum / nq;
    yRatio = ( gNum+yNum ) / nq;
    results.performance.gRatio = gRatio;
    results.performance.yRatio = yRatio;
    results.performance.gQuery = gQuery;
    results.performance.yQuery = yQuery;
    results.performance.rQuery = rQuery;
    results.performance.winner = winner;
    results.performance.gNum = gNum;
    results.performance.yNum = yNum;
    results.performance.rNum = rNum;
    results.performance.minScore = minScore;
    results.performance.maxScore = maxScore;
    
end

save(file,'results')