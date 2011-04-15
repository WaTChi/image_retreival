function [match,poss] = evaluateScores(vote,score,cand,gt,wv,ws,rrN)

% [match,possible] = post_process(score,vote,cand,gt,ws,wv,rrN)
% 
% N = maximum number of candidate images
% 
% DESCRIPTION
%   This function processes the votes and scores of each candidate image to
%   determine whether post processing has succeeded or not.
% 
% INPUTS
%   vote:   Vector of candidate image vote totals
%   score:  Vector of candidate image ratio scores
%   cand:   Candidate image filenames
%   gt:     Vector of ground truth matches
%   wv:     Weight on vote total
%   ws:     Weight on ratio score
%   rrN:    Boolean flag for reranking the top N or filtering the top N
% 
% OUTPUTS
%   match:  Nx1 vector detailing the performance of this query for post
%           processing on the top 1 through the top N (1=success,0=fail).
%           If rrN = true then this details the performance of this query
%           on a rerank of the top 10 with post processing.
%   poss:   Nx1 vector telling us if a match is contained in the candidate
%           images for top 1 through top N (1=yes,0=no). If rrN = true then
%           this output is meaningless.

N = 10; % maximum number of candidate images

if nargin < 5
    ws = 1;
    wv = 0;
elseif nargin < 6
    ws = 0;
end
if nargin < 7
    rrN = false;
end
nruns = length(wv);

match = zeros(nruns,N);
poss = zeros(nruns,N);
if rrN
    v = max(vote);
    s = max(score);
    for j=1:nruns
        final_score = ( wv(j)*(vote/v) + ws(j)*(score/s) ) / (wv(j)+ws(j));
        [~,idx] = sort(final_score,'descend');
        rerank = cand(idx);
        for k=1:N
            match(j,k) = textMatch(rerank(1:k),gt);
        end
    end
else
    for k=(N:-1:1)
        
        if k<length(vote) && vote(k)>vote(end)
            cand = cand(1:k);
            vote = vote(1:k);
            score = score(1:k);
        end
        poss(:,k) = textMatch(cand,gt);
        v = max(vote);
        s = max(score);

        for j=1:nruns
            final_score = ( wv(j)*(vote/v) + ws(j)*(score/s) ) / (wv(j)+ws(j));
            [~,idx] = max(final_score);
            idx = min(idx);
            winner = cand(idx);
            match(j,k) = textMatch(winner,gt);
        end
        
    end
end