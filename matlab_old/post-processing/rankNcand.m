function [cand,param,score] = rankNcand(cand,param,method,ncand,N)

% [cand] = rankNcand(method,N,cand,dist,vote,score)
% 
% P = number of parameters
% K = number of input candidate images
% 
% DESCRIPTION
%   This function takes in a list of candidates with corresponding
%   parameters and ranks and returns the top N results based on the values
%   of each parameter and the ranking method
% 
% INPUTS
%   cand:   Kx1 cell array of candidate images
%   param:  KxP matrix of parameter values. The type of parameter in each
%           column is determined by the contents of method
%   method: The query method structure (same as post_process)
%   ncand:  Number of candidates for query
%   N:      The number of candidates to return (optional)
%           - if unspecified, N = K
%   
% 
% OUTPUTS
%   cand:   Nx1 cell array of ranked candidates
%   score:  Nx1 array of output scores, if desired

% Parse method decision
idx = strfind(method.decision,'-');
decision = method.decision(1:idx-1);
decis_prm = method.decision(idx+1:end);

% Set size variables
[K,P] = size(param);
if nargin < 5
    N = K;
else
    N = min(N,K);
end

if strcmp(decision,'vote')
% --------------------------------
% Vote Decision
% --------------------------------
    score = param(:,2);
    
else % if strcmp(decision,'bayes') || strcmp(decision,'nbayes')
% --------------------------------
% Bayes Decision
% --------------------------------
    % Load the classifier
    cls = getClassifier(method,0);
    cls_idx = 1 + sum( ncand > method.canddiv );
    classifier = cls{cls_idx};
    
    % Set up the feature matrix
    if P<2
        param(:,end+1:2) = NaN;
    end
    if isempty(strfind(decis_prm,'d'))
        param(:,1) = nan(K,1);
    end
    if isempty(strfind(decis_prm,'v'))
        param(:,2) = nan(K,1);
    end
    
    % Get the conditional probabilities and sort
    [~,score] = classifybayes(param,classifier);
     
end

score = [score(:,1),param(:,2),param(:,1)];
[score,rank_idx] = sortrows(score,[-1,-2,+3]);
cand = cand(rank_idx(1:N));
param = param(rank_idx(1:N),:);