function [class,condP] = classifybayes(features,classifier,dists)

% [class,condP] = classifybayes(features,classifier,dists)
% 
% First coded 12 Dec 2010 by Aaron Hallquist.
% Latest revision 28 Jan 2011 by Aaron Hallquist.
%
% K = number of samples to classify
% M = number of features
% C = number of classes
% 
% DESCRIPTION:
%   This function classifies a list of feature samples using a pre-trained
%   Naive Bayes classifier with a Gamma assumption on the distibution of
%   the features. It takes as input a list of K samples to classify as well
%   as a structure containing the pre-trained classifier information, and
%   gives as output a vector containing the predicted classes. If a feature
%   is missing, this is denoted by the value NaN.
% 
% INPUT:
%   features:   KxM matrix of features
%   classifier: Structure with the following elements...
%       .dists:     1xM cell array of strings containing the assumed
%                       distribution for each feature
%       .priors:	Cx1 vector of prior probabilities for the classes
%                       1 == sum(priors)
%       .nsamps:    Cx1 vector containing the number of training samples in
%                       each class:     prios = nsamps / sum(nsamps)
%       .means:     CxM matrix of feature means
%       .vars:      CxM matrix of feature variances
%       .nfeats:    CxM matrix containing the number of features for each
%                       class used to compute means and variances
%   dists:      Optional cell array input to override the training assumed
%               distributions (note that training is still valid)
% 
% OUTPUT:
%   class:      Kx1 vector of predicted classes for each feature sample
%   condP:      KxC vector of conditional probabilities P(C|F1..Fn) for
%               feature sample : ones(K,1) == sum(condP,2)
%                   - This is not computed if nargin < 2

% Get size parameters
[J,M] = size(features);
C = length(classifier.priors);

% Set distribution types and initialize priors
if nargin < 3
    if isfield(classifier,'dists')
        dists = classifier.dists;
    else
        dists = repmat( cellstr('Gamma') , [1,M] );
    end
end
priors = classifier.priors';

% ------------------------------------------------------------
% Compute conditional probabilities and predicted classes
% ------------------------------------------------------------
% Iterate through each feature sample and compute its class and likelihoods
class = ones(J,1);
condP = zeros(J,C);
for j=1:J

    % Initialize the log likelihood and denominator prod( p(Fi) )
    LogL = log(priors); % accounting for priors
    
    % Iterate through each *known* feature and update the likelihoods
    feat_idx = find( isfinite(features(j,:)) ); % NaN indicates unknown
    for f = feat_idx
        
        % Get the trained statistics and feature value
        x = features(j,f);
        m = classifier.means(:,f)'; % mean
        v = classifier.vars(:,f)'; % variance

        % Compute log p(Fi|Cj) depending on the distribution
        if strcmp(dists{f},'Gamma')
            k = m.^2 / v; % shape parameter k
            t = v ./ m; % scale parameter theta
            logpfc = (k-1)*log(x) - x./t - log(gamma(k)) - k.*log(t);
        else % if strcmp(dists{f},'Normal')
            logpfc = (-1/2)*log(2*pi*v) - (x-m).^2./(2*v);
        end

        % Update log likelihoods and denominator sum log( p(Fi) )
        LogL = LogL + logpfc;

    end
    
    % Compute maximum likelihood class
    [~,class(j)] = max(LogL);

    % Compute conditional probabilities
    total = log( sum( exp(LogL) ) );
    condP(j,:) = exp(LogL-total);

end