function [class,condP] = classifybayes(features,classifier,min_samp)

% [class,condP] = classifybayes(features,classifier,min_samp)
% 
% First coded 12 Dec 2010 by Aaron Hallquist.
% Latest revision 12 Feb 2011 by Aaron Hallquist.
%
% J = number of samples to classify
% M = number of features
% C = number of classes
% B(m) = number of bins for a given feature
% 
% DESCRIPTION:
%   This function classifies a list of feature samples using a pre-trained
%   Naive Bayes classifier using bins to learn the distributions. In
%   classifying a feature, if the feature bin does not have an adequate 
%   number of training samples for a class, the bins are expanded for each 
%   class until a minimum number of samples is found in each class bin.
% 
% INPUT:
%   features:   KxM matrix of features
%   classifier: Structure with the following elements...
%       .spacing:   1xM vector of bin spacings for each feature
%       .nsamps:    1xC vector containing the number of training samples in
%                       each class:     prios = nsamps / sum(nsamps)
%       .bins:      1xM cell array of distribution counts. Cells contain...
%                   - B(m) x C vector containing the counts of each
%                     training sample in each bin for each class
%                   - bins start at 0 and are separated by spacing(m)
%   min_samp:   Minimum number of samples to stop bin growth (default = 10)
% 
% OUTPUT:
%   class:      Jx1 vector of predicted classes for each feature sample
%   condP:      JxC vector of conditional probabilities P(C|F1..Fn) for
%               feature sample : ones(J,1) == sum(condP,2)

% minimum trained samples; bins expand until each class sees min_samp
% number of training samples in its associated bin
if nargin < 3
    min_samp = 100;
end

% Get size parameters
[J,M] = size(features);
C = length(classifier.nsamps);

% ------------------------------------------------------------
% Compute conditional probabilities and predicted classes
% ------------------------------------------------------------
% Iterate through each feature sample and compute its class and likelihoods
class = ones(J,1);
condP = zeros(J,C);
for j=1:J
    
    % Initialize often used variables
    nsamps = classifier.nsamps;
    
    % Initialize the likelihoods with the priors
    prob = nsamps / sum(nsamps);
    
    % Iterate through each known feature and update the likelihoods
    feat_idx = find(isfinite(features(j,:))); % NaN indicates unknown
    for f = feat_idx
        
        % Get the spacing and bins
        bins = classifier.bins{f};
        B = size(bins,1);
        sp = classifier.spacing(f);
        
        % Get the feature value and associated bin
        x = features(j,f);
        b = max(min(ceil(x/sp),B),1);
        
        % Find the size of the bin necessary to satisfy min_samp
        sb = b; lb = b; % small and large bins
        count = sum( bins(sb:lb,:) , 1 );
        while sum(count) < min_samp
            sb = max(1,sb-1);
            lb = min(B,lb+1);
            count = sum( bins(sb:lb,:) , 1 );
        end
        
        % Factor in feature conditional probability
        prob = prob .* ( count ./ nsamps );
        
    end
    
    % Get the class and normalize the likelihoods to get condP
    [~,class(j)] = max(prob);
    condP(j,:) = prob / sum(prob);
    
end