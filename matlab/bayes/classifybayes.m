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
% number of training samples in its associated bin (as fraction of total
% samples)
if nargin < 3
    min_samp = 1e-7;
end
total = sum(classifier.nsamps);

% Get bin information
bins = classifier.bins;
sp = classifier.spacing;
B = size(bins); B = B(1:2);

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
    
    % Get features and associated bins
    f1 = features(j,1);
    if isnan(f1)
        sb1 = 1; lb1 = B(1);
    else
        sb1 = max(min(ceil(f1/sp(1)),B(1)),1); lb1 = sb1;
    end
    f2 = features(j,2);
    if isnan(f2)
        sb2 = 1; lb2 = B(2);
    else
        sb2 = max(min(ceil(f2/sp(2)),B(2)),1); lb2 = sb2;
    end
    
    % Expand bin region to satisfy min_samp
    frct1 = (lb1-sb1+1)/B(1);
    frct2 = (lb2-sb2+1)/B(2);
    count = reshape( sum( sum( bins(sb1:lb1,sb2:lb2,:) , 1 ) , 2 ) , [1,2] );
    while sum(count) < min_samp * total / (frct1*frct2)
        frct1 = (lb1-sb1+1)/B(1);
        frct2 = (lb2-sb2+1)/B(2);
        if frct2 > frct1
            sb1 = max(1,sb1-1);
            lb1 = min(B(1),lb1+1);
        else
            sb2 = max(1,sb2-1);
            lb2 = min(B(2),lb2+1);
        end
        count = reshape( sum( sum( bins(sb1:lb1,sb2:lb2,:) , 1 ) , 2 ) , [1,2] );
    end
    
    % Factor in feature conditional probability
    prob = prob .* ( count ./ nsamps );
    
    % Get the class and normalize the likelihoods to get condP
    [~,class(j)] = max(prob);
    condP(j,:) = prob / sum(prob);
    
end