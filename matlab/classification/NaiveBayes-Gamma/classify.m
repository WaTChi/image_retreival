function [classes,LogL] = classify(features,classifier)

% [classes] = classify(features,classifier)
% 
% First coded 12 Dec 2010 by Aaron Hallquist.
% Latest revision 6 Jan 2011 by Aaron Hallquist.
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
%       .priors:	Cx1 vector of prior probabilities for the classes
%                       1 = sum(priors)
%       .means:     (C+1)xM matrix of feature means
%       .vars:      (C+1)xM matrix of feature variances
% 
% OUTPUT:
%   classes:    Kx1 vector of predicted classes
%   LogL:       KxC vector of log likelihoods

% Get size parameters
[K,M] = size(features);
C = length(classifier.priors);
        
% Iterate through each sample and feature to compute the log likelihood 
% (LogL) that each sample is in each class (KxC matrix)
LogL = repmat( log(classifier.priors)' , [K,1] ); % accounting for priors
for k=1:K
    feat_idx = isfinite(features(k,:));
    % Add in conditional log p(Fi|C)
    x = repmat( features(k,feat_idx) , [C,1] );
    m = classifier.means(1:C,feat_idx);
    v = classifier.vars(1:C,feat_idx);
    prmK = m.^2 ./ v;
    prmT = v ./ m;
    dLogL = (prmK-1).*log(x) - (x./prmT) - ...
            log(gamma(prmK)) - prmK.*log(prmT);
    LogL(k,:) = LogL(k,:) + sum(dLogL,2)';
    % Add in log p(Fi)
    x = features(k,feat_idx);
    m = classifier.means(end,feat_idx);
    v = classifier.vars(end,feat_idx);
    prmK = m.^2 ./ v;
    prmT = v ./ m;
    dLogL = (prmK-1).*log(x) - (x./prmT) - ...
            log(gamma(prmK)) - prmK.*log(prmT);
    LogL(k,:) = LogL(k,:) - repmat( sum(dLogL) , [1,C] );
end

% Compare the log likelihoods and choose the class with the largest value
[~,classes] = max(LogL,[],2);