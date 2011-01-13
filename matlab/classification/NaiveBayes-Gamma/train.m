function [classifier] = train(features,classes)

% [classifier] = train(features,classes)
% 
% First coded 12 Dec 2010 by Aaron Hallquist.
% Latest revision 12 Dec 2010 by Aaron Hallquist.
%
% K = number of training samples
% M = number of features
% C = number of classes
% 
% DESCRIPTION:
%   This function trains a Naive Bayes classifier with a Gamma distribution
%   assumption on the features. It takes as input a list of K training 
%   samples with their features and classes and outputs the priors and 
%   feature distributions necessary to classify new samples. If a feature
%   is absent for a training sample, this is denoted with the value NaN.
% 
% INPUT:
%   features:   KxM matrix of features
%   classes:    Kx1 vector of classes
%                   C = max(classes)
%                   classes must contain integers ranging from 1 to C
% 
% OUTPUT:
%   classifier: Structure with the following elements...
%       .priors:	Cx1 vector of prior probabilities for the classes
%                       1 = sum(priors)
%       .means:     (C+1)xM matrix of feature means
%                       means(end) = mean for all training samples
%       .vars:      (C+1)xM matrix of feature variances
%                       vars(end) = variance for all training samples

% Get size parameters
[K,M] = size(features);
C = max(classes);

% Initialize classifier
classifier.priors = zeros(C,1);
classifier.means = zeros(C+1,M);
classifier.vars = zeros(C+1,M);

% Iterate through classes to compute classifier information
for c=1:C
    class_idx = ( c == classes );
    class_features = features(class_idx,:);
    classifier.priors(c) = sum(class_idx) / K;
    % Iterate through each feature to compute its mean and variance
    for m=1:M
        feat = class_features( isfinite(class_features(:,m)) , m );
        classifier.means(c,m) = mean(feat);
        classifier.vars(c,m) = var(feat);
    end
end

% Compute the mean and variance of the features regardless of class
for m=1:M
    feat = features( isfinite(features(:,m)) , m );
    classifier.means(end,m) = mean(feat);
    classifier.vars(end,m) = var(feat);
end % not necessary, but useful to have