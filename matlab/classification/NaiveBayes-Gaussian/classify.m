function [classes] = classify(features,classifier)

% [classes] = classify(features,classifier)
% 
% First coded 12 Dec 2010 by Aaron Hallquist.
% Latest revision 12 Dec 2010 by Aaron Hallquist.
%
% K = number of samples to classify
% M = number of features
% C = number of classes
% 
% DESCRIPTION:
%   This function classifies a list of feature samples using a pre-trained
%   Naive Bayes classifier with a Gaussian assumption on the distibution of
%   the features. It takes as input a list of K samples to classify as well
%   as a structure containing the pre-trained classifier information, and
%   gives as output a vector containing the predicted classes.
% 
% INPUT:
%   features:   KxM matrix of features
%   classifier: Structure with the following elements...
%       .priors:	Cx1 vector of prior probabilities for the classes
%                       1 = sum(priors)
%       .means:     CxM matrix of feature means
%       .vars:      CxM matrix of feature variances
% 
% OUTPUT:
%   classes:    Kx1 vector of classes

% Get size parameters
[K,M] = size(features);
C = length(classifier.priors);

% Iterate through each sample and predict its class
for k=1:K
    % Iterate through each class to compute its LLR
    LLR = log(classifier.priors);
    for c=1:C
        LLR(
        
        
% Iterate through each feature to compute the log likelihood (LogL)
% that each sample is in each class (KxC matrix)
LogL = repmat( log(classifier.priors)' , [K,1] ); % accounting for priors
for k=1:M
    x = repmat( features(:,1) , [1,C] ); % KxC matrix of feature values
    m = repmat( classifier.means(:,k)' , [K,1] ); % KxC matrix of means
    v = repmat( classifier.vars(:,k)' , [K,1] ); % KxC matrix of variances
    LogL = LogL - (1/2) * (x-m).^2  ./  v ;
end

% Compare the log likelihoods and choose the class with the largest value
[~,classes] = max(LogL,[],2);