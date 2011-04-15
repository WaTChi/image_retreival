function [classifier] = trainbayes(features,classes,classifier_info,weights)

% [classifier] = trainbayes(features,classes,classifier_info,weights)
% 
% First coded 12 Dec 2010 by Aaron Hallquist.
% Latest revision 12 Feb 2011 by Aaron Hallquist.
%
% K = number of training samples
% M = number of features
% C = number of classes
% B(m) = number of bins for a given feature
% 
% DESCRIPTION:
%   This function trains a Naive Bayes classifier. It takes as input a list
%   of K training samples with their features and classes and outputs the 
%   priors and feature distributions necessary to classify new samples. If 
%   a feature is absent for a training sample, this is denoted with the 
%   value NaN. The feature distributions are trained using evenly spaced 
%   bins, with spacing specified by by the input classifier_info.
% 
% INPUT:
%   features:   KxM matrix of features
%   classes:    Kx1 vector of classes
%                 	C == max(classes)
%                 	classes must contain integers ranging from 1 to C
%   bayes_info: Input which can be a structure or an Mx1 vector
%               - If it is a vector, the input is assumed to be an Mx1
%                 vector containing the bin spacings for each feature.
%               - If it is a structure, the input is assumed to be a
%                 classifier which we are further training.
%   weights:    Kx1 vector of sample weights (influence of sample)
%               - This is optional. Default = 1 for all samples.
% 
% OUTPUT:
%   classifier: Structure with the following elements...
%       .spacing:   1xM vector of bin spacings for each feature
%       .nsamps:    1xC vector containing the number of training samples in
%                       each class:     prios = nsamps / sum(nsamps)
%       .bins:      1xM cell array of distribution counts. Cells contain...
%                   - B(m) x C vector containing the counts of each
%                     training sample in each bin for each class
%                   - bins start at 0 and are separated by spacing(m)

% Get size parameters
[K,M] = size(features);
C = max(classes);

% Set weights if unspecified
if nargin < 4
    weights = ones(K,1);
end

% Set classifier initial state
if isstruct(classifier_info)
    classifier = classifier_info;
    if ~isfield(classifier,'nsamps')
        classifier.nsamps = zeros(1,C);
        classifier.bins = cell(1,M);
        classifier.bins(:) = {zeros(0,C)};
    end
    C = length(classifier.nsamps);
    M = length(classifier.spacing);
else % numberic vector
    classifier.spacing = classifier_info;
    classifier.nsamps = zeros(1,C);
    classifier.bins = cell(1,M);
    classifier.bins(:) = {zeros(0,C)};
end

% Iterate through features to update classifier
for m=1:M
    
    % Get feature specifics
    spacing = classifier.spacing(m);
    try
    bins = classifier.bins{m};
    catch e
        classifier.bins
        throw(e)
    end
    B = size(bins,1);
    feat = features(:,m);
    
    % Update the size of bins if necessary
    if max(feat) > B*spacing
        B = ceil(max(feat)/spacing);
        bins(end+1:B,:) = 0;
    end
    
    % Iterate through each class to update classifier
    for c=1:C
        
        % Get class features
        cidx = ( c == classes );
        f = feat(cidx);
        
        % Update classifier if there are features associated with class
        if ~isempty(f)
            
            % Get weights and number of features
            w = weights(cidx);
            nf = length(f);

            % Get feature bins and accumulate results
            fbin = ceil(f/spacing);
            b = 1:B;
            bin_add = double( repmat(fbin,[1,B]) == repmat(b,[nf,1]) )' * w;
            bins(:,c) = bins(:,c) + bin_add;

            % Add to nsamps if the first iteration
            if m==1
                classifier.nsamps(c) = classifier.nsamps(c) + sum(w);
            end
            
        end
        
    end
    
    % Update classifier
    classifier.bins{m} = bins;
        
end