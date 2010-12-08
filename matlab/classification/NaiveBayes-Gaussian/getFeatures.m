function [features] = getFeatures(distributions)

% [features] = getFeatures(distributions)
% 
% First coded 12 Dec 2010 by Aaron Hallquist.
% Latest revision 12 Dec 2010 by Aaron Hallquist.
%
% K = number of samples from which to retrieve features
% M = number of features
% N = vote threshold for feature (2)
% 
% DESCRIPTION:
%   This function reads the distributions from K samples and computes the
%   features used in classification for those samples. The features
%   computed in this version are:
%       1) Total number of votes
%       2) Number of images with less than N votes.
%       3) Max vote as a percentage of total votes
%       4) Mean vote as a percentage of total votes
%       5) Standard deviation of votes as a percentage of total votes
%       6) Skew of votes as a percentage of total votes
%       7) Kurtosis of votes as a percentage of total votes

%       
% 
% INPUT:
%   distributions:  Kx1 cell array with each cell containing the sample's
%                   distribution of votes
% 
% OUTPUT:
%   features:       KxM matrix of features for each sample

% Get and set parameters
K = length(distributions);
M = 7;
N = 10;

% Initialize the feature matrix
features = zeros(K,M);

% Iterate through the distributions to compute features
for k=1:K
    votes = distributions{k};
    features(k,1) = sum( votes ); %---- Total number of votes ----%
    features(k,2) = sum( votes < N ); %---- Quantity of low votes ----%
    votes = votes / features(k,1); % scale by number of votes
    features(k,3) = max(votes); %---- Max vote ----%
    stat = getStatistics(votes); % get distribution statistics
    features(k,4) = stat.mean; %---- Mean vote ----%
    features(k,5) = stat.stdv; %---- Standard deviation of votes ----%
    features(k,6) = stat.skew; %---- Skewness of votes ----%
    features(k,7) = stat.kurt; %---- Kurtosis of votes ----%
end
    