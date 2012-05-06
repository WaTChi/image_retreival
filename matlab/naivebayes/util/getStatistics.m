function [stat] = getStatistics(samples)

% [stat] = getStatistics(samples)
% 
% First coded 12 Dec 2010 by Aaron Hallquist.
% Latest revision 12 Dec 2010 by Aaron Hallquist.
%
% N = number of samples from which to retrieve statistics
% 
% DESCRIPTION:
%   This function computes the statistics of N samples of a random
%   variable. It retrives the mean, standard deviation, skewness, and
%   kurtosis for the sample distribution.
% 
% INPUT:
%   samples:    Nx1 array of samples
% 
% OUTPUT:
%   stat:       Structure containing the following statistics...
%       mean:   Sample mean.
%       stdv:   Sample standard deviation.
%       skew:   Sample skewness normalized by variance.
%       kurt:   Sample kurtosis normalized by variance.

ns = length(samples);
mu = mean(samples);
sd = std(samples,1);

stat.mean = mu;
stat.stdv = sd;
stat.skew = (1/ns) * sum( (samples-mu).^3 ) / sd^3;
stat.kurt = -3 + (1/ns) * sum( (samples-mu).^4 ) / sd^4;