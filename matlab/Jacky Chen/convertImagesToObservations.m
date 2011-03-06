function [] = convertImagesToObservations(dataset, threshold)

%-----------------------------------------------------------------------------
% convertImagesToObservations
% coded by Jacky Chen
%
% August 5, 2010
%
% This script is based on the following paper:
%
%   M. Cummins and P. Newman. "Probabilistic Appearance Based Navigation 
%   and Loop Closing." International Conference on Robotics and Automation,
%   2007.
%
% Please refer to the above paper (should be included as a PDF in the same
% directory as this script) for reading this script. Note that this
% script tries to stick fairly close to the loop detection algorithm in the
% paper.
%
% Input:
% 
%   dataset: string that specifies the directory and the file names.
%   Typically of the form "date_dataset"
%   threshold: number of features that all images in the dataset should
%   have. This program will prune/delete those that have fewer than this.
%
% Output:
%
%   A .mat file containing observations for each image in the directory.
%   Filename is typically of the form "obs'date'_'dataset'"
%

dir = ['.\', dataset, '\']; % must end in file separator (e.g. '\'); this 
            % is an input directory that must consist of *.mat files 
            % specifying SIFT features according to a format that Matt 
            % Carlberg was using

siftFeatures = ls(fullfile(dir, '*.feat')); %   assumes files containing 
            % SIFT features are .feat files
siftFeatures = [repmat(dir, size(siftFeatures, 1), 1), siftFeatures];

N = size(siftFeatures, 1); 

fprintf('Current threshold requires images to have at least %d features\n', threshold);
for i = 1:N
    if numOfFeatures(siftFeatures(i,:)) < threshold
        fprintf('Excluding file %s\n', siftFeatures(i,:));
        delete(siftFeatures(i, :));
    end
end

siftFeatures = ls(fullfile(dir, '*.feat')); 
siftFeatures = [repmat(dir, size(siftFeatures, 1), 1), siftFeatures];
N = size(siftFeatures, 1); % N = total number of documents/images after throwing images with few features

M = 100;                   % M = number of documents/images to train on

%% First build a vocabulary

if M > N
    M = N;
end

samples       = randsample(N, M);
clusters      = zeros(0, 128); % the average of features within the cluster
clusterRadius = 300;
clusterSizes  = zeros(0, 1);

fprintf('Training on %d (out of %d) images...\n', M, N);

for i = 1:M
    if rem(i,10) == 0
        fprintf(' Adding images...%d images done\n', i);
    end

    featureDescriptors = extractDescriptors(siftFeatures(samples(i), :));

    for j = 1:size(featureDescriptors, 1)
        featureDescriptor = featureDescriptors(j, :);
        done              = false;

        for k = 1:size(clusters, 1)
            % if the current feature belongs to this cluster
            if norm(featureDescriptor - clusters(k, :)) < clusterRadius 
                clusters(k, :)  = clusterSizes(k) * clusters(k, :) + ...
                                  featureDescriptor;
                clusterSizes(k) = clusterSizes(k) + 1;
                clusters(k, :)  = clusters(k, :) / clusterSizes(k);
                done            = true;
                break
            end
        end

        if ~done
            clusters     = [clusters; featureDescriptor];
            clusterSizes = [clusterSizes; 1];
        end
    end
end

K = size(clusters, 1); % K = number of clusters

fprintf('Vocabulary size: %d\n', K);

%% Now go through images and find binary features for all images
observations   = zeros(N, K);

fprintf('Finding binary features for all images (%d images)...\n', N);

p = 10;
for i = 1:N
    if rem(i,10) == 1
        if i+9 > N
            fprintf(' Processing images %d-%d/%d...\n', i, N, N);
        else
            fprintf(' Processing images %d-%d/%d...\n', i, i+9, N);
        end
    end
    if p*N < 100*i
        fprintf(' -- %d%% done -- \n', p);
        p = p+10;
    end

    featureDescriptors = extractDescriptors(siftFeatures(i, :));

    for j = 1:size(featureDescriptors, 1)
        featureDescriptor = featureDescriptors(j, :);

        for k = 1:K
            if norm(featureDescriptor - clusters(k, :)) < clusterRadius
                observations(i,k) = 1;
                break
            end
        end
    end
end

save(['obs', dataset], 'observations');

cleanseData(['obs', dataset, '.mat']);
