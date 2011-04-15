function [ratio] = nnDist(i, j, dataset)

dir = ['.\', dataset, '\']; % must end in file separator (e.g. '\'); this 
            % is an input directory that must consist of *.mat files 
            % specifying SIFT features according to a format that Matt 
            % Carlberg was using

siftFeatures = ls(fullfile(dir, '*.feat')); 
siftFeatures = [repmat(dir, size(siftFeatures, 1), 1), siftFeatures];

queryDes = extractDescriptors(siftFeatures(i, :));
candDes = extractDescriptors(siftFeatures(j, :));

Q = size(queryDes, 1); 
ratio = zeros(1, Q); 
for k = 1:Q
    featureDescriptor = queryDes(k, :);
    nn = 1000;
    nn2 = 1000;
    for l = 1:size(candDes, 1)
        dist = norm(featureDescriptor - candDes(l,:));
        if dist < nn
            nn2 = nn;
            nn = dist;
        else if dist < nn2
                nn2 = dist;
            end
        end
    end
    ratio(k) = nn/nn2;
end

hist(ratio);