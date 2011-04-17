function [ratio] = thresDist(i, dataset)

dir = ['.\', dataset, '\']; % must end in file separator (e.g. '\'); this 
            % is an input directory that must consist of *.mat files 
            % specifying SIFT features according to a format that Matt 
            % Carlberg was using

siftFeatures = ls(fullfile(dir, '*.feat')); 
N = size(siftFeatures, 1);
siftFeatures = [repmat(dir, N, 1), siftFeatures];

queryDes = extractDescriptors(siftFeatures(i, :));
Q = size(queryDes, 1); 
ratio = zeros(1, N-1); 

x = 1;
for j = 1:N
    score = 0;
    fprintf('on image %d of %d\n', j, N);
    if j == i
        continue
    end
    candDes = extractDescriptors(siftFeatures(j, :));
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
        if nn/nn2 < 0.6
            score = score + 1;
        end
    end
    ratio(x) = score/Q;
    x = x+1;
end

hist(ratio);