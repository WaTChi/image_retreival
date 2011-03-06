function [score, Q] = nnRatioN(p1, p2, dataset)

dir = ['.\', dataset, '\']; % must end in file separator (e.g. '\'); this 
            % is an input directory that must consist of *.mat files 
            % specifying SIFT features according to a format that Matt 
            % Carlberg was using

queryDes = extractDescriptors([dir, p1]);
candDes = extractDescriptors([dir, p2]);

Q = size(queryDes, 1);
score = 0;
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
        score = score+1;
    end
end
%score = score/size(queryDes, 1);