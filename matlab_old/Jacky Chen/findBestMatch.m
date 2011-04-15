% threshold - does not consider images with too low of a threshold

function [u, v] = findBestMatch(candidate, query, range, datasetC, threshold)

u = candidate{1}; % to extract char from cell

dir = ['.\', datasetC, '\'];

max = 0;
maxFeat = 0;
R = generateRange(query, range);
for r = 1:2*range+1
    fprintf('range %d\n', r);
    tmp = [dir, R(r,:)];
    if exist(tmp, 'file')
        [score1, Q1] = nnRatioN(u, R(r,:), datasetC);
        if Q1 < threshold
            continue
        end
        [score2, Q2] = nnRatioN(R(r,:), u, datasetC);
        if Q2 < threshold
            continue
        end
        if score1/Q1 > max
            max = score1/Q1;
            maxFeat = R(r,:);
        end
        if score2/Q2 > max
            max = score2/Q2;
            maxFeat = R(r,:);
        end
    end
end
v = maxFeat;