function [] = findBestMatches(dataset, range, threshold)

fileName = strcat('unsortedMatches', dataset);
fileName = strcat(fileName, 'c.txt');

fid = fopen(fileName, 'r');
match = textscan(fid, '%s %s', 'delimiter', '-');
i = match{1};
j = match{2};

newFile = strcat('bm', dataset);
newFile = strcat(newFile, '.txt');
nid = fopen(newFile, 'wt');
datasetC = [dataset, 'r'];
for k = 1:size(i)
    fprintf('Processing %d out of %d pairs\n', k, size(i));
    if k == 5
        fprintf('debug');
    end
    [u, v] = findBestMatch(i(k), j(k), range, datasetC, threshold);
    fprintf(nid, '%s-%s', u, v);
    fprintf(nid, '\n');
end
fclose(fid);