function [num] = numOfFeatures(fileName)

fid = fopen(fileName, 'r');
numSift = textscan(fid,'%f',1);
num = numSift{1};
fclose(fid);