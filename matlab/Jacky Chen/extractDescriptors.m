function [featureDescriptors] = extractDescriptors(fileName)

%   extractDescriptors(fileName) returns an array of SIFT descriptors from
%   .sift files.
%
%   Input:
%       - fileName: name of the .sift file
%
%   Output:
%       - featureDescriptors: matrix with SIFT descriptors in rows
%
SIFT_LENGTH = 128;
ROW_LENGTH = 20;

fid = fopen(fileName, 'r');
numSift = textscan(fid,'%f',1);
featureDescriptors = zeros(numSift{1}, SIFT_LENGTH);
textscan(fid,'%f',1);

for i = 1:numSift{1}
    descriptor = zeros(1, 0);
    textscan(fid,'%f',4);
    for j = 1:floor(SIFT_LENGTH/ROW_LENGTH)
        tmp = textscan(fid,'%f',ROW_LENGTH);
        descriptor = [descriptor tmp{1}'];
    end
    tmp = textscan(fid,'%f',rem(SIFT_LENGTH, ROW_LENGTH));
    descriptor = [descriptor tmp{1}'];
    featureDescriptors(i,:) = descriptor;
end
fclose(fid);