function [ outputn outputc outputa ] = rank( query, direc )
%RANK Summary of this function goes here
%   Detailed explanation goes here

files = dir(strcat(direc,'*.pgm'));
outputc=zeros(1,numel(files));
outputa=zeros(1,numel(files));
counts= zeros(1,numel(files));
avgs = zeros(1,numel(files));

[im, des, loc] = sift(query);
for i = 1:numel(files)
    [a,b] = match(des,strcat(direc,files(i).name));
    counts(i)=a;
    avgs(i)=b;
end
[b,idx]=sort(counts);
for i = 1:numel(idx)
    outputn(i)=cellstr(files(idx(i)).name);
    outputc(i)=counts(idx(i));
    outputa(i)=avgs(idx(i));
end