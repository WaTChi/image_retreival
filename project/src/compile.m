%% 
%  Copyright (c) 2011  Chen Feng (cforrest[at]umich[dot]edu)
%   and the University of Michigan
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.

%% compile
mex -O -output lsd lsd_matlab.c lsd.c
%% test
clear all; close all; clc;
chairs = imread('SnowCity.jpeg');
imshow(chairs); hold on;
tic
lines = lsd(double(chairs));
%%lines = LineSegmentDetection(double(chairs), 0.8, 0.6, 2.0, 22.5, 0.0, 0.7, 1024, 255.0, 0.0);
%%disp(lines);
t = toc;
disp(['[lsd] ',num2str(t),' seconds elapsed.']);
nl = size(lines,1);
disp(nl);
for i=1:nl
    plot(lines(i,1:2:4),lines(i,2:2:4),'r-');
end