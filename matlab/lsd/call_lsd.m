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

function call_lsd(imgpath,outpath)

%% compile
mex -O -output lsd lsd_matlab.c lsd.c
%% test
tic
img = rgb2gray(imread(imgpath));
lines = lsd(double(img));
save(outpath,'lines','-ASCII')
t = toc;
disp(['[lsd] ',num2str(t),' seconds elapsed.']);
% nl = size(lines,1);
% imshow(chairs); hold on;
% disp(nl);
% for i=1:nl
%     plot(lines(i,1:2:4),lines(i,2:2:4),'r-');
% end