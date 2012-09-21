function correct_image(path, name)
    close all;
%     path = '~/Desktop/query/project/src/tutorial/testpatch/'
%     name = 'test2.jpg'
    fullpath = fullfile(path, name);
    I_O = imread(fullpath);
    I = rgb2gray(I_O);
    I = medfilt2(I, [4 4]);
    [BW3, thresh] = edge(I,'canny');

    SE = strel('disk', 1);
    SE2 = strel('rectangle', [3 2]);
    BW3 = imdilate(BW3,SE);
    s = getBoundingBoxes(BW3);
%     dispBoundingBoxes(BW3, s);
    
    skewangle = getHoughAngle(BW3, s);
    
    BW3 = imrotate(BW3,-skewangle, 'bicubic');
    s = getBoundingBoxes(BW3);
%     dispBoundingBoxes(BW3, s);
    
    slantangle = getSlantAngle(BW3, s);
    tform = maketform('affine',[1 0 0; slantangle 1 0; 0 0 1]);
    BW3 = imtransform(BW3,tform, 'bicubic');
    s = getBoundingBoxes(BW3);
%     dispBoundingBoxes(BW3, s);
    
    ratioarray = zeros(size(s));
    for i = 1:numel(s)
        ratioarray(i) = s(i).BoundingBox(4)/s(i).BoundingBox(3);
    end
    
    myStdev = std(ratioarray);
    myMean = mean(ratioarray);
    
    avgRatio = 0;
    counter = 0;
    
    for i = 1:numel(ratioarray)
        if abs((ratioarray(i) - myMean)/myStdev) < 1.281
            avgRatio = avgRatio + ratioarray(i);
            counter = counter + 1;
        end
    end
    
    avgRatio = avgRatio/counter;
    disp(strcat('Average Bounding Box Ratio: ', num2str(avgRatio)));
    I_O = imrotate(I_O,-skewangle);
    I_O = imtransform(I_O, tform);
    
    if avgRatio == 0
        avgRatio = 1.4
    end
    
    if avgRatio > 1.4
        I_O = imresize(I_O, 'Scale', [1 avgRatio/1.4]);
    end
    
    I_O = rgb2gray(I_O);
    t = graythresh(I_O);
%     I_O = adaptivethreshold(I_O, 11, 0.03, 0);
%     disp(t);
    I_O = im2bw(I_O , t);
    
%     figure, imshow(I_O);
    writePath = fullfile(path, strcat('corrected', name));
    imwrite(I_O, writePath,'jpg');
    
end

function avgAng = getSlantAngle(BW3, s)
    ccwrotate45 = [cosd(45) -sind(45); sind(45) cosd(45)];
    cwrotate45 = [cosd(45) -sind(45); sind(45) cosd(45)];
    avgAng = 0;
    for i = 1:numel(s)
        minx = bitmax;
        miny = bitmax;
        maxx = -bitmax;
        maxy = -bitmax;
        for j = 1:size(s(i).PixelList, 1)
            temp = s(i).PixelList(j,:) * ccwrotate45;
            if temp(1) < minx
                minx = temp(1);
            end
            if temp(2) < miny
                miny = temp(2);
            end
            if temp(1) > maxx
                maxx = temp(1);
            end
            if temp(2) > maxy
                maxy = temp(2);
            end
            
        end

        top = [minx maxy] * cwrotate45; 
        bottom = [maxx miny] * cwrotate45;
        
        myB = top(1) - bottom(1);
        myA = s(i).BoundingBox(4);
        avgAng = avgAng + atan(myB/myA); 
%         disp(atan(myB/myA));
    end
    if numel(s) ~= 0
        avgAng = avgAng/numel(s);
    else
        avgAng = 0;
    end
    disp(strcat('Slant Angle: ', num2str(avgAng)));
end

function angle = getHoughAngle(BW3, s)
    houghmat = zeros(size(BW3));
    for i = 1:numel(s)
        x = s(i).BoundingBox;
        houghmat(round(x(1)+ x(3)/2), round(x(2) + x(4)/2)) = 1;
    end
    [H T R] = hough(houghmat, 'RhoResolution',3, 'Theta',-45:.1:45);
    
%     figure, imshow(imadjust(mat2gray(H)), [], 'XData',T, 'YData',R, ...
%       'InitialMagnification','fit')
%     xlabel('\theta (degrees)'), ylabel('\rho')
%     axis on, axis normal, hold on
%     colormap(hot), colorbar

    P  = houghpeaks(H, 50, 'threshold',ceil(max(H(:))));
%     plot(T(P(:,2)), R(P(:,1)), 'gs', 'LineWidth',2);
    
    
    mySize = size(H);
    offset = mySize(2)/2;
    temp = 0;
    for i = 1:size(P,1)
        temp = temp + P(i,2);
    end
    temp = temp/size(P,1);
    if numel(s) <= 1
        angle = 0;
    else
        angle = (temp-offset)/offset * 45;
    end
    disp(strcat('Skew Angle: ', num2str(angle)));
end

function dispBoundingBoxes(BW3, s)
    L = bwlabel(BW3, 8);
    RGB = label2rgb(L, 'jet', 'w', 'shuffle'); 
    figure, imshow(RGB);

    hold on
    for i = 1:numel(s)
        x = s(i).BoundingBox;
        
        hnd1=text('Position',[x(1)+ x(3)/2 x(2) + x(4)/2], 'String', i);
        set(hnd1,'FontSize',12);
        rectangle('Position',s(i).BoundingBox());
    end
    hold off
end

function s = getBoundingBoxes(BW3)
    s  = regionprops(BW3, 'Centroid', 'BoundingBox', 'Area', 'FilledArea', 'PixelList', 'Extrema');
    totalarea = size(BW3, 1) * size(BW3, 2);
    i = 1;
    while i <= numel(s)
        if (s(i).FilledArea/totalarea) < 0.000208
            s(i) = [];

            continue;
        end
        w = s(i).BoundingBox(3);
        h = s(i).BoundingBox(4);
        if (w/h > 1.2) || (h/w > 8) || (s(i).FilledArea/(w*h) < 0.25) || (s(i).FilledArea/(w*h) > 0.9)
            s(i) = [];
            continue;
        end
        i = i + 1;
    end

    j = 1;
    while j <= numel(s)
        if line_test(s(j).PixelList, 8) < 0.625
            s(j) = [];
            continue;
        end
        j = j + 1;
    end
end

function ratio = line_test(pL, numlines)
    for i = 1:size(pL, 1)
        temp = pL(i, 2);
        pL(i, 2) = pL(i, 1);
        pL(i, 1) = temp;
    end
    pL = sortrows(pL, [1, 2]);
    %disp(pL);
    low = (pL(1, 1));
    high = (pL(end, 1));
    int = floor((high - low) / (numlines + 1));
    linearray = zeros(numlines, 1);
    temp = low;
    for i = 1:numlines
        temp = temp + int;
        linearray(i) = temp;
    end
    %disp(linearray);
    numlinespass = 0;
    firstofline = 1;
    numintersects = -1;
    prevNumber = -1;
    for i = 1:size(pL, 1)
        if any(linearray == pL(i, 1))
            if firstofline == 1
                firstofline = 0;
                numintersects = 1;
                prevNumber = pL(i, 2);
            else
                if prevNumber ~= pL(i, 2) - 1
                    numintersects = numintersects + 1;
                end
                prevNumber = pL(i, 2);
            end
        else
            if numintersects ~= -1 && numintersects > 1 && numintersects < 5
                numlinespass = numlinespass + 1;
            end
            numintersects = -1;
            firstofline = 1;
            prevNumber = -1;
        end
    end
    ratio = numlinespass/numlines;
end

function bw=adaptivethreshold(IM,ws,C,tm)
%ADAPTIVETHRESHOLD An adaptive thresholding algorithm that seperates the
%foreground from the background with nonuniform illumination.
%  bw=adaptivethreshold(IM,ws,C) outputs a binary image bw with the local
%   threshold mean-C or median-C to the image IM.
%  ws is the local window size.
%  tm is 0 or 1, a switch between mean and median. tm=0 mean(default); tm=1 median.
%
%  Contributed by Guanglei Xiong (xgl99@mails.tsinghua.edu.cn)
%  at Tsinghua University, Beijing, China.
%
%  For more information, please see
%  http://homepages.inf.ed.ac.uk/rbf/HIPR2/adpthrsh.htm

    if (nargin<3)
        error('You must provide the image IM, the window size ws, and C.');
    elseif (nargin==3)
        tm=0;
    elseif (tm~=0 && tm~=1)
        error('tm must be 0 or 1.');
    end

    IM=mat2gray(IM);

    if tm==0
        mIM=imfilter(IM,fspecial('average',ws),'replicate');
    else
        mIM=medfilt2(IM,[ws ws]);
    end
    sIM=mIM-IM-C;
    bw=im2bw(sIM,0);
    bw=imcomplement(bw);
end

function [ optimalThreshold, J ] = kittlerMinimimErrorThresholding( img )
%KITTLERMINIMIMERRORTHRESHOLDING Compute an optimal image threshold.
%   Computes the Minimum Error Threshold as described in
%   
%   'J. Kittler and J. Illingworth, "Minimum Error Thresholding," Pattern
%   Recognition 19, 41-47 (1986)'.
%   
%   The image 'img' is expected to have integer values from 0 to 255.
%   'optimalThreshold' holds the found threshold. 'J' holds the values of
%   the criterion function.

%Initialize the criterion function
J = Inf * ones(255, 1);

%Compute the relative histogram
histogram = double(histc(img(:), 0:255)) / size(img(:), 1);

%Walk through every possible threshold. However, T is interpreted
%differently than in the paper. It is interpreted as the lower boundary of
%the second class of pixels rather than the upper boundary of the first
%class. That is, an intensity of value T is treated as being in the same
%class as higher intensities rather than lower intensities.
for T = 1:255

    %Split the hostogram at the threshold T.
    histogram1 = histogram(1:T);
    histogram2 = histogram((T+1):end);

    %Compute the number of pixels in the two classes.
    P1 = sum(histogram1);
    P2 = sum(histogram2);

    %Only continue if both classes contain at least one pixel.
    if (P1 > 0) && (P2 > 0)

        %Compute the standard deviations of the classes.
        mean1 = sum(histogram1 .* (1:T)') / P1;
        mean2 = sum(histogram2 .* (1:(256-T))') / P2;
        sigma1 = sqrt(sum(histogram1 .* (((1:T)' - mean1) .^2) ) / P1);
        sigma2 = sqrt(sum(histogram2 .* (((1:(256-T))' - mean2) .^2) ) / P2);

        %Only compute the criterion function if both classes contain at
        %least two intensity values.
        if (sigma1 > 0) && (sigma2 > 0)

            %Compute the criterion function.
            J(T) = 1 + 2 * (P1 * log(sigma1) + P2 * log(sigma2)) ...
                     - 2 * (P1 * log(P1) + P2 * log(P2));

        end
    end

end

%Find the minimum of J.
[~, optimalThreshold] = min(J);
optimalThreshold = optimalThreshold - 0.5;
end
