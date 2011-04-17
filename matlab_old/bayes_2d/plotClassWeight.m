function plotClassWeight(classifier,class1,class2,min_samp,fignum)

% plotClassWeight(classifier,plotmode,class1,class2,fignum)
% 
% M = number of features
% C = number of classes
% 
% DESCRIPTION
%   This function plots the relative weighing of the two classes specified
%   (the first two if left unspecified) for a range of feature values and
%   for each feature individually.
% 
% INPUTS
%   classifier: Structure containing the classifier information
%   class1:     The first class to compare (optional; default = 1)
%   class2:     The second class to compare (optional; default = class1+1)
%       - A weight > 1 in the plot means that the given feature value
%         class1, while a weight < 1 in the plot favors class 2.
%   fignum:     Optional figure start. Plots in figure(fignum+1:fignum+M)

% Size parameters
C = length(classifier.nsamps);
M = length(classifier.bins);

% Set classes
if nargin < 2
    cls = [1,2];
elseif nargin < 3
    cls = [class1,mod(class1+1,C)];
else
    cls = [class1,class2];
end

% Set min_samp
if nargin < 4
    min_samp = 20;
end
total = sum(classifier.nsamps);

% Set plot flag
plotflag = (nargin<5);

% Load necessary information
nsamps = classifier.nsamps;
sp = classifier.spacing;
bins = classifier.bins;
B = size(bins); B = B(1:2);
x = sp(1)/2 : sp(1) : (B(1)-1/2)*sp(1);
y = sp(2)/2 : sp(2) : (B(2)-1/2)*sp(2);
prob1 = zeros(B);
prob2 = zeros(B);

% Get relative weighing via classify function bin expansion method
for j=1:B(1)
    
    fprintf([num2str(j),'/',num2str(B(1)),'\n'])
    
    for k=1:B(2)
        
        % Expand bin region to satisfy min_samp
        sb1 = j; lb1 = j;
        sb2 = k; lb2 = k;
        frct1 = (lb1-sb1+1)/B(1);
        frct2 = (lb2-sb2+1)/B(2);
        count = reshape( sum( sum( bins(sb1:lb1,sb2:lb2,:) , 1 ) , 2 ) , [1,2] );
        while sum(count) < min_samp
            frct1 = (lb1-sb1+1)/B(1);
            frct2 = (lb2-sb2+1)/B(2);
            if frct2 > frct1
                sb1 = max(1,sb1-1);
                lb1 = min(B(1),lb1+1);
            else
                sb2 = max(1,sb2-1);
                lb2 = min(B(2),lb2+1);
            end
            count = reshape( sum( sum( bins(sb1:lb1,sb2:lb2,:) , 1 ) , 2 ) , [1,2] );
        end
        
        prob1(j,k) = count(1) / nsamps(1);
        prob2(j,k) = count(2) / nsamps(2);
        
    end
    
end

weight = prob1 ./ prob2;

% Create plots and label them
mx = max(max(max(prob1,prob2)));
if plotflag
    figure
else
    figure(fignum+1)
end
title('Feature distribution for class 1')
xlabel('Feature 1 values')
ylabel('Feature 2 values')
imagesc(x,y,prob1)
colorbar
caxis([0,mx])
axis xy

if plotflag
    figure
else
    figure(fignum+2)
end
title('Feature distribution for class 2')
xlabel('Feature 1 values')
ylabel('Feature 2 values')
imagesc(x,y,prob2)
colorbar
caxis([0,mx])
axis xy


if plotflag
    figure
else
    figure(fignum+3)
end
title('Relative weight between class 1 and 2')
xlabel('Feature 1 values')
ylabel('Feature 2 values')
imagesc(x,y,log(weight))
colorbar
axis xy


