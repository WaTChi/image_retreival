function plotClassWeight(classifier,plotmode,class1,class2,fignum)

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
%   plotmode:   String for plot mode:
%       - 'xy' : Linear in x and y
%       - 'Xy' : Logarithmic in x, linear in y
%       - 'xY' : Linear in x, logarithmic in y
%       - 'XY' : Logarithmic in x and y
%   class1:     The first class to compare (optional; default = 1)
%   class2:     The second class to compare (optional; default = class1+1)
%       - A weight > 1 in the plot means that the given feature value
%         class1, while a weight < 1 in the plot favors class 2.
%   fignum:     Optional figure start. Plots in figure(fignum+1:fignum+M)

% Size parameters
[C,M] = size(classifier.means);
npts = 100;

% Set classes
if nargin < 3
    cls = [1,2];
elseif nargin < 4
    cls = [class1,mod(class1+1,C)];
else
    cls = [class1,class2];
end

% Set plot flag
plotflag = (nargin<5);

% Iterate through features
for f=1:M
    
    % Create plot and label it
    if plotflag
        figure
    else
        figure(fignum+f)
    end
    title(['Feature ',num2str(f),' weighing'])
    xlabel(['Feature ',num2str(f),' values'])
    ylabel(['Relative weiging : Class ',num2str(cls(1)), ...
        ' / Class ',num2str(cls(2))])
    
    if strcmp(classifier.dists{f},'Gamma')
        
        means = classifier.means(cls,f);
        vars = classifier.vars(cls,f);
        k = means.^2 ./ vars;
        t = vars ./ means;
        
        if strcmp(plotmode(1),'x') % linear in x
%             maxX = max(means+3*sqrt(vars));
%             x = 0:maxX/npts:maxX;
            x = 0:1/npts:1;
        else % logarithmic in x
            extent = max(log(k));
            minX = min(log(means)-2*extent);
            maxX = max(log(means)+extent);
            x = exp( minX:(maxX-minX)/npts:maxX );
        end
        
        log1 = (k(1)-1)*log(x) - x/t(1) - log(gamma(k(1))) - k(1)*log(t(1));
        log2 = (k(2)-1)*log(x) - x/t(2) - log(gamma(k(2))) - k(1)*log(t(2));
        y = exp(log1-log2);
        
    else % if strcmp(classifier.dists{f},'Normal')
        
        means = classifier.means(cls,f);
        vars = classifier.vars(cls,f);
        
        if strcmp(plotmode(1),'x') % linear in x
            minX = min(means-2*sqrt(vars));
            maxX = max(means+2*sqrt(vars));
            x = minX:(maxX-minX)/npts:maxX;
        else % logarithmic in x
            extent = max(log(means.^2./vars));
            minX = min(log(means)-extent);
            maxX = max(log(means)+extent);
            x = exp( minX:(maxX-minX)/npts:maxX );
        end
        
        log1 = (-1/2)*log(2*pi*vars(1)) - (x-means(1)).^2/(2*vars(1));
        log2 = (-1/2)*log(2*pi*vars(2)) - (x-means(2)).^2/(2*vars(2));
        y = exp(log1-log2);
        
    end
    
    if strcmp(plotmode,'xy') % linear in x and y
        plot(x,y)
    elseif strcmp(plotmode,'Xy') % logarithmic in x and linear in y
        semilogx(x,y)
    elseif strcmp(plotmode,'xY') % linear in x and logarithmic in y
        semilogy(x,y)
    else % logrithmic in x and y
        loglog(x,y)
    end
    
end