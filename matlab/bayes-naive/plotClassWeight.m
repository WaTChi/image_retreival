function plotClassWeight(classifier,class1,class2,fignum)

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

% Set plot flag
plotflag = (nargin<4);

% Iterate through features
for f=1:M
    
    % Create plots and label it
    if plotflag
        figure
    else
        figure(fignum+f)
    end
    subplot(211)
    title(['Feature ',num2str(f),' distribution'])
    xlabel(['Feature ',num2str(f),' values'])
    ylabel('Probability density')
    subplot(212)
    title(['Feature ',num2str(f),' weighing'])
    xlabel(['Feature ',num2str(f),' values'])
    ylabel(['Relative weiging : Class ',num2str(cls(1)), ...
        ' / Class ',num2str(cls(2))])
    
    % Load necessary information
    nsamps = classifier.nsamps;
    spacing = classifier.spacing(f);
    bins = classifier.bins{f};
    prob1 = (1/spacing) * bins(:,cls(1)) / nsamps(cls(1));
    prob2 = (1/spacing) * bins(:,cls(2)) / nsamps(cls(2));
    B = length(prob1);
    x = spacing/2 : spacing : (B-1/2)*spacing;
    
    % Get relative weighing via classify function bin expansion method
    min_samp = 100;
    weight = zeros(B,1);
    for b=1:B
        
        % Find the size of the bin necessary to satisfy min_samp
        sb = b; lb = b; % small and large bins
        count = sum( bins(sb:lb,:) , 1 );
        while sum(count) < min_samp
            sb = max(1,sb-1);
            lb = min(B,lb+1);
            count = sum( bins(sb:lb,:) , 1 );
        end
        weight(b) = count(cls(1)) / count(cls(2));
        
    end
        
    
    % Plot curves
    subplot(211)
    hold off
    plot(x,prob1,'b')
    hold on
    plot(x,prob2,'r')
    subplot(212)
    semilogy(x,weight,'k')
    
end