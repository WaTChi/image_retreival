function PP_LocErr(lrnum,errFlag)

nerrs=2;
if ~exist('errFlag','var')
    errFlag = 0;
end
if ~exist('lrnum','var')
    lrnum = 0;
end
fig_offset = 10*lrnum;
ncomp = 1;
if lrnum == 0
    file = 'Z:\ah\lastrun\reproj.txt';
    ncomp = 2;
else
    file = ['Z:\ah\lastrun',num2str(lrnum),'\reproj.txt'];
end

% Remove specified queries
fid = fopen(file);
names = textscan(fid,'%s%f%f%f%f%f%f%f%f%f%f%f%f');
names = names{1};
idx = [find(~cellfun('isempty',strfind(names,'2011-04-04_14-56-08_926'))), ...
       find(~cellfun('isempty',strfind(names,'2011-04-04_15-06-58_888'))), ...
       find(~cellfun('isempty',strfind(names,'2011-04-04_15-19-11_784'))), ...
       find(~cellfun('isempty',strfind(names,'2011-04-04_15-34-01_368')))];
fclose(fid);

results = importdata(file);
results(idx,:) = [];

dx = 2; xmax = 30;
r = results(:,2:ncomp+1);
r(r>100) = 100;

for k=1:ncomp
    figure(fig_offset+k)
    bins = dx/2:dx:xmax+dx/2;
    err = r(:,k);
    count = cumsum(hist(err,bins));
    count = 100 * count / count(end);
    hold off
    bar(bins(1:end-1),count(1:end-1),'g')
    hold on
    count(bins<10) = 0;
    bar(bins(1:end-1),count(1:end-1),'r')
    xlabel('Location error (meters)')
    ylabel('Percentage of queries')
    title('Percentage of queries localized within x meters')
    xlim([0,xmax])
    ylim([0,100])
end

% Constrained homography analysis
% Remove extreme errors for better correlation
if errFlag
    loc = results(:,2);
    tra = results(:,2+nerrs);
    pna = results(:,3+nerrs);
    loc = (loc<5) + (loc<10) + (loc<20);
    pna = (pna<5) + (pna<10) + (pna<20);
    tra = (tra<5) + (tra<10) + (tra<20);
else
    mask = results(:,2) < 50;
    loc = results(mask,2);
    tra = results(mask,2+nerrs);
    pna = results(mask,3+nerrs);
    fprintf('\n')
    disp(['Median location error = ',num2str(round(median(loc))),' meters'])
    disp(['Fraction under 5 meters: ',num2str(100*sum(loc<5)/length(loc),2),'%'])
    disp(['Fraction under 10 meters: ',num2str(100*sum(loc<10)/length(loc),2),'%'])
    disp(['Median translation angle error = ',num2str(median(tra)),' degrees'])
    disp(['Fraction under 10 degrees: ',num2str(100*sum(tra<10)/length(loc),2),'%'])
    disp(['Median normal angle error = ',num2str(median(pna)),' degrees'])
    disp(['Fraction under 10 degrees: ',num2str(100*sum(pna<10)/length(loc),2),'%'])
end
CC = corrcoef([loc tra pna]);
fprintf('\n')
disp('Correlation coefficients: Loc TrA PnA')

disp(CC)