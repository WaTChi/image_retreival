function pose_analysis(runnum,setnum,gps)

% function pose_analysis(runnum,setnum,gps)
% runnum - run number (0 = last, >0 is number appended to folder)
% setnum - set number (10+i for set i in Oakland, i for set i in Berkeley)
% gps    - flag for displaying gps error

if ~exist('runnum','var')
    runnum = 0;
end
if ~exist('setnum','var')
    setnum = 11;
end
if ~exist('gps','var')
    gps = false;
end
if ~exist('yaw','var')
    yaw = false;
end
fig_offset = 100*(setnum<10) + 10*runnum;
if setnum < 10
    if runnum == 0
        file = 'Z:\ah\pose_runs\berkeley\pose_results.txt';
%         file = '/media/DATAPART2/ah/pose_runs/berkeley/pose_results.txt';
%         gps = true; yaw = true;
    else
        file = ['Z:\ah\pose_runs\berkeley',num2str(runnum),'\pose_results.txt'];
%         file =
%         ['/media/DATAPART2/ah/pose_runs/berkeley',num2str(runnum),'/pose_results.txt'];
    end
else
    if runnum == 0
        file = 'Z:\ah\pose_runs\oakland\pose_results.txt';
%         file = '/media/DATAPART2/ah/pose_runs/oakland/pose_results.txt';
%         gps = true; yaw = true;
    else
        file = ['Z:\ah\pose_runs\oakland',num2str(runnum),'\pose_results.txt'];
%         file =
%         ['/media/DATAPART2/ah/pose_runs/oakland',num2str(runnum),'/pose_results.txt'];
    end
end
% % Remove specified queries
% fid = fopen(file);
% names = textscan(fid,'%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
% names = names{1};
% idx = [find(~cellfun('isempty',strfind(names,'2011-10-28_13-00-41_047'))), ...
%        find(~cellfun('isempty',strfind(names,'2011-04-04_15-17-45_763'))), ...
%        find(~cellfun('isempty',strfind(names,'2011-04-04_15-22-31_544'))), ...
%        find(~cellfun('isempty',strfind(names,'2011-10-28_12-23-24_368'))), ...
%        find(~cellfun('isempty',strfind(names,'2011-10-28_12-13-41_589'))), ...
%        find(~cellfun('isempty',strfind(names,'2011-10-28_12-12-29_132'))), ...
%        find(~cellfun('isempty',strfind(names,'2011-10-28_12-12-20_282'))), ...
%        find(~cellfun('isempty',strfind(names,'2011-10-28_12-11-42_811'))), ...
%        find(~cellfun('isempty',strfind(names,'2011-10-28_12-09-45_828'))), ...
%        find(~cellfun('isempty',strfind(names,'2011-10-28_12-09-35_578'))), ...
%        find(~cellfun('isempty',strfind(names,'2011-10-28_11-59-53_593'))), ...
%        find(~cellfun('isempty',strfind(names,'2011-10-28_12-01-35_821'))), ...
%        find(~cellfun('isempty',strfind(names,'2011-04-04_14-56-08_926'))), ...
%        find(~cellfun('isempty',strfind(names,'2011-04-04_15-06-58_888'))), ...
%        find(~cellfun('isempty',strfind(names,'2011-04-04_15-19-11_784'))), ...
%        find(~cellfun('isempty',strfind(names,'2011-04-04_15-32-55_864'))), ...
%        find(~cellfun('isempty',strfind(names,'2011-10-28_12-45-09_775'))), ...
%        find(~cellfun('isempty',strfind(names,'2011-04-04_15-34-01_368')))];
% fclose(fid);

% Import pose results
disp(file)
results = importdata(file);
% results(idx,:) = [];

% Pose analysis parameters
dx = 1; xtick = 5; xmax = 25;
bins = dx/2:dx:xmax+dx/2;
if setnum == 11
    mincount = 94; % number of queries to analyze for query 1 Oakland
else % if setnum == 5
    mincount = 79; % number of queries to analyze for query 5 Berkeley
end
mincount = 0;

% Plot location error
err = results(:,2); expv = err(err<25);
err(err>100) = 100;
figure(fig_offset+1), clf, hold on
count = cumsum(hist(err,bins));
count = 100 * count / max(mincount,count(end));
bar(bins(1:end-1),count(1:end-1),'r')
bar(bins(bins<10),count(bins<10),'g')
xlabel('Location error (meters)')
ylabel('Percentage of queries')
if setnum < 10, title('Percent of Berkeley queries localized within x meters (homography)')
else            title('Percent of Oakland queries localized within x meters (homography)')
end
xlim([0,xmax])
ylim([0,100])
set(gca,'XTick',[0:xtick:xmax])

expv = sqrt(mean(expv.^2))

% Plot GPS error if runnum = 0 or gps flag is set
if gps
    err = results(:,3);
    err(err>100) = 100;
    figure(fig_offset+9), clf, hold on
    count = cumsum(hist(err,bins));
    count = 100 * count / max(mincount,count(end));
    bar(bins(1:end-1),count(1:end-1),'r')
    bar(bins(bins<10),count(bins<10),'g')
    xlabel('Location error (meters)')
    ylabel('Percentage of queries')
    if setnum < 10, title('Percent of Berkeley queries localized within x meters (GPS)')
    else        	title('Percent of Oakland queries localized within x meters (GPS)')
    end
    xlim([0,xmax])
    ylim([0,100])
    set(gca,'XTick',[0:xtick:xmax])
end
% 
% if yaw
%     err = results(:,end);
%     figure(fig_offset+8)
%     yaw_bins = 1:2:31;
%     count = cumsum(hist(err,yaw_bins));
%     count = 100 * count / max(mincount,count(end));
%     hold off
%     bar(yaw_bins(1:end-1),count(1:end-1),'g')
%     hold on
%     count(yaw_bins<10) = 0;
%     bar(yaw_bins(1:end-1),count(1:end-1),'r')
%     xlabel('Cell phone yaw error (degrees)')
%     ylabel('Percentage of queries')
%     title('Percentage of queries with yaw error within x degrees')
%     xlim([0,30])
%     ylim([0,100])
%     set(gca,'XTick',[0:5:30])
%     hold off
% end
% 
% % Constrained homography analysis
% loc = results(:,2);
% close = loc < 10;
% tra = results(:,5);
% pna = results(:,6);
% disp(['Median location error = ',num2str(round(median(loc))),' meters'])
% disp(['Percentage under 5 meters: ',num2str(100*sum(loc<5)/length(loc),2),'%'])
% disp(['Percentage under 10 meters: ',num2str(100*sum(loc<10)/length(loc),2),'%'])
% disp(['Median translation angle error = ',num2str(median(tra)),' degrees'])
% disp(['Percentage under 10 degrees: ',num2str(100*sum(tra<10)/length(loc),2),'%'])
% disp(['Median normal angle error = ',num2str(median(pna)),' degrees'])
% disp(['Percentage under 10 degrees: ',num2str(100*sum(pna<10)/length(loc),2),'%'])
% CC = corrcoef([close tra pna]);
% fprintf('\n')
% disp('Correlation coefficients: Loc TrA PnA')
% disp(CC)
% 
% % Evaluate different indicators for pose estimation performance
% numq = length(close);
% ind1 = results(:,7)./results(:,8);
% ind2 = 1./results(:,10);
% CC = corrcoef([close ind1 ind2]);
% fprintf('\n')
% disp('Correlation coefficients of indicators: Loc Ind1 Ind2')
% disp(CC)
% 
% % Set up ROC curve analysis
% num_thresholds = 10000;
% P = sum(close); N = numq - P;
% zto = [0,1];
% 
% % Plot indicator 1 ROC curve: numi / numd
% minind = min(ind1); maxind = max(ind1);
% dthreshold = (maxind - minind) / (num_thresholds-3);
% thresholds = repmat( minind-dthreshold : dthreshold : maxind+dthreshold , [numq,1] );
% ind1mat = repmat( ind1 , [1,num_thresholds] ); closemat = repmat( close , [1,num_thresholds] );
% TPR = sum( closemat & (ind1mat>thresholds) , 1 ) / P;
% FPR = sum( ~closemat & (ind1mat>thresholds) , 1 ) / N;
% figure(fig_offset+2)
% hold off
% plot(FPR,TPR)
% hold on
% plot(zto,zto,'r-.')
% xlim([0,1])
% ylim([0,1])
% xlabel('False positive rate')
% ylabel('True positive rate')
% title('ROC Curve for Indicator = numi / numq')
% hold off
% 
% % Plot indicator 1 ROC curve: numi / numd
% minind = min(ind2); maxind = max(ind2);
% dthreshold = (maxind - minind) / (num_thresholds-3);
% thresholds = repmat( minind-dthreshold : dthreshold : maxind+dthreshold , [numq,1] );
% ind2mat = repmat( ind2 , [1,num_thresholds] ); closemat = repmat( close , [1,num_thresholds] );
% TPR = sum( closemat & (ind2mat>thresholds) , 1 ) / P;
% FPR = sum( ~closemat & (ind2mat>thresholds) , 1 ) / N;
% figure(fig_offset+3)
% hold off
% plot(FPR,TPR)
% hold on
% plot(zto,zto,'r-.')
% xlim([0,1])
% ylim([0,1])
% xlabel('False positive rate')
% ylabel('True positive rate')
% title('ROC Curve for Indicator = Ind1 / Reproj Error')
% hold off
