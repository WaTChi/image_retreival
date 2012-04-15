function prehom_analysis(runflag,setnum)

if nargin < 2
    setnum = 11;
end
if nargin < 1
    runflag = 0;
end
oakland = setnum > 10;

% import data
if setnum < 10
    if runflag == 0
        file = 'Z:\ah\pose_runs\berkeley\extras.txt';
%         file = '/media/DATAPART2/ah/pose_runs/berkeley/extras.txt';
    else
        file = ['Z:\ah\pose_runs\berkeley',num2str(runflag),'\extras.txt'];
%         file = ['/media/DATAPART2/ah/pose_runs/berkeley',num2str(runflag),'/extras.txt'];
    end
else
    if runflag == 0
        file = 'Z:\ah\pose_runs\oakland\extras.txt';
%         file = '/media/DATAPART2/ah/pose_runs/oakland/extras.txt';
    else
        file = ['Z:\ah\pose_runs\oakland',num2str(runflag),'\extras.txt'];
%         file = ['/media/DATAPART2/ah/pose_runs/oakland',num2str(runflag),'/extras.txt'];
    end
end
prehom = importdata(file);

% add rows that went missing
idx = 1;
while idx+7 <= size(prehom,1)
    if prehom(idx+7,1) < 1000 && prehom(idx+4,1) < 1000
        prehom = [ prehom(1:idx,:);
                   NaN*zeros(6,7);
                   prehom(idx+1:end,:) ] ;
    elseif prehom(idx+7,1) < 1000 && prehom(idx+3,1) < 0
        prehom = [ prehom(1:idx+3,:);
                   NaN*zeros(3,7);
                   prehom(idx+4:end,:) ] ;
    elseif prehom(idx+7,1) < 1000 % && prehom(idx+3,1) > 0
        prehom = [ prehom(1:idx,:);
                   NaN*zeros(3,7);
                   prehom(idx+1:end,:) ] ;
    end
    idx = idx+7;
end
if size(prehom,1) < idx+6
    prehom = [ prehom(1:idx,:);
               NaN*zeros(3,7);
               prehom(idx+1:end,:) ] ;
end

% Remove specified queries
times = [ 130041 , 122324 , 121341 , 121229 , 121220 , ...
          121142 , 120945 , 120935 , 115953 , 120135 , ...
          145608 , 150658 , 151911 , 153255 , 153401 , ...
          151745 , 152231 ];
for i=1:length(times)
    rmidx = find(prehom(1:7:end,1)==times(i));
    prehom(7*rmidx-6:7*rmidx,:) = [];
end

% parse prehomography data
tnorms = prehom(1:7:end,2);
tyaws  = prehom(1:7:end,3);
cyaws  = prehom(1:7:end,4);
vyaws  = prehom(1:7:end,5);
vconfs = prehom(1:7:end,6);
vyerrs = prehom(1:7:end,7);
vnconf = prehom(2:7:end,:);
vnorms = prehom(3:7:end,:);
dnorms = prehom(6:7:end,:);
dnerr  = prehom(7:7:end,:);

% % Analyze yaw estimate
% numq = length(tnorms);
% close_angle = 10;
% thresholds = 0:.01:1;
% truth_idx = vyerrs <= close_angle;
% P = sum(truth_idx); N = numq - sum(truth_idx);
% tprs = zeros(numq,1); fprs = zeros(numq,1);
% for i=1:length(thresholds)
%     estim_idx = vconfs > thresholds(i);
%     TPR = sum( estim_idx &  truth_idx ) / P;
%     FPR = sum( estim_idx & ~truth_idx ) / N;
%     tprs(i) = TPR; fprs(i) = FPR;
% end
% dec_thresh = 0.45;
% tpr = sum( (vconfs>dec_thresh) &  truth_idx ) / P;
% fpr = sum( (vconfs>dec_thresh) & ~truth_idx ) / N;
% clear 300
% figure(300)
% hold on
% plot(fprs,tprs)
% plot(fpr,tpr,'ro')
% plot([0,1],[0,1],'r-')
% xlim([0,1])
% ylim([0,1])
% xlabel('False positive rate')
% ylabel('True positive rate')
% title('ROC Curve for VP Yaw Confidence')
% 
% figure(1000)
% plot(vconfs,vyerrs,'.')
% figure(1001)
% bins = 0:1:11;
% yebins = hist(vyerrs,bins);
% yebins = cumsum(yebins);
% bar(bins(1:end-1),yebins(1:end-1))
% xlim([-0.5,10.5])
% ylim([0,length(vyerrs)])
% 
% % Analyze vp normal estimates
% tnorms_copy = tnorms;
% rmidx = isnan(tnorms);
% tnorms(rmidx) = []; vnorms(rmidx,:) = []; vnconf(rmidx,:) = [];
% vn = zeros(0,1); vnc = zeros(0,1);
% i = 1; j = 1;
% while i <= length(tnorms)
%     numv = sum(~isnan(vnorms(j,:)));
%     vn = [vn,vnorms(j,1:numv)];
%     vnc = [vnc,vnconf(j,1:numv)];
%     tnorms = [tnorms(1:i-1);tnorms(i)*ones(numv,1);tnorms(i+1:end)];
%     i = i + numv; j = j + 1;
% end
% vnorms = vn'; vnconf = vnc';
% normerr = mod(vnorms-tnorms,90);
% normerr(normerr>45) = 90-normerr(normerr>45);
% clear 301
% figure(301)
% hold on
% plot(vnconf,normerr,'.')

% Analyze db normal estimates
if oakland
    rmidx = isnan(tnorms);
    tnorms(rmidx) = []; dnorms(rmidx,:) = []; dnerr(rmidx,:) = [];
    %dnerr
    dn_diff = mod( dnorms - repmat(tnorms,[1,7]) , 360 );
    dn_diff = dn_diff + (360-2*dn_diff) .* (dn_diff>180);
    dnerr = min(dn_diff,[],2);
    dnerr(isnan(dnerr)) = [];
    xmax  = 25;
    bins  = 0.5:1:xmax+0.5;
    dbins = hist(dnerr+0.1,bins);
    dsum  = cumsum(dbins);
    figure(14), clf, hold on
    bar(bins(1:end-1),100*dsum(1:end-1)/dsum(end),'r')
    bar(bins(bins<10),100*dsum(bins<10)/dsum(end),'g')
    xlabel('Bearing error (degrees)')
    ylabel('Percent of queries')
    title('Percent of Oakland yaw values within x degrees (database planes)')
    xlim([0,xmax])
    ylim([0,100])
end


% dn = zeros(0,1); dnc = zeros(0,1);
% i = 1; j = 1;
% while i <= length(tnorms)
%     numd = sum(~isnan(dnorms(j,:)));
%     dn = [dn,dnorms(j,1:numd)];
%     dnc = [dnc,dnerr(j,1:numd)];
%     tnorms = [tnorms(1:i-1);tnorms(i)*ones(numd,1);tnorms(i+1:end)];
%     i = i + numd; j = j + 1;
% end
% dnorms = dn'; dnerr = dnc';
% normerr = mod(dnorms-tnorms,90);
% normerr(normerr>45) = 90-normerr(normerr>45);
% clear 302
% figure(302)
% hold on
% plot(dnerr,normerr,'.')

cyerrs = mod(cyaws-tyaws,360);
cyerrs = cyerrs + (cyerrs>180) .* (360-2*cyerrs);
vyerrs = mod(vyaws-tyaws,360);
vyerrs = vyerrs + (vyerrs>180) .* (360-2*vyerrs);

xmax  = 25;
bins  = 0.5:1:xmax+0.5;
vbins = hist(vyerrs+0.1,bins);
cbins = hist(cyerrs+0.1,bins);
vsum  = cumsum(vbins);
csum  = cumsum(cbins);

figure(5+10*oakland), clf, hold on
bar(bins(1:end-1),100*vsum(1:end-1)/vsum(end),'r')
bar(bins(bins<10),100*vsum(bins<10)/vsum(end),'g')
xlabel('Bearing error (degrees)')
ylabel('Percent of queries')
if oakland, title('Percent of Oakland yaw values within x degrees (vanishing points)')
else        title('Percent of Berkeley yaw values within x degrees (vanishing points)')
end    
xlim([0,xmax])
ylim([0,100])

figure(6+10*oakland), clf, hold on
bar(bins(1:end-1),100*csum(1:end-1)/csum(end),'r')
bar(bins(bins<10),100*csum(bins<10)/csum(end),'g')
xlabel('Bearing error (degrees)')
ylabel('Percent of queries')
if oakland, title('Percent of Oakland yaw values within x degrees (compass)')
else        title('Percent of Berkeley yaw values within x degrees (compass)')
end    
xlim([0,xmax])
ylim([0,100])