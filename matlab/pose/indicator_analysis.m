function indicator_analysis

file = 'indicator_results.txt';


% Remove specified queries
fid = fopen(file);
names = textscan(fid,'%s\t%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
names = names{1};
idx = [find(~cellfun('isempty',strfind(names,'2011-10-28_13-00-41_047'))), ...
       find(~cellfun('isempty',strfind(names,'2011-04-04_15-17-45_763'))), ...
       find(~cellfun('isempty',strfind(names,'2011-04-04_15-22-31_544'))), ...
       find(~cellfun('isempty',strfind(names,'2011-10-28_12-23-24_368'))), ...
       find(~cellfun('isempty',strfind(names,'2011-10-28_12-13-41_589'))), ...
       find(~cellfun('isempty',strfind(names,'2011-10-28_12-12-29_132'))), ...
       find(~cellfun('isempty',strfind(names,'2011-10-28_12-12-20_282'))), ...
       find(~cellfun('isempty',strfind(names,'2011-10-28_12-11-42_811'))), ...
       find(~cellfun('isempty',strfind(names,'2011-10-28_12-09-45_828'))), ...
       find(~cellfun('isempty',strfind(names,'2011-10-28_12-09-35_578'))), ...
       find(~cellfun('isempty',strfind(names,'2011-10-28_11-59-53_593'))), ...
       find(~cellfun('isempty',strfind(names,'2011-10-28_12-01-35_821'))), ...
       find(~cellfun('isempty',strfind(names,'2011-04-04_14-56-08_926'))), ...
       find(~cellfun('isempty',strfind(names,'2011-04-04_15-06-58_888'))), ...
       find(~cellfun('isempty',strfind(names,'2011-04-04_15-19-11_784'))), ...
       find(~cellfun('isempty',strfind(names,'2011-04-04_15-32-55_864'))), ...
       find(~cellfun('isempty',strfind(names,'2011-04-04_15-34-01_368')))];
fclose(fid);
names(idx) = [];
2011-04-04_15-17-45_763

% Import pose results
results = importdata(file);
results(idx,:) = [];

rperr = results(:,2);
inl   = results(:,3);
yerr  = results(:,4);
yerr(yerr>50) = 50;
clrs = zeros(length(yerr),3);
clrs(yerr>2.5,1) = 1; clrs(yerr<25,2) = 0.5; clrs(yerr<10,2) = 1;
% clrs  = [yerr/50, 1-yerr/50, zeros([length(yerr),1])];

figure(200)
clear 200
hold on
xlim([0,2])
% ylim([0,2])
while ~isempty(yerr)
    idx = find(~cellfun('isempty',strfind(names,names{1})));
    rp = rperr(idx); in = inl(idx); ye = yerr(idx); cl = clrs(idx,:);
    [~,i] = min(rp); rp = rp / rp(i); [~,i] = max(in); in = in / in(i);
%     for j=1:length(rp)
%         plot(in(j),rp(j),'.','Color',cl(j,:))
%     end
    [~,i] = min(ye); plot(in(i),rp(i),'g.')
    names(idx) = []; rperr(idx) = []; inl(idx) = []; yerr(idx) = []; clrs(idx,:) = [];
end
% levels = .005:.005:0.1;
% inlrange = .01:.01:1;
% for level=levels
%     rprange = level*sqrt(inlrange);
%     plot(inlrange,rprange,'k-')
% end
hold off