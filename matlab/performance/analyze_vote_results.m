function analyze_vote_results(fignum,vote_dir,match_file)

if nargin<2
    vote_dir = 'E:\q3results\';
end
if nargin<3
    match_file = 'E:\Research\g=100,r=d=236.6cells,query3matches.txt';
end

[qNum,qName,qMatch] = vote_results(vote_dir,match_file);
nq = length(qNum);

qVote = cell(nq,1);
for k=1:nq
    vote = importdata([vote_dir,qName(k,:)]);
    qVote{k} = vote(:,1);
end

gNum = [];
gName = [];
gVote = cell(0,1);
bNum = [];
bName = [];
bVote = cell(0,1);
for k=1:nq
    if strcmp(qMatch(k),'G') || strcmp(qMatch(k),'Y')
        gNum  = [gNum; qNum(k)];
        gName = [gName;qName(k,:)];
        gVote = [gVote;qVote(k)];
    else
        bNum  = [bNum; qNum(k)];
        bName = [bName;qName(k,:)];
        bVote = [bVote;qVote(k)];
    end
end
ng = length(gNum);
nb = length(bNum);

% -------------------------------

cdfBins = 0 : .02 : .98 ;
histBins = .01 : .02 : .99 ;

qFeat = zeros(nq,1);
qMaxVote = zeros(nq,1);
qRatio = zeros(nq,1);
qCDF = zeros(size(cdfBins));
for k=1:nq
    qFeat(k) = sum(qVote{k});
    vote_norep = sort(unique(100*qVote{k}/qFeat(k)),1,'descend');
    qMaxVote(k) = vote_norep(1);
    qRatio(k) = vote_norep(2)/qMaxVote(k);
    histVote = (100*qVote{k}/qFeat(k))/qMaxVote(k);
    histVote = hist(histVote,histBins);
    qCDF = qCDF + fliplr(cumsum(fliplr(histVote)));
end
qCDF = qCDF / nq;

bFeat = zeros(nb,1);
gMaxVote = zeros(ng,1);
gRatio = zeros(ng,1);
gCDF = zeros(size(cdfBins));
for k=1:ng
    gFeat(k) = sum(gVote{k});
    vote_norep = sort(unique(100*gVote{k}/gFeat(k)),1,'descend');
    gMaxVote(k) = vote_norep(1);
    gRatio(k) = vote_norep(2)/gMaxVote(k);
    histVote = (100*gVote{k}/gFeat(k))/gMaxVote(k);
    histVote = hist(histVote,histBins);
    gCDF = gCDF + fliplr(cumsum(fliplr(histVote)));
end
gCDF = gCDF / ng;

bFeat = zeros(nb,1);
bMaxVote = zeros(nb,1);
bRatio = zeros(nb,1);
bCDF = zeros(size(cdfBins));
for k=1:nb
    bFeat(k) = sum(bVote{k});
    vote_norep = sort(unique(100*bVote{k}/bFeat(k)),1,'descend');
    bMaxVote(k) = vote_norep(1);
    bRatio(k) = vote_norep(2)/bMaxVote(k);
    histVote = (100*bVote{k}/bFeat(k))/bMaxVote(k);
    histVote = hist(histVote,histBins);
    bCDF = bCDF + fliplr(cumsum(fliplr(histVote)));
end
bCDF = bCDF / nb;

% Plot cdf of vote ratios
plotratio = 0;
for k=length(cdfBins):-1:1
    if cdfBins(k) >= plotratio
        plotidx = k;
    end
end
maxPhotos = max([qCDF(plotidx),gCDF(plotidx),bCDF(plotidx)]);
figure(3*fignum-2)
% All Queries
subplot(3,1,1)
bar(cdfBins(plotidx:end),qCDF(plotidx:end),'FaceColor','b')
axis([plotratio,1,0,maxPhotos])
title('CDF: all queries')
xlabel('Fraction of max vote number')
ylabel('Number of photos')
% Good Queries
subplot(3,1,2)
bar(cdfBins(plotidx:end),gCDF(plotidx:end),'FaceColor','g')
axis([plotratio,1,0,maxPhotos])
title('CDF: good queries')
xlabel('Fraction of max vote number')
ylabel('Number of photos')
% Bad Queries
subplot(3,1,3)
bar(cdfBins(plotidx:end),bCDF(plotidx:end),'FaceColor','r')
axis([plotratio,1,0,maxPhotos])
title('CDF: bad queries')
xlabel('Fraction of max vote number')
ylabel('Number of photos')

% Plot pdf of max vote number
nbins = 100;
minVote = min(qMaxVote);
maxVote = max(qMaxVote);
dbin = (maxVote-minVote)/nbins;
voteBins = minVote+dbin/2 : dbin : maxVote-dbin/2;
qmVote = hist(qMaxVote,voteBins)/nq;
gmVote = hist(gMaxVote,voteBins)/ng;
bmVote = hist(bMaxVote,voteBins)/nb;
maxY = max([max(qmVote),max(gmVote),max(bmVote)]);
figure(3*fignum-1)
% All Queries
subplot(3,1,1)
bar(voteBins,qmVote,'FaceColor','b')
axis([minVote,maxVote,0,maxY])
title('Max vote number: all queries')
xlabel('Max Vote Number')
% Good Queries
subplot(3,1,2)
bar(voteBins,gmVote,'FaceColor','g')
axis([minVote,maxVote,0,maxY])
title('Max vote number: good queries')
xlabel('Max Vote Number')
% Bad Queries
subplot(3,1,3)
bar(voteBins,bmVote,'FaceColor','r')
axis([minVote,maxVote,0,maxY])
title('Max vote number: bad queries')
xlabel('Max Vote Number')

% Plot pdf of ratio between 2nd highest vote number and highest
nbins = 50;
minRatio = min(qRatio);
dbin = (1-minRatio) / nbins;
ratioBins = minRatio + dbin/2 : dbin : 1-dbin/2;
qrVote = hist(qRatio,ratioBins)/nq;
grVote = hist(gRatio,ratioBins)/ng;
brVote = hist(bRatio,ratioBins)/nb;
maxY = max([max(qrVote),max(grVote),max(brVote)]);
figure(3*fignum)
% All Queries
subplot(3,1,1)
bar(ratioBins,qrVote,'FaceColor','b')
axis([minRatio,1,0,maxY])
title('Vote ratio: all queries')
xlabel('Second / First Highest Vote Number')
% Good Queries
subplot(3,1,2)
bar(ratioBins,grVote,'FaceColor','g')
axis([minRatio,1,0,maxY])
title('Vote ratio: good queries')
xlabel('Second / First Highest Vote Number')
% Bad Queries
subplot(3,1,3)
bar(ratioBins,brVote,'FaceColor','r')
axis([minRatio,1,0,maxY])
title('Vote ratio: bad queries')
xlabel('Second / First Highest Vote Number')

% figure(301)
% hist(qFeat)
figure(302)
hist(gFeat,50)
figure(303)
hist(bFeat,50)

end
