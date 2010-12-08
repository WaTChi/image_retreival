function qRemove = analyzeVotes(fignum,exclude,goodMask,qFilter,vote_dir)

% fignum = 0 for no plots

if nargin<4
    vote_dir = 'E:\q3results\';
end

[qNum,qName] = parseQuery(vote_dir);
qNum(exclude) = [];
qName(exclude,:) = [];

% All queries
nq = length(qNum);
qVote = cell(nq,1);
qFeat = zeros(nq,1);
qMaxVote = zeros(nq,1);
qMaxPerc = zeros(nq,1);
for k=1:nq
    vote = importdata([vote_dir,qName(k,:)]);
    qVote{k} = vote(:,1);
    qFeat(k) = sum(qVote{k});
    qMaxVote(k) = max(qVote{k});
    qMaxPerc(k) = 100 * qMaxVote(k)/qFeat(k);
end

% Good queries
ng = sum(goodMask);
gNum = qNum(goodMask);
gName = qName(goodMask,:);
gVote = qVote(goodMask);
gFeat = zeros(ng,1);
gMaxVote = zeros(ng,1);
gMaxPerc = zeros(ng,1);
for k=1:ng
    gFeat(k) = sum(gVote{k});
    gMaxVote(k) = max(gVote{k});
    gMaxPerc(k) = 100 * gMaxVote(k)/gFeat(k);
end

% Bad queries
nb = sum(~goodMask);
bNum = qNum(~goodMask);
bName = qName(~goodMask,:);
bVote = qVote(~goodMask);
bFeat = zeros(nb,1);
bMaxVote = zeros(nb,1);
bMaxPerc = zeros(nb,1);
for k=1:nb
    bFeat(k) = sum(bVote{k});
    bMaxVote(k) = max(bVote{k});
    bMaxPerc(k) = 100 * bMaxVote(k)/bFeat(k);
end

% -------------------------------

% FILTER QUERIES
% all queries
qfMask = zeros(nq,1);
for k=1:nq
    if qMaxPerc(k) >= qFilter(1)
        qfMask(k) = 1;
    end
end
nqf = sum(qfMask);
pqf = dRound(100*nqf/nq,0);
disp([num2str(pqf),'% of all queries kept.'])
% good queries
gfMask = zeros(ng,1);
for k=1:ng
    if gMaxPerc(k) >= qFilter(1)
        gfMask(k) = 1;
    end
end
ngf = sum(gfMask);
pgf = dRound(100*ngf/ng,0);
disp([num2str(pgf),'% of good queries kept.'])
% bad queries
bfMask = zeros(nb,1);
for k=1:nb
    if bMaxPerc(k) >= qFilter(1)
        bfMask(k) = 1;
    end
end
nbf = sum(bfMask);
pbf = dRound(100*nbf/nb,0);
disp([num2str(100-pbf),'% of bad queries removed.'])

if logical(fignum)
    
    % Plot pdf of max vote
    nbins = 20;
    maxFeat = max(qFeat);
    dbin = maxFeat/nbins;
    featBins = dbin/2 : dbin : maxFeat-dbin/2;
    qnFeat = hist(qFeat,featBins)/nq;
    gnFeat = hist(gFeat,featBins)/ng;
    bnFeat = hist(bFeat,featBins)/nb;
    maxY = max([max(qnFeat),max(gnFeat),max(bnFeat)]);
    figure(fignum+1)
    % All Queries
    subplot(3,1,1)
    bar(featBins,qnFeat,'FaceColor','b')
    axis([0,maxFeat,0,maxY])
    title('Number of features: all queries')
    xlabel('Max Vote Number')
    % Good Queries
    subplot(3,1,2)
    bar(featBins,gnFeat,'FaceColor','g')
    axis([0,maxFeat,0,maxY])
    title('Number of features: all queries')
    xlabel('Max Vote Number')
    % Bad Queries
    subplot(3,1,3)
    bar(featBins,bnFeat,'FaceColor','r')
    axis([0,maxFeat,0,maxY])
    title('Number of features: all queries')
    xlabel('Max Vote Number')

    % Plot pdf of max vote
    nbins = 20;
    maxVote = 100;
	% maxVote = max(qMaxVote);
    dbin = maxVote/nbins;
    voteBins = dbin/2 : dbin : maxVote-dbin/2;
    qmVote = hist(qMaxVote,voteBins)/nq;
    gmVote = hist(gMaxVote,voteBins)/ng;
    bmVote = hist(bMaxVote,voteBins)/nb;
    maxY = max([max(qmVote),max(gmVote),max(bmVote)]);
    figure(fignum+2)
    % All Queries
    subplot(3,1,1)
    bar(voteBins,qmVote,'FaceColor','b')
    axis([0,maxVote,0,maxY])
    title('Max vote number: all queries')
    xlabel('Max vote number')
    % Good Queries
    subplot(3,1,2)
    bar(voteBins,gmVote,'FaceColor','g')
    axis([0,maxVote,0,maxY])
    title('Max vote number: good queries')
    xlabel('Max vote number')
    % Bad Queries
    subplot(3,1,3)
    bar(voteBins,bmVote,'FaceColor','r')
    axis([0,maxVote,0,maxY])
    title('Max vote number: bad queries')
    xlabel('Max vote number')

    % Plot pdf of max percent
    nbins = 20;
    maxPerc = 5;
    % maxPerc = max(qMaxPerc);
    dbin = maxPerc/nbins;
    percBins = dbin/2 : dbin : maxPerc-dbin/2;
    qmPerc = hist(qMaxPerc,percBins)/nq;
    gmPerc = hist(gMaxPerc,percBins)/ng;
    bmPerc = hist(bMaxPerc,percBins)/nb;
    maxY = max([max(qmPerc),max(gmPerc),max(bmPerc)]);
    figure(fignum+3)
    % All Queries
    subplot(3,1,1)
    hold off
    bar(percBins,qmPerc,'FaceColor','b')
    axis([0,maxPerc,0,maxY])
    hold on
    plot([qFilter(1),qFilter(1)],[0,maxY],'k')
    title('Max vote percent: all queries')
    xlabel('Max vote percent')
    % Good Queries
    subplot(3,1,2)
    hold off
    bar(percBins,gmPerc,'FaceColor','g')
    axis([0,maxPerc,0,maxY])
    hold on
    plot([qFilter(1),qFilter(1)],[0,maxY],'k')
    title('Max vote percent: good queries')
    xlabel('Max vote percent')
    % Bad Queries
    subplot(3,1,3)
    hold off
    bar(percBins,bmPerc,'FaceColor','r')
    axis([0,maxPerc,0,maxY])
    hold on
    plot([qFilter(1),qFilter(1)],[0,maxY],'k')
    title('Max vote percent: bad queries')
    xlabel('Max vote percent')
    
end

qRemove = ~logical(qfMask);

end
