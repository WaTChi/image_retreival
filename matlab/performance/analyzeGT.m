function [goodMask,noMatch,tmp] = ...
    analyzeGT(fignum,exclude,mFilter,vote_dir,gt_dir)

% fignum = 0 for no plots

if nargin<4
    vote_dir='E:\q3results\';
end
if nargin<5
    gt_dir='E:\Research\app\code\';
end

[qNum,qName] = parseQuery(vote_dir);
[gtIdx,gtFile] = parseGT('GY',gt_dir);
nq = length(gtIdx);


qNumVote = zeros(nq,1);
qRank = zeros(nq,1);
qRatio = zeros(nq,1);
qMaxVote = zeros(nq,1);
qMaxPerc = zeros(nq,1);
qMatch = zeros(nq,1);
goodMask = zeros(nq,1);
noMatch = [];
for k=1:nq
    [qVote,qPhoto] = parseVotes(qName(k,:),vote_dir);
    qNumVote(k) = sum(qVote);
    matches = gtFile(gtIdx(k,2):gtIdx(k,3));
    nm = length(matches);
    mVote = zeros(nm,1);
    for j=1:nm
        idx = find(strcmp(qPhoto,matches{j}));
        if length(idx)==1
            mVote(j) = qVote(idx);
        else
            mVote(j) = 0;
        end
    end
    if isempty(mVote)
        noMatch = [noMatch,k];
    else
        qMaxVote(k) = max(mVote);
        qRatio(k) = qMaxVote(k) / qVote(1);
        qMaxPerc(k) = 100 * qMaxVote(k) / qNumVote(k);
        if qMaxVote(k) ~= 0
            qRank(k) = find(qVote==qMaxVote(k),1,'last');
        else
            qRank(k) = 1000;
        end
        if (qRank(k)<=mFilter(1)) && (qRatio(k)>=mFilter(2)) && ...
                (qMaxVote(k)>=mFilter(3)) && (qMaxPerc(k)>=mFilter(4))
            goodMask(k) = 1;
        end
        % Find number of matches from filter
        ratioIdx = find( qVote/qVote(1) >= mFilter(2) , 1 , 'last' );
        voteIdx = find( qVote >= mFilter(3) , 1 , 'last' );
        percIdx = find( 100*qVote/qNumVote(k) >= mFilter(3) , 1 , 'last');
        qMatch(k) = min([ratioIdx,voteIdx,percIdx]);
        if qMatch(k) > mFilter(1)
            qMatch(k) = 0;
            tmp = qVote(1);
            for j=2:mFilter(1)+1
                if qVote(j)<tmp
                    qMatch(k) = j-1;
                end
            end
        end
    end
%     if k==48
%         qNum(k)
%         mVote
%         qMaxVote(k)
%         qRatio(k)
%         qMaxPerc(k)
%         qRank(k)
%         goodMask(k)
%         qRank(k)<=mFilter(1)
%     end
end
% Remove queries with no ground truth matches
qNum(noMatch) = [];
gtIdx(noMatch) = [];
qNumVote(noMatch) = [];
qRank(noMatch) = [];
qRatio(noMatch) = [];
qMaxVote(noMatch) = [];
qMaxPerc(noMatch) = [];
goodMask(noMatch) = [];
qMatch(noMatch) = [];
nq = length(qNum);
% Remove excluded queries
qNum(exclude) = [];
gtIdx(exclude) = [];
qNumVote(exclude) = [];
qRank(exclude) = [];
qRatio(exclude) = [];
qMaxVote(exclude) = [];
qMaxPerc(exclude) = [];
goodMask(exclude) = [];
qMatch(exclude) = [];
nq = length(qNum)

% Print ground truth performance
goodMask = logical(goodMask);
ngood = sum(goodMask);
pgood = dRound(100*ngood/nq,0);
disp([num2str(pgood),'% of queries find a ground truth match.'])
avgMatch = dRound(sum(qMatch)/nq,-1);
disp(['On average, ',num2str(avgMatch),' images are retained.'])

tmp = qNum(goodMask)
bad = qNum(~goodMask)

if logical(fignum)

    % Plot rank results : mFilter(1)
    figure(fignum+1)
    maxBin = 40;
    rankBins = 1:maxBin;
    histRank = hist(qRank,rankBins)/nq;
    cdfRank = cumsum(histRank);
    hold off
    bar(rankBins,cdfRank,'FaceColor','g');
    axis([0.5,maxBin+0.5,0,1])
    hold on
    plot([mFilter(1),mFilter(1)],[0,1],'r')
    xlabel('Maximum number of images retained')
    ylabel('Fraction of queries with a good match')

    % Plot ratio results : mFilter(2)
    figure(fignum+2)
    nbins = 40;
    dbin = 1 / nbins;
    ratioBins = dbin/2 : dbin : 1-dbin/2;
    histRatio = hist(qRatio,ratioBins)/nq;
    cdfRatio = fliplr(cumsum(fliplr(histRatio)));
    hold off
    bar(ratioBins,cdfRatio,'FaceColor','g');
    axis([0,1,0,1])
    hold on
    plot([mFilter(2),mFilter(2)],[0,1],'r')
    xlabel('Ratio threshold')
    ylabel('Fraction of queries with a good match')

    % Plot vote results : mFilter(3)
    figure(fignum+3)
    nbins = 40;
    mxVote = max(qMaxVote);
    dbin = mxVote / nbins;
    mvBins = dbin/2 : dbin : mxVote - dbin/2;
    histMV = hist(qMaxVote,mvBins)/nq;
    cdfMV = fliplr(cumsum(fliplr(histMV)));
    hold off
    bar(mvBins,cdfMV,'FaceColor','g');
    axis([0,mxVote,0,1])
    hold on
    plot([mFilter(3),mFilter(3)],[0,1],'r')
    xlabel('Threshold for number of total votes')
    ylabel('Fraction of queries with a good match')

    % Plot percent results : mFilter(4)
    figure(fignum+4)
    nbins = 40;
    mxPerc = max(qMaxPerc);
    dbin = mxPerc / nbins;
    mpBins = dbin/2 : dbin : mxPerc - dbin/2;
    histMP = hist(qMaxPerc,mpBins)/nq;
    cdfMP = fliplr(cumsum(fliplr(histMP)));
    hold off
    bar(mpBins,cdfMP,'FaceColor','g');
    axis([0,mxPerc,0,1])
    hold on
    plot([mFilter(4),mFilter(4)],[0,1],'r')
    xlabel('Threshold for percentage of total votes')
    ylabel('Fraction of queries with a good match')
    
end

end