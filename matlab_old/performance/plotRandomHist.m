function plotRandomHist(fignum,exclude,goodMask,nGood,nBad,vote_dir)

if nargin<5
    vote_dir = 'E:\q3results\';
end

[qNum,qName] = parseVotes(vote_dir);
nq = length(qNum);

gNum = [];
gName = [];
bNum = [];
bName = [];
for k=1:nq
    if strcmp(qMatch(k),'G') || strcmp(qMatch(k),'Y') || ...
            strcmp(qMatch(k),'R') || strcmp(qMatch(k),'B') || ...
            strcmp(qMatch(k),'O')
        gNum  = [gNum; qNum(k)];
        gName = [gName;qName(k,:)];
    else
        bNum  = [bNum; qNum(k)];
        bName = [bName;qName(k,:)];
    end
end
ng = length(gNum);
nb = length(bNum);

gPerm = randperm(ng);
bPerm = randperm(nb);

for k=1:nGood
    j=gPerm(k);
    vote = importdata([vote_dir,gName(j,:)]);
    vote = vote(:,1);
    size(vote)
    maxvote = max(vote);
    [n,x]=hist(vote,1:maxvote);
    figure(100+k)
    bar(x,n,'FaceColor','g','EdgeColor','k')
    title(['Photo #',num2str(gNum(j))])    
end

for k=1:nBad
    j=bPerm(k);
    vote = importdata([vote_dir,bName(j,:)]);
    vote = vote(:,1);
    maxvote = max(vote);
    [n,x]=hist(vote,1:maxvote);
    figure(200+k)
    bar(x,n,'FaceColor','r','EdgeColor','k')
    title(['Photo #',num2str(bNum(j))])    
end

end