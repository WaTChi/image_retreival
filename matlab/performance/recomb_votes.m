gt_dir = 'E:\Research\app\code\matlab\';
cand_dir = 'E:\Research\query3top2\';
dbDir = 'E:\Research\collected_images\earthmine-new,culled\37.871955,-122.270829\';
qDir = 'E:\Research\collected_images\query\query3\';

run 'C:\vlfeat-0.9.9\toolbox\vl_setup'

[gtIdx,gtFile] = parseGT('GYRBO',gt_dir);
[candIdx,candFile] = parseCand(cand_dir);
nq = length(candIdx);
nc = length(candFile);

tic
% 10 is hard coded
scores = zeros(nq,2);
rIdx = [];
yIdx = [];
gIdx = [];
disp(' ')
for k=1:nq
    
    disp(['Recombination on query ',num2str(k)])
    
    querySift = getQuerySift(candIdx(k,1),qDir);
    [~,da] = vl_ubcread(strcat(qDir,querySift));
    
    topN = 1;
    cand = candFile( candIdx(k,2) : candIdx(k,3) );
    for j=1:length(cand)
        candSift = [cand{j},'sift.txt'];
        [~,db] = vl_ubcread(strcat(dbDir,candSift));
        scores(k,j) = length(vl_ubcmatch(da,db));
    end
    [~,perm] = sort(scores(k,:),2,'descend');
    winners = cand( perm(1:topN) );
    
    idx = find(candIdx(k,1)==gtIdx(:,1));
    gt = gtFile( gtIdx(idx,2) : gtIdx(idx,3) );
    if textMatch(winners,gt)
        gIdx = [gIdx,candIdx(k)];
    elseif textMatch(cand,gt)
        yIdx = [yIdx,candIdx(k)];
    else
        rIdx = [rIdx,candIdx(k)];
    end
    
end
disp(' ')
ng = length(gIdx);
ny = length(yIdx);
nr = length(rIdx);
disp(['Successful recombination for ',num2str(ng),' queries out of ',...
    num2str(nq),' (',num2str(dRound(100*ng/nq,0)),'%).'])
disp(['Filter lost ', num2str(ny),' queries out of ',...
    num2str(ng+ny),' (',num2str(dRound(100*ny/(ny+ng),0)),'%).'])
elapsed_time = toc;
disp(['Average match time of ',dRound(num2str(elapsed_time/nc),-1),' seconds.'])