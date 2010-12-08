load dec2_perf

gtIdx2 = gtIdx(:,1);
idx=[];
for k=rIdx
    idx = [idx,find(gtIdx2==k,1)];
end
gtIdx2(idx) = [];


bi=[];
for k=yIdx
    bi = [bi,find(gtIdx2==k,1)];
end

gi=[];
for k=gIdx
    gi = [gi,find(gtIdx2==k,1)];
end

minG = min(scores(gi,:),[],2);
maxG = max(scores(gi,:),[],2);

minB = min(scores(bi,:),[],2);
maxB = max(scores(bi,:),[],2);

figure(1)
hold on
bar(gi,maxG,0.8,'FaceColor','g')
bar(gi,minG,0.8,'FaceColor','w','EdgeColor','w')
bar(bi,maxB,0.8,'FaceColor','r')
bar(bi,minB,0.8,'FaceColor','w','EdgeColor','w')