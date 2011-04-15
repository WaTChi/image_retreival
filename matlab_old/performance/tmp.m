idx=[];
for k=1:rIdx
    idx = [idx,find(gtIdx==k)];
end
gtIdx(idx) = [];

bi=[];
for k=yIdx
    bi = [bi,find(gtIdx==k,1)];
end

gi=[];
for k=gIdx
    gi = [gi,find(gtIdx==k,1)];
end

bi
gi

minG = min(scores(gi,:),[],2);
maxG = max(scores(gi,:),[],2);

minB = min(scores(bi,:),[],2);
maxB = max(scores(bi,:),[],2);

figure(1)
hold on
bar(gi,maxG,0.4,'FaceColor','g')
bar(gi,minG,0.4,'FaceColor','w')
bar(bi,maxB,0.4,'FaceColor','r')
bar(bi,minB,0.4,'FaceColor','w')