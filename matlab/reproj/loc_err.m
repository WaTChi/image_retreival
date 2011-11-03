function PP_LocErr(ncomp,lrnum,fig_offset)

if ~exist('lrnum','var')
    lrnum = 1;
end

if ~exist('fig_offset','var')
    fig_offset = 10*(lrnum-1);
end

dx = 5; xmax = 100; ymax = 30;
if lrnum == 1
    r = importdata('Z:\ah\lastrun\reproj.txt');
else
    r = importdata(['Z:\ah\lastrun',num2str(lrnum),'\reproj.txt']);
end
r = r(:,2:ncomp+1);
r(r>100) = 100;
disp(num2str(r))

for k=1:ncomp
    figure(fig_offset+k)
    bins = dx/2:dx:xmax-dx/2;
    err = r(:,k);
    gidx = err < 10;
    gcount = hist(err(gidx),bins);
    bcount = hist(err(~gidx),bins);
    bar(bins,gcount,'g')
    hold on
    bar(bins,bcount,'r')
    xlim([0,xmax])
    ylim([0,ymax])
end