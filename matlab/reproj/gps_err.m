function loc_err(ncomp,fig_offset)

if ~exist('fig_offset','var')
    fig_offset = 0;
end

dx = 5; xmax = 100; ymax = 30;
r = importdata('Z:\ah\lastrun2\reproj.txt');
r = r(:,2:ncomp+1);
r(r>100) = 100;
disp(str(r))

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