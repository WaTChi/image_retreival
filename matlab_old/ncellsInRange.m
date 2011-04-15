function ncellsInRange(r,g,d,npts)

% assumes r = d
if nargin < 3
    d = r;
end
if nargin < 4
    npts = 101;
end

x = 0:d/(npts-1):d;
y = 0:d/(npts-1):d;

reach = d+r+g;

reachX = ceil( reach/d );
reachY = ceil( reach/(sqrt(3)*d/2) );

cellptX = [];
cellptY = [];
for pty = -reachY:reachY
    ptY = pty * sqrt(3)*d/2;
    if mod(pty,2) == 0
        for ptx = -reachX:reachX
            ptX = ptx * d;
            cellptX = [cellptX;ptX];
            cellptY = [cellptY;ptY];
        end
    else
        for ptx = -reachX+0.5:1:reachX-0.5
            ptX = ptx * d;
            cellptX = [cellptX;ptX];
            cellptY = [cellptY;ptY];
        end
    end
end

cell_mesh = zeros(npts,npts);
for j=1:npts
    for k=1:npts
        dists = sqrt( (cellptX-x(j)).^2 + (cellptY-y(k)).^2 );
        cell_mesh(j,k) = sum( (dists < r+g) );
    end
end

figure
imagesc(x,y,cell_mesh')
colorbar
axis equal
axis xy
x1 = 0; y1 = 0;
x2 = d; y2 = 0;
x3 = d/2; y3 = sqrt(3)*d/2;
hold on
plot([x1 x2 x3 x1],[y1,y2,y3 y1])
title(['g = ',num2str(g),' meters'])