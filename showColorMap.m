function showColorMap(x,y,z,titleName,xLabel,yLabel,fontSize)
xSplit = diff(x)/2;                 % Find edge points
ySplit = diff(y)/2;
xEdges = [x(1)-xSplit(1) x(2:end)-xSplit x(end)+xSplit(end)];
yEdges = [y(1)-ySplit(1) y(2:end)-ySplit y(end)+ySplit(end)];
[XGrid, YGrid] = meshgrid(xEdges,yEdges);
%YGrid = flipud(YGrid);              % To match expected behavior
z = [[z zeros(size(z,1),1)] ; zeros(1,size(z,2)+1)];% Last row/col ignored

figure('DefaultAxesFontSize',fontSize);
theFigure = pcolor(XGrid,YGrid,z);
set(theFigure, 'EdgeColor', 'none');
% hold on                             % Plot original data points
% [X,Y] = meshgrid(x,y);
% %Y = flipud(Y);
% plot(X,Y,'or')
% xticks(unique(x));yticks(unique(y));
% caxis([0 1]);
colorbar;title(titleName);
xlabel(xLabel);ylabel(yLabel);
end