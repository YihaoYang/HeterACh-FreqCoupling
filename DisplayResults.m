% Scrips to generate figures in the paper

%% Figure-2 - pre
clear all;
load('random_gks_data.mat')
close all;i = 3;j = 3;k = 1;
myFont = 16;
%% Figure-2a
fig = figure('DefaultAxesFontSize',myFont);
gKsHeatmapMatrix = ...
    reshape(net(i,j,k).gKs(1:4*sideNumI^2),[2*sideNumI,2*sideNumI]);
pcolor(gKsHeatmapMatrix);caxis([0 1.5]);
pbaspect([1 1 1]);
set(gca,'xtick',[1 20])
set(gca,'ytick',[1 20])
xlb = xlabel('x');
xlb.Position(2) = xlb.Position(2)+1;
ylb = ylabel('y');
ylb.Position(1) = ylb.Position(1)+1;
% pbaspect([1 1 1]);
colormap parula;
tempCB = colorbar;
caxis([0 1.5])
title(tempCB,'g_{Ks}')

%% Figure-2b
fig = figure('DefaultAxesFontSize',myFont);
net(i,j,k).showColorCodedRasterPlot;
set(gca,'xtick',[4000 5000])
set(gca,'ytick',[1 400 500])
%             hold on; plot(0:5000,400*ones(1,5001));hold off;
xlim([4000 5000]);
% xlabel('time/ms')
xlb = xlabel("t/ms");
xlb.Position(2) = -30;
ylb = ylabel('neuron ID');
ylb.Position(1) = ylb.Position(1)+30;
plot(0:5000,400*ones(1,5001));
box on
grid on
fig.OuterPosition(3) = 2*fig.OuterPosition(3);
axis fill

%% Figure-2c spiking num map
fig = figure('DefaultAxesFontSize',myFont);
spikeTimes = net(i,j,k).spikeTimes';
spikeTimes = spikeTimes(spikeTimes(:,2) <= net(i,j,k).numE ,:);
spikeNum = zeros(net(i,j,k).numE,1);
for neuronId = 1:net(i,j,k).numE
spikeNum(neuronId) = ...
       sum(spikeTimes(:,2) == neuronId);
end
heatmapMatrix = reshape(spikeNum,[2*sideNumI,2*sideNumI]);
pcolor(heatmapMatrix/4);

tempCB = colorbar;
title(tempCB,'$\overline{f}(Hz)$','interpreter','latex')
set(gca,'xtick',[1 20])
set(gca,'ytick',[1 20])
xlb = xlabel('x');
xlb.Position(2) = xlb.Position(2)+1;
ylb = ylabel('y');
ylb.Position(1) = ylb.Position(1)+1;
pbaspect([1 1 1]);
colormap parula;

%% Figure-2d type map
fig = figure('DefaultAxesFontSize',myFont);
net(i,j,k).mapRhythmVec1();
RhythmHeatmapMatrix = ...
    reshape(net(i,j,k).neuronRhythmVec(1:4*sideNumI^2),[2*sideNumI,2*sideNumI]);
pcolor(RhythmHeatmapMatrix);caxis([0,3]);
pbaspect([1 1 1]);
colormap parula;
% colorbar;
set(gca,'xtick',[1 20])
set(gca,'ytick',[1 20])
xlb = xlabel('x');
xlb.Position(2) = xlb.Position(2)+1;
ylb = ylabel('y');
ylb.Position(1) = ylb.Position(1)+1;

%% Figure-2e
figure('DefaultAxesFontSize',myFont);
net(i,j,k).detectDynamics;
grid on
box on
disp(['theta = ',num2str(net(i,j,k).thetaPeak),newline, ...
    'gamma = ', num2str(net(i,j,k).gammaPeak)]);
fig = gcf;
fig.OuterPosition(4) = fig.OuterPosition(4)/2;
axis fill
%% Figure-2f
fig = figure('DefaultAxesFontSize',myFont);
noneGks = gKsHeatmapMatrix(RhythmHeatmapMatrix == 0);
noneGksCount = hist(noneGks,0:0.1:1.5)';
thetaGks = gKsHeatmapMatrix(RhythmHeatmapMatrix == 1);
thetaGksCount = hist(thetaGks,0:0.1:1.5)';
gammaGks = gKsHeatmapMatrix(RhythmHeatmapMatrix == 2);
gammaGksCount = hist(gammaGks,0:0.1:1.5)';
mixedGks = gKsHeatmapMatrix(RhythmHeatmapMatrix == 3);
mixedGksCount = hist(mixedGks,0:0.1:1.5)';
barArray = bar(0:0.1:1.5,[noneGksCount,thetaGksCount,gammaGksCount,mixedGksCount],'stacked');
h0 = barArray(1);h1 = barArray(2);h2 = barArray(3);h3 = barArray(4);
% legend('none','theta','gamma','mixed','Location','best')
xlabel('g_{Ks}');ylabel('counts');
cmap = colormap;
[tempVal,~] = size(cmap);
h0.FaceColor = cmap(1,:);
h1.FaceColor = cmap(round(tempVal/3),:);
h2.FaceColor = cmap(round(2*tempVal/3),:);
h3.FaceColor = cmap(tempVal,:);
fig.OuterPosition(4) = fig.OuterPosition(4)/2;
axis fill

%% Figure-4
clear all
load('double_peaked_gks_data.mat')
close all;
myFont = 16;

%% Figure-4a
i = 1;j = 4;k = 3;
fig = figure('DefaultAxesFontSize',myFont);
gKsHeatmapMatrix = ...
    reshape(net(i,j,k).gKs(1:4*sideNumI^2),[2*sideNumI,2*sideNumI]);
pcolor(gKsHeatmapMatrix);caxis([0 1.5]);
set(gca,'xtick',[1 20])
set(gca,'ytick',[1 20])
xlb = xlabel('x');
xlb.Position(2) = xlb.Position(2)+1;
ylb = ylabel('y');
ylb.Position(1) = ylb.Position(1)+1;
pbaspect([1 1 1]);
colormap parula;
tempCB = colorbar;
title(tempCB,'g_{Ks}')


%% Figure-4b
i = 1;j = 4;k = 3;
figure('DefaultAxesFontSize',myFont);
net(i,j,k).mapRhythmVec();
RhythmHeatmapMatrix = ...
    reshape(net(i,j,k).neuronRhythmVec(1:4*sideNumI^2),[2*sideNumI,2*sideNumI]);
pcolor(RhythmHeatmapMatrix);caxis([0,3]);
pbaspect([1 1 1]);
set(gca,'xtick',[1 20])
set(gca,'ytick',[1 20])
xlb = xlabel('x');
xlb.Position(2) = xlb.Position(2)+1;
ylb = ylabel('y');
ylb.Position(1) = ylb.Position(1)+1;
colormap parula;


%% Figure-4cde pre
thetaHeatmap = zeros(numX,numY,numRepetition);
gammaHeatmap = zeros(numX,numY,numRepetition);
mixedHeatmap = zeros(numX,numY,numRepetition);
thetaFreqMat = zeros(numX,numY,numRepetition);
gammaFreqMat = zeros(numX,numY,numRepetition);
for k = 1:numRepetition
    for i = 1:numX
        for j = 1:numY
%             net(i,j,k).run2DNet();
            net(i,j,k).mapRhythmVec();
            RhythmHeatmapMatrix = ...
                reshape(net(i,j,k).neuronRhythmVec(1:4*sideNumI^2),[2*sideNumI,2*sideNumI]);            
            thetaHeatmap(i,j,k) = sum((RhythmHeatmapMatrix == 1),'all');
            gammaHeatmap(i,j,k) = sum((RhythmHeatmapMatrix == 2),'all');
            mixedHeatmap(i,j,k) = sum((RhythmHeatmapMatrix == 3),'all');
            net(i,j,k).detectDynamics();
            thetaFreqMat(i,j,k) = net(i,j,k).thetaFreq;
            gammaFreqMat(i,j,k) = net(i,j,k).gammaFreq;
            thetaMat(i,j,k) = net(i,j,k).thetaPeak;
            gammaMat(i,j,k) = net(i,j,k).gammaPeak;
            relativeStrengthMat(i,j,k) = net(i,j,k).relativeStrength;
        end
    end
end
% double count
thetaHeatmap = thetaHeatmap + mixedHeatmap;
gammaHeatmap = gammaHeatmap + mixedHeatmap;

%% Supplementary: theta frequency map
thetaFreqMat(1:4,1,:) = nan;
thetaFreqMat(1,2,:) = nan;

showColorMap(radiusPoints, ...
    distancePoints, ...
    mean(thetaFreqMat,3)', ...
    'theta freq','r','d',myFont)
caxis([3 8])

%% Supplementary: gamma frequency map
gammaFreqMat(1,1,:) = nan;

showColorMap(radiusPoints, ...
    distancePoints, ...
    mean(gammaFreqMat,3)', ...
    'gamma freq','r','d',myFont)


%% Figure-4c
showColorMap(radiusPoints, ...
    distancePoints, ...
    mean(thetaHeatmap,3)', ...
    '# theta neurons','r','d',myFont)


%% Figure-4d
showColorMap(radiusPoints, ...
    distancePoints, ...
    mean(gammaHeatmap,3)', ...
    '# gamma neurons','r','d',myFont)

caxis([27.5 52.5])

%% Figure-4e
showColorMap(radiusPoints, ...
    distancePoints, ...
    mean(mixedHeatmap,3)', ...
    '# mixed neurons','r','d',myFont)

caxis([0 51])

%% Figure-4f1
i = 1;j = 2;k = 2;
fig = figure('DefaultAxesFontSize',myFont);
% title('(f) raster plot')
net(i,j,k).showColorCodedRasterPlot;
% set(gca,'xtick',[])
% set(gca,'ytick',[])
plot(0:5000,400*ones(1,5001));
xlim([4500 5000]);
box on
fig.OuterPosition(4) = fig.OuterPosition(4)/2;

%             title('(b)', 'Position', [-1, 100, 1]);
% tempCB = colorbar;
caxis([0 1.5])
% title(tempCB,'gKs')
axis fill

%% Figure-4f2
i = 1;j = 2;k = 2;
RhythmHeatmapMatrix = ...
                reshape(net(i,j,k).neuronRhythmVec(1:4*sideNumI^2),[2*sideNumI,2*sideNumI]);            
gKsHeatmapMatrix = ...
    reshape(net(i,j,k).gKs(1:4*sideNumI^2),[2*sideNumI,2*sideNumI]);
fig = figure('DefaultAxesFontSize',myFont);
% noneGks = gKsHeatmapMatrix(RhythmHeatmapMatrix == 0);
% noneGksCount = hist(noneGks,0:0.1:1.5)';
thetaGks = gKsHeatmapMatrix(RhythmHeatmapMatrix == 1);
thetaGksCount = hist(thetaGks,0:0.1:1.5)';
gammaGks = gKsHeatmapMatrix(RhythmHeatmapMatrix == 2);
gammaGksCount = hist(gammaGks,0:0.1:1.5)';
mixedGks = gKsHeatmapMatrix(RhythmHeatmapMatrix == 3);
mixedGksCount = hist(mixedGks,0:0.1:1.5)';
barArray = bar(0:0.1:1.5,[thetaGksCount,gammaGksCount,mixedGksCount],'stacked');
% h0 = barArray(1);
h1 = barArray(1);h2 = barArray(2);h3 = barArray(3);
legend('theta','gamma','mixed','Location','best')
xlabel('gKs');ylabel('counts');
cmap = colormap;
[tempVal,~] = size(cmap);
% h0.FaceColor = cmap(1,:);
h1.FaceColor = cmap(round(tempVal/3),:);
h2.FaceColor = cmap(round(2*tempVal/3),:);
h3.FaceColor = cmap(tempVal,:);
fig.OuterPosition(4) = fig.OuterPosition(4)/2;
axis fill

%% Figure-4g1
i = 2;j = 2;k = 2;
fig = figure('DefaultAxesFontSize',myFont);
% title('(g) raster plot')
net(i,j,k).showColorCodedRasterPlot;

plot(0:5000,400*ones(1,5001));
xlim([4500 5000]);
box on
fig.OuterPosition(4) = fig.OuterPosition(4)/2;
axis fill
%             title('(b)', 'Position', [-1, 100, 1]);
% tempCB = colorbar;
caxis([0 1.5])
% title(tempCB,'gKs')


%% Figure-4g2
i = 2;j = 2;k = 2;
RhythmHeatmapMatrix = ...
                reshape(net(i,j,k).neuronRhythmVec(1:4*sideNumI^2),[2*sideNumI,2*sideNumI]);            
gKsHeatmapMatrix = ...
    reshape(net(i,j,k).gKs(1:4*sideNumI^2),[2*sideNumI,2*sideNumI]);
fig = figure('DefaultAxesFontSize',myFont);
% noneGks = gKsHeatmapMatrix(RhythmHeatmapMatrix == 0);
% noneGksCount = hist(noneGks,0:0.1:1.5)';
thetaGks = gKsHeatmapMatrix(RhythmHeatmapMatrix == 1);
thetaGksCount = hist(thetaGks,0:0.1:1.5)';
gammaGks = gKsHeatmapMatrix(RhythmHeatmapMatrix == 2);
gammaGksCount = hist(gammaGks,0:0.1:1.5)';
mixedGks = gKsHeatmapMatrix(RhythmHeatmapMatrix == 3);
mixedGksCount = hist(mixedGks,0:0.1:1.5)';
barArray = bar(0:0.1:1.5,[thetaGksCount,gammaGksCount,mixedGksCount],'stacked');
% h0 = barArray(1);
h1 = barArray(1);h2 = barArray(2);h3 = barArray(3);
legend('theta','gamma','mixed','Location','best')
xlabel('gKs');ylabel('counts');
cmap = colormap;
[tempVal,~] = size(cmap);
% h0.FaceColor = cmap(1,:);
h1.FaceColor = cmap(round(tempVal/3),:);
h2.FaceColor = cmap(round(2*tempVal/3),:);
h3.FaceColor = cmap(tempVal,:);
fig.OuterPosition(4) = fig.OuterPosition(4)/2;
axis fill



%% Figure-4h1
i = 8;j = 4;k = 1;
fig = figure('DefaultAxesFontSize',myFont);
% title('(h) raster plot')
net(i,j,k).showColorCodedRasterPlot;
plot(0:5000,400*ones(1,5001));
xlim([4500 5000]);
box on
fig.OuterPosition(4) = fig.OuterPosition(4)/2;
axis fill
%             title('(b)', 'Position', [-1, 100, 1]);
% tempCB = colorbar;
caxis([0 1.5])
% title(tempCB,'gKs')


%% Figure-4h2
i = 8;j = 4;k = 1;
RhythmHeatmapMatrix = ...
                reshape(net(i,j,k).neuronRhythmVec(1:4*sideNumI^2),[2*sideNumI,2*sideNumI]);            
gKsHeatmapMatrix = ...
    reshape(net(i,j,k).gKs(1:4*sideNumI^2),[2*sideNumI,2*sideNumI]);
fig = figure('DefaultAxesFontSize',myFont);
% noneGks = gKsHeatmapMatrix(RhythmHeatmapMatrix == 0);
% noneGksCount = hist(noneGks,0:0.1:1.5)';
thetaGks = gKsHeatmapMatrix(RhythmHeatmapMatrix == 1);
thetaGksCount = hist(thetaGks,0:0.1:1.5)';
gammaGks = gKsHeatmapMatrix(RhythmHeatmapMatrix == 2);
gammaGksCount = hist(gammaGks,0:0.1:1.5)';
mixedGks = gKsHeatmapMatrix(RhythmHeatmapMatrix == 3);
mixedGksCount = hist(mixedGks,0:0.1:1.5)';
barArray = bar(0:0.1:1.5,[thetaGksCount,gammaGksCount,mixedGksCount],'stacked');
% h0 = barArray(1);
h1 = barArray(1);h2 = barArray(2);h3 = barArray(3);
legend('theta','gamma','mixed','Location','best')
xlabel('gKs');ylabel('counts');
cmap = colormap;
[tempVal,~] = size(cmap);
% h0.FaceColor = cmap(1,:);
h1.FaceColor = cmap(round(tempVal/3),:);
h2.FaceColor = cmap(round(2*tempVal/3),:);
h3.FaceColor = cmap(tempVal,:);
fig.OuterPosition(4) = fig.OuterPosition(4)/2;
axis fill


%% Figure-3a
i = 5;j = 1;k = 1;
fig = figure('DefaultAxesFontSize',myFont);
gKsHeatmapMatrix = ...
    reshape(net(i,j,k).gKs(1:4*sideNumI^2),[2*sideNumI,2*sideNumI]);
pcolor(gKsHeatmapMatrix);caxis([0 1.5]);
set(gca,'xtick',[1 20])
set(gca,'ytick',[1 20])
xlb = xlabel('x');
xlb.Position(2) = xlb.Position(2)+1;
ylb = ylabel('y');
ylb.Position(1) = ylb.Position(1)+1;
pbaspect([1 1 1]);
colormap parula;
tempCB = colorbar;
caxis([0 1.5])
title(tempCB,'g_{Ks}')
% title('(a) gks map');

%% Figure-3b
% reverse double count
thetaHeatmap = thetaHeatmap - mixedHeatmap;
gammaHeatmap = gammaHeatmap - mixedHeatmap;
figure('DefaultAxesFontSize',myFont);
hold all
box on
xlabel('r')

e2 = errorbar(radiusPoints,...
mean(gammaHeatmap(:,1,:),3),...
std(gammaHeatmap(:,1,:),0,3),...
'linewidth',2,'DisplayName','# gamma only cells');

e2.Color = '#77AC30';
hold on
e3 = errorbar(radiusPoints,...
    mean(thetaHeatmap(:,1,:),3),...
    std(thetaHeatmap(:,1,:),0,3),...
    'linewidth',2,'DisplayName','# theta only cells');

e3.Color = '#0072BD';

e1 = errorbar(radiusPoints,...
mean(mixedHeatmap(:,1,:),3),...
std(mixedHeatmap(:,1,:),0,3),...
'linewidth',2,'DisplayName','# mixed cells');

legend('Location','best');



%% Figure-3c : type map
i = 5;j = 1;k = 1;
fig = figure('DefaultAxesFontSize',myFont);
net(i,j,k).mapRhythmVec();
RhythmHeatmapMatrix = ...
    reshape(net(i,j,k).neuronRhythmVec(1:4*sideNumI^2),[2*sideNumI,2*sideNumI]);
pcolor(RhythmHeatmapMatrix);caxis([0,3]);
pbaspect([1 1 1]);
colormap parula;
% colorbar;
% title('(e) neuron type');
set(gca,'xtick',[1 20])
set(gca,'ytick',[1 20])
xlb = xlabel('x');
xlb.Position(2) = xlb.Position(2)+1;
ylb = ylabel('y');
ylb.Position(1) = ylb.Position(1)+1;



%% Figure-3d new- type map
i = 7;j = 1;k = 1;
fig = figure('DefaultAxesFontSize',myFont);
net(i,j,k).mapRhythmVec();
RhythmHeatmapMatrix = ...
    reshape(net(i,j,k).neuronRhythmVec(1:4*sideNumI^2),[2*sideNumI,2*sideNumI]);
pcolor(RhythmHeatmapMatrix);caxis([0,3]);
pbaspect([1 1 1]);
colormap parula;
% colorbar;
% title('(e) neuron type');
set(gca,'xtick',[1 20])
set(gca,'ytick',[1 20])
xlb = xlabel('x');
xlb.Position(2) = xlb.Position(2)+1;
ylb = ylabel('y');
ylb.Position(1) = ylb.Position(1)+1;

%% Figure-5
clear all;
load('random_gks_data.mat')
close all;
k = 1;i = 8;j = 2;  % star one
i=1;j=2;            % square one
myFont = 16;
%% Figure-5a-1-4
for k = 1:4
fig = figure('DefaultAxesFontSize',myFont);
gKsHeatmapMatrix = ...
    reshape(net(i,j,k).gKs(1:4*sideNumI^2),[2*sideNumI,2*sideNumI]);
pcolor(gKsHeatmapMatrix);caxis([0 1.5]);
pbaspect([1 1 1]);
set(gca,'xtick',[])
set(gca,'ytick',[])
% pbaspect([1 1 1]);
colormap parula;
tempCB = colorbar;
caxis([0 1.5])
title(tempCB,'g_{Ks}')
% title('(a) gks map');

end
%% Figure-5b-1-4
for k=1:4
fig = figure('DefaultAxesFontSize',myFont);
net(i,j,k).showColorCodedRasterPlot;
set(gca,'xtick',[3000 3600])
set(gca,'ytick',[])
%             hold on; plot(0:5000,400*ones(1,5001));hold off;
xlim([3000 3600]);
% xlb = xlabel("[spike time - t/ms]");
% xlb.Position(2) = -30;
ylabel('neuron ID');
plot(0:5000,400*ones(1,5001));
box on
grid on
fig.OuterPosition(3) = 2*fig.OuterPosition(3);
axis fill

end
%% Figure-5c-1-4
for k=1:4
figure('DefaultAxesFontSize',myFont);
net(i,j,k).detectDynamics;
grid on
box on
%             set(gca,'YAxisLocation','right')
%             set(gca,'XAxisLocation','origin')
% title('(c) rhythm spectrum')

disp(['theta = ',num2str(net(i,j,k).thetaPeak),newline, ...
    'gamma = ', num2str(net(i,j,k).gammaPeak)]);
fig = gcf;
fig.OuterPosition(4) = fig.OuterPosition(4)/2;
axis fill
end

%% Figure-6B
load('double_peaked_gks_data.mat')
load('adjMatForNearest12.mat')
i = 8;j=4;
rng("default")
distancePoints = 0:6;
correlationPoints = zeros(numRepetition,length(distancePoints));
MIpoints = zeros(numRepetition,length(distancePoints));
close all;
for k = 1:numRepetition
        net(i,j,k).run2DNet();
        net(i,j,k).detectDynamics();
        figure();
        net(i,j,k).showColorCodedRasterPlot();
        [tPoints,traces] = net(i,j,k).computeSpikeTraces();
    for whichDistance = 1:length(distancePoints)
        distanceToHotSpot = distancePoints(whichDistance);
        id = sub2ind([20 20],7+distanceToHotSpot,7);
        LFP = adjMatForNearest12(id,1:400) * traces;        
        % power spectrum part
        fs = 1000/(tPoints(2)-tPoints(1)); %sampling freq
        freq = 1:100; %power spectrum freq      
        
        % bandpass %%%%%%%%%%%%%%%%% bandpass %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        thetaFilteredLFP = bandpass(LFP,[net(i,j,k).thetaFreq-0.5 net(i,j,k).thetaFreq+0.5],fs);
        gammaFilteredLFP = bandpass(LFP,[net(i,j,k).gammaFreq-2 net(i,j,k).gammaFreq+2],fs);
        fig = figure('DefaultAxesFontSize',26,'name','bandpass');
        subplot(3, 1, 1)        
        plot(tPoints,LFP,'lineWidth',2);
        title(['distance to hotspot = ',num2str(distanceToHotSpot)])
        subplot(3, 1, 2)
        hold all
        plot(tPoints,thetaFilteredLFP,'lineWidth',2)
        plot(tPoints,gammaFilteredLFP,'lineWidth',1)
        gammaEnv = abs(hilbert(gammaFilteredLFP));
        plot(tPoints,gammaEnv,'lineWidth',2)  
        subplot(3,1,3)
        % MI-modulation Index %%%%%%%%%%%%%%%%%%%%%%%
        thetaPhase = angle(hilbert(thetaFilteredLFP));
        plot(tPoints,thetaPhase)
        numPhaseBins = 18;
        meanAmp_gamma = zeros(1,numPhaseBins);
        for ii=1:numPhaseBins
            meanAmp_gamma(ii) = mean(gammaEnv(thetaPhase>(ii-1)*pi/numPhaseBins & thetaPhase <= ii*pi/numPhaseBins));
        end
        ampDistribution = meanAmp_gamma/sum(meanAmp_gamma);
        D_KL_p2u = log(numPhaseBins) + ampDistribution * log(ampDistribution)';
        MI = D_KL_p2u/log(numPhaseBins);
        MIpoints(k,whichDistance) = MI;
        correlationPoints(k,whichDistance) = xcorr(thetaFilteredLFP,gammaEnv,0,'coeff'); % ,'coeff'
    end    
end
fig = figure('name','correlation','DefaultAxesFontSize',26);
box on
grid on
errorbar(distancePoints,mean(correlationPoints,1),std(correlationPoints,0,1)/sqrt(numRepetition),'linewidth',2);
xlabel('distance[pixel]')
ylabel('correlation')
xlim([distancePoints(1) distancePoints(end)])
xticks(distancePoints)
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_LFP_corr','-dpdf','-r800')

fig = figure('name','MI','DefaultAxesFontSize',26);
errorbar(distancePoints,mean(MIpoints,1),std(MIpoints,0,1)/sqrt(numRepetition),'linewidth',2);
xlabel('distance[pixel]')
ylabel('MI')
xlim([distancePoints(1) distancePoints(end)])
xticks(distancePoints)
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_LFP_MI','-dpdf','-r800')
%% Figure-6 D-E
k = 1;
[tPoints,traces] = net(i,j,k).computeSpikeTraces();
for whichDistance = 1:length(distancePoints)
    distanceToHotSpot = distancePoints(whichDistance);
    id = sub2ind([20 20],7+distanceToHotSpot,7);
    LFP = adjMatForNearest12(id,1:400) * traces;
    % power spectrum part
    fs = 1000/(tPoints(2)-tPoints(1)); %sampling freq
    freq = 1:100; %power spectrum freq    
    % bandpass %%%%%%%%%%%%%%%%% bandpass %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    thetaFilteredLFP = bandpass(LFP,[net(i,j,k).thetaFreq-0.5 net(i,j,k).thetaFreq+0.5],fs);
    gammaFilteredLFP = bandpass(LFP,[net(i,j,k).gammaFreq-2 net(i,j,k).gammaFreq+2],fs);
    fig = figure('DefaultAxesFontSize',26,'name','bandpass');
    box on
    grid on
    hold all
    plot(tPoints,thetaFilteredLFP,'lineWidth',2)
    plot(tPoints,gammaFilteredLFP,'lineWidth',1)
    gammaEnv = abs(hilbert(gammaFilteredLFP));
    plot(tPoints,gammaEnv,'lineWidth',2)
    legend('Filtered Theta','Filtered Gamma','Gamma envolope')
    fig.OuterPosition(3) = 2*fig.OuterPosition(3);
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig,['fig_LFP_filtered',num2str(whichDistance)],'-dpdf','-r800')
end

%% Figure-7 - Pre
clear all
rng("default")
load('double_peaked_gks_data.mat')
load('adjMatForNearest12.mat')
i = 8;

dPtMat = zeros(4,3,2);

dPgMat = zeros(4,3,2);

dFtMat = zeros(4,3,2);

dFgMat = zeros(4,3,2);

dFmeanMat = zeros(4,3,2);

dFmean4allMat = zeros(4,3,2);
for j = [2 1 4]
    if j == 2
        computing = 0;
    elseif j == 1
        computing = 1;
    else
        computing = 2;
    end
for k = 1:2
    
    if computing == 0
        net(i,j,k).gKs(:) = 1.5;
    end
    
    net(i,j,k).run2DNet();
    net(i,j,k).detectDynamics();
    Pt = net(i,j,k).thetaPeak;
    Pg = net(i,j,k).gammaPeak;
    Ft = net(i,j,k).thetaFreq;
    Fg = net(i,j,k).gammaFreq;
    
    
    spikeTimes = net(i,j,k).spikeTimes';
    spikeTimes = spikeTimes(spikeTimes(:,2) <= net(i,j,k).numE ,:);
    rawSpikeNum = zeros(net(i,j,k).numE,1);
    for neuronId = 1:net(i,j,k).numE
    rawSpikeNum(neuronId) = ...
           sum(spikeTimes(:,2) == neuronId);
    end
    
    
    distanceToHotSpot = 0;%10; %6 for outside for one hotspot
    if computing == 2        
        id = sub2ind([20 20],7+distanceToHotSpot,7);%11 for single hotspot        
    else
        id = sub2ind([20 20],11+distanceToHotSpot,11+distanceToHotSpot);%11 for single hotspot
    end
    
    net(i,j,k).driveI  = 0.5*net(i,j,k).driveI .* [adjMatForNearest12(id,1:400),zeros(1,100)]' + net(i,j,k).driveI;
        figure();
        pcolor(reshape(net(i,j,k).driveI(1:400),[2*sideNumI,2*sideNumI]))
        colorbar()
    
    net(i,j,k).run2DNet();
    net(i,j,k).detectDynamics();
    dPtIn = net(i,j,k).thetaPeak - Pt;
    dPgIn = net(i,j,k).gammaPeak - Pg;
    dFtIn = net(i,j,k).thetaFreq - Ft;
    dFgIn = net(i,j,k).gammaFreq - Fg;
    
    spikeTimes = net(i,j,k).spikeTimes';
    spikeTimes = spikeTimes(spikeTimes(:,2) <= net(i,j,k).numE ,:);
    spikeNum = zeros(net(i,j,k).numE,1);
    for neuronId = 1:net(i,j,k).numE
    spikeNum(neuronId) = ...
           sum(spikeTimes(:,2) == neuronId);
    end
    dFmeanIn = sum(adjMatForNearest12(id,1:400)' .* spikeNum)/sum(adjMatForNearest12(id,1:400)' .* rawSpikeNum) - 1;
    dFmean4allIn = mean(spikeNum-rawSpikeNum)/5;
    if computing == 0
        dPtOut = dPtIn;
        dPgOut = dPgIn;
        dFtOut = dFtIn;
        dFgOut = dFgIn;
        dPtIn = nan;
        dPgIn = nan;
        dFtIn = nan;
        dFgIn = nan;
        dFmeanOut = dFmeanIn;
        dFmeanIn = nan;
        dFmean4allOut = dFmean4allIn;
        dFmean4allIn = nan;
    else
        net(i,j,k).driveI(:) = 3; % set the drive back to 3
        if computing == 1
            distanceToHotSpot = 6;
            id = sub2ind([20 20],11+distanceToHotSpot,11+distanceToHotSpot);%11 for single hotspot
        else
            distanceToHotSpot = 10; %6 for outside for one hotspot
            id = sub2ind([20 20],7+distanceToHotSpot,7);%11 for single hotspot
        end
        net(i,j,k).driveI  = 0.5*net(i,j,k).driveI .* [adjMatForNearest12(id,1:400),zeros(1,100)]' + net(i,j,k).driveI;   
        net(i,j,k).run2DNet();
        net(i,j,k).detectDynamics();
        spikeTimes = net(i,j,k).spikeTimes';
        spikeTimes = spikeTimes(spikeTimes(:,2) <= net(i,j,k).numE ,:);
        spikeNum = zeros(net(i,j,k).numE,1);
        for neuronId = 1:net(i,j,k).numE
        spikeNum(neuronId) = ...
               sum(spikeTimes(:,2) == neuronId);
        end
        dFmeanOut = sum(adjMatForNearest12(id,1:400)' .* spikeNum)/sum(adjMatForNearest12(id,1:400)' .* rawSpikeNum) - 1;
        dFmean4allOut = mean(spikeNum-rawSpikeNum)/5;
        dPtOut = net(i,j,k).thetaPeak - Pt;
        dPgOut = net(i,j,k).gammaPeak - Pg;
        dFtOut = net(i,j,k).thetaFreq - Ft;
        dFgOut = net(i,j,k).gammaFreq - Fg;
    end    
    dPtMat(k,(computing+1),1) = dPtIn;
    dPgMat(k,(computing+1),1) = dPgIn;
    dFtMat(k,(computing+1),1) = dFtIn;
    dFgMat(k,(computing+1),1) = dFgIn;
    dFmeanMat(k,(computing+1),1) = dFmeanIn;
    dFmean4allMat(k,(computing+1),1) = dFmean4allIn;
    dPtMat(k,(computing+1),2) = dPtOut;
    dPgMat(k,(computing+1),2) = dPgOut;
    dFtMat(k,(computing+1),2) = dFtOut;
    dFgMat(k,(computing+1),2) = dFgOut;    
    dFmeanMat(k,(computing+1),2) = dFmeanOut;
    dFmean4allMat(k,(computing+1),2) = dFmean4allOut;
end
end

save externalDrive
%% Figure-7E 
load('externalDrive.mat')
for ii = 1:6
    if ii == 1
        usingMat = dPtMat;
        usingName = 'dPt';
    elseif ii == 2
        usingMat = dPgMat;
        usingName = 'dPg';
    elseif ii == 3
        usingMat = dFtMat;
        usingName = 'dFt';
    elseif ii == 4
        usingMat = dFgMat;
        usingName = 'dFg';
    elseif ii == 5
        usingMat = dFmeanMat;
        usingName = 'dFreq';
    else
        usingMat = dFmean4allMat;
        usingName = 'dFmean4All';
    end
model_series = squeeze(mean(usingMat(:,2:3,:)-usingMat(:,1,2)));
model_error = 0.5 *squeeze(std(usingMat(:,2:3,:)-usingMat(:,1,2)));
fig = figure('DefaultAxesFontSize',20);
b = bar(model_series,'grouped');
b(1).FaceColor = [1.0 1.0 .0];
b(2).FaceColor = [.0 .0 .9];
ylabel(usingName)
legend('Inside','Outside')
ylim([-20,60])
hold on

% Find the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);

% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none','markersize',28,'HandleVisibility','off');
end
xticklabels({'Single','Double'})
fig.OuterPosition(4) = fig.OuterPosition(4)/2;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,['fig_' usingName],'-dpdf','-r800')
end

