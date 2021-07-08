%% Simulations for double-peaked gks distribution
close all;clear all;
load('adjMat.mat')
numX = 8; % number of radius points
numY = 5; % number of distance points
numRepetition = 4;
sideNumI = 10;
radiusPoints = round(linspace(5,25,numX));
halfDistancePoints = 1:numY; % control the distance between two spots

net(numX,numY,numRepetition) = EI2DNet;

parfor k = 1:numRepetition
for i = 1:numX
    for j = 1:numY
        net(i,j,k).synapticStrength = 1; % this is a multiplier for synaptic weights
        net(i,j,k).I2IWeight = 0.04;
        net(i,j,k).E2IWeight = 0.05;
        net(i,j,k).E2EWeight = 0.01;
        net(i,j,k).I2EWeight = 0.04;
        net(i,j,k).gKsMin = 0.2;
        
        net(i,j,k).dffsConst = 0.1;
        net(i,j,k).decayConst = 0.02;
        net(i,j,k).releaseNum = 2;
        net(i,j,k).halfDistance = halfDistancePoints(j);
        net(i,j,k).releaseRadius = radiusPoints(i);
        net(i,j,k).releaseAmp = 6/net(i,j,k).releaseNum;
        net(i,j,k).releaseDuration = 1000;
        net(i,j,k).setBasics(sideNumI,1000);        
        net(i,j,k).adjacencyMatrix = adjMat;
        
        net(i,j,k).run2DNet;        
    end
end
end

% uncomment the next line to save the simulation data
% save double_peaked_gks_data.mat
%% Remapping radius parameters
% this section is to map original radius parameters to displayed radius
% parameters since the diffusion in the gks map generation caused
% distortion.
clear all
radiusPoints = [1 3 5 7 8 9 11 13 14 15 16 19 22 25];
numX = length(radiusPoints);
numY = 1;
numRepetition = 1;
sideNumI = 10;
gKsLowerBoundPoints = linspace(0,1.4,numX);
inputDCPoints = linspace(1,6,numY);

net(numX,numY,numRepetition) = EI2DNet;
remappedR = zeros(numX,numY,numRepetition);

for k = 1:numRepetition
for i = 1:numX
    for j = 1:numY
        net(i,j,k).synapticStrength = 1; % this is a multiplier for synaptic weights
        net(i,j,k).I2IWeight = 0.04;
        net(i,j,k).E2IWeight = 0.05;
        net(i,j,k).E2EWeight = 0.01;
        net(i,j,k).I2EWeight = 0.04;
        net(i,j,k).gKsMin = 0.2;
        
        net(i,j,k).dffsConst = 0.1;
        net(i,j,k).decayConst = 0.02;
        net(i,j,k).releaseNum = 1;

        net(i,j,k).releaseRadius = radiusPoints(i);
        net(i,j,k).releaseAmp = 6/net(i,j,k).releaseNum;
        net(i,j,k).releaseDuration = 1000;
        net(i,j,k).setBasics(sideNumI,1000);
        tempGks = net(i,j,k).gKs(1:4*sideNumI^2);
        remappedR(i) = sqrt(sum(1.5-tempGks)/(1.5*pi));
    end
end
end

%% Simulations for randomly generated spatial gks distributions
close all;clear all;
load('adjMat.mat')
numX = 8; % number of radius points
numY = 5; % number of spot variables
numRepetition = 4;
sideNumI = 10;
radiusPoints = round(linspace(1,15,numX));
spotNumPoints = round(linspace(1,20,numY));
net(numX,numY,numRepetition) = EI2DNet;

thetaMat = zeros(numX,numY,numRepetition);
gammaMat = zeros(numX,numY,numRepetition);
relativeStrengthMat = zeros(numX,numY,numRepetition);

parfor k = 1:numRepetition
for i = 1:numX
    for j = 1:numY
        net(i,j,k).synapticStrength = 1; % this is a multiplier for synaptic weights
        net(i,j,k).I2IWeight = 0.04;
        net(i,j,k).E2IWeight = 0.05;
        net(i,j,k).E2EWeight = 0.01;
        net(i,j,k).I2EWeight = 0.04;
        net(i,j,k).gKsMin = 0.2;
        
        net(i,j,k).dffsConst = 0.1;
        net(i,j,k).decayConst = 0.02;
        net(i,j,k).releaseNum = spotNumPoints(j);
        net(i,j,k).releaseRadius = radiusPoints(i);
        net(i,j,k).releaseAmp = 6/net(i,j,k).releaseNum;
        net(i,j,k).releaseDuration = 1000;
        net(i,j,k).setBasics(sideNumI,1000);
        
        net(i,j,k).adjacencyMatrix = adjMat;
        net(i,j,k).run2DNet;
        net(i,j,k).detectDynamics;        
        thetaMat(i,j,k) = net(i,j,k).thetaPeak;
        gammaMat(i,j,k) = net(i,j,k).gammaPeak;
        relativeStrengthMat(i,j,k) = net(i,j,k).relativeStrength;
    end
end
end

% uncomment the next line to save data file
% save random_gks_data.mat;