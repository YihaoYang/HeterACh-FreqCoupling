classdef EI2DNet < NeuralNet
    % E-I two layered Network Object created by Yihao Yang Apr 2019
    properties
        tau_z = 75.0
        numI
        sideNumI
        numE
        sideNumE
        numG %gKs pixel number
        sideNumG
        gKsMatrix
        dffsConst = 0.01
        decayConst = 0.05
        releaseNum = 10
        releaseAmp
        releaseRadius
        releaseDuration
        sqrtE2Iratio = 2
        sqrtG2Iratio = 5 % 10 fold to make a better resolution?
        
        I2IWeight = 1
        E2EWeight = 1
        I2EWeight = 1
        E2IWeight = 1
        relativeStrength
        thetaPeak
        gammaPeak
        thetaFreq
        gammaFreq
        gammaPeakAbsVal
        thetaMagnitude
        gammaMagnitude
        neuronRhythmVec
        
        halfDistance
        ruggedness
    end
    
    methods   
        function [x,traces] = extractSpikeTraces(self)                  
            
            cutoffTime = 1000; %ms
            
            % ignore spikes in initial transient for time
%             spiketimes = spiketimes(spiketimes(:,1)>cutoffTime,:);
            % get the traces
            x = cutoffTime:self.deltaT:self.stopT;%temp x points
            traces = self.vPoints(1:self.numE,self.tPoints>=cutoffTime);
            
        end
        
        function [x,traces] = computeSpikeTraces(self)                  
            spiketimes = self.spikeTimes';
            cutoffTime = 1000; %ms
            sigma = 1.5;
            % ignore spikes in initial transient for time
            spiketimes = spiketimes(spiketimes(:,1)>cutoffTime,:);
            % get the traces
            x = cutoffTime:self.deltaT:self.stopT;%temp x points
            traces = zeros(self.numE,length(x));
            for i = 1:self.numE
                mu = spiketimes(spiketimes(:,2)==i,1);                
                if ~isempty(mu) 
                    if length(mu) == 1
                        traces(i,:) = (normpdf(x,mu,sigma));
                    else
                        traces(i,:) = sum(normpdf(x,mu,sigma));
                    end
                end
%                 figure(666)
%                 subplot(2,1,1)
%                 hold on
%                 plot(x,traces(i,:))
            end            
        end
        
        function mapRhythmVec1(self)
            cutTime = self.stopT/5;
            endTime = self.stopT;
            freqUpperLimit = 100;
            step = 0.5*1000/freqUpperLimit;%ms
            
            thresholdGamma = 0.5;
            thresholdTheta = 12;
            strictness = 5;     % was 0.3       
            spiketimes = self.spikeTimes';
            % ignore spikes in initial transient for time
            spiketimes = spiketimes(spiketimes(:,1)>cutTime,:);

            % get the traces
            x = cutTime:step:endTime;%temp x points
            
            traces = zeros(self.neuronNum,length(x));  
            
            for i = 1:self.neuronNum
                mu = spiketimes(spiketimes(:,2)==i,1);
                ISIs = diff(mu);
                if length(ISIs)>= 3*(endTime - cutTime)/1000 % >3Hz
                    sigma = ([ISIs(1);ISIs] + [ISIs;ISIs(end)]) ...
                        /(2*strictness);
                    y = normpdf(x,mu,sigma); %temp y points
                    y = sum(y./max(y,[],2));
                    traces(i,:) = y/max(y);
                end
%                 figure(666)
%                 subplot(2,1,1)
%                 hold on
%                 plot(x,traces(i,:))
            end
            T = step/1000; %s
            Fs = 1/T; %Hz
            Y = fft(traces,[],2);    
            signalLength = size(traces,2);
            P2 = abs(Y/signalLength);    % two-sided spectrum P2
            P1 = P2(:,1:(signalLength/2+1)); % single-sided spectrum P1
            P1(:,2:end-1) = 2*P1(:,2:end-1);
            f = Fs*(0:(signalLength/2))/signalLength;

            % output
            baseline = mean(P1(:,2:end),2);
            singleThetaPeak = max(P1(:,f>4 & f<20),[],2) ./ baseline - 1;
            singleGammaPeak = max(P1(:,f>40 & f<70),[],2) ./ baseline - 1;
            self.neuronRhythmVec = zeros(self.neuronNum,1); % 0 for NULL type
            self.neuronRhythmVec(singleThetaPeak > thresholdTheta) = 1;      % 1 for theta 
            self.neuronRhythmVec(singleGammaPeak > thresholdGamma) = 2;      % 2 for Gamma, 3 for mixed
            self.neuronRhythmVec(singleThetaPeak > thresholdTheta & singleGammaPeak > thresholdGamma) = 3;
        end        
        
        function mapRhythmVec(self)
            thresholdTheta = 1;
            thresholdGamma = 1;
            numUsedHere = self.neuronNum;
            cutTime = self.stopT/5;
            endTime = self.stopT;
            freqUpperLimit = 100;
            step = 0.5*1000/freqUpperLimit;%ms
            
            tempSpikeTimes = self.spikeTimes;
            tempSpikeTimes = tempSpikeTimes(:,tempSpikeTimes(1,:)>cutTime);
            if numUsedHere == self.numE
                % only considering E-cells
                tempSpikeTimes = ...
                    tempSpikeTimes(:,tempSpikeTimes(2,:)<=self.numE);
            elseif numUsedHere == self.numI
                % only considering I-cells
                tempSpikeTimes = ...
                    tempSpikeTimes(:,tempSpikeTimes(2,:)>self.numE);
                tempSpikeTimes(2,:) = tempSpikeTimes(2,:) - self.numE;
            end         
            
            timePoints = cutTime:step:endTime;
            iterationTimes = length(timePoints);
            
            logicalRasterMat = zeros(iterationTimes,numUsedHere);
            for i = 1:iterationTimes
                temp1 = [abs(tempSpikeTimes(1,:)-timePoints(i)); ...
                    tempSpikeTimes(2,:)];
                temp2 = temp1(2,temp1(1,:) <= step/2);
                logicalRasterMat(i,temp2) = 1;
            end
            logicalRasterMat = logicalRasterMat';
            
            halfL = iterationTimes - 1;
            tempCorr = zeros(numUsedHere,halfL);
            for ii = 1:numUsedHere
                c = xcorr(logicalRasterMat(ii,:));
                temp3 = c(1:halfL) + c((end + 1 - halfL):end);
                tempCorr(ii,:) = temp3;
            end
         
            T = step/1000; %s
            Fs = 1/T; %Hz
            Y = fft(tempCorr,[],2);
            signalLength = halfL; % length of the signal (phase corr)
            P2 = abs(Y/signalLength);    % two-sided spectrum P2
            P1 = P2(:,1:signalLength/2+1); % single-sided spectrum P1
            P1(:,2:end-1) = 2*P1(:,2:end-1);
            f = Fs*(0:(signalLength/2))/signalLength;
            
            baseline = 1;
%             baseline = mean(P1(:,2:end),2);
            singleThetaPeak = max(P1(:,f>2.5 & f<20),[],2) ./ baseline;
            singleGammaPeak = max(P1(:,f>25 & f<100),[],2) ./ baseline;
            self.neuronRhythmVec = zeros(self.neuronNum,1); % 0 for NULL type
            self.neuronRhythmVec(singleThetaPeak > thresholdTheta) = 1;      % 1 for theta 
            self.neuronRhythmVec(singleGammaPeak > thresholdGamma) = 2;      % 2 for Gamma, 3 for mixed
            self.neuronRhythmVec(singleThetaPeak > thresholdTheta & singleGammaPeak > thresholdGamma) = 3;
        end
        
        function detectDynamics(self)
            numUsedHere = self.neuronNum;
            cutTime = self.stopT/5;
            endTime = self.stopT;
            freqUpperLimit = 100;
            step = 0.5*1000/freqUpperLimit;%ms
            
            tempSpikeTimes = self.spikeTimes;
            tempSpikeTimes = tempSpikeTimes(:,tempSpikeTimes(1,:)>cutTime);
            if numUsedHere == self.numE
                % only considering E-cells
                tempSpikeTimes = ...
                    tempSpikeTimes(:,tempSpikeTimes(2,:)<=self.numE);
            elseif numUsedHere == self.numI
                % only considering I-cells
                tempSpikeTimes = ...
                    tempSpikeTimes(:,tempSpikeTimes(2,:)>self.numE);
                tempSpikeTimes(2,:) = tempSpikeTimes(2,:) - self.numE;
            end         
            
            timePoints = cutTime:step:endTime;
            iterationTimes = length(timePoints);
            
            logicalRasterMat = zeros(iterationTimes,numUsedHere);
            for i = 1:iterationTimes
                temp1 = [abs(tempSpikeTimes(1,:)-timePoints(i)); ...
                    tempSpikeTimes(2,:)];
                temp2 = temp1(2,temp1(1,:) <= step/2);
                logicalRasterMat(i,temp2) = 1;
            end
            logicalRasterMat = logicalRasterMat';
            
            halfL = iterationTimes - 1;
            tempCorr = zeros(numUsedHere,halfL);
            for ii = 1:numUsedHere
                c = xcorr(logicalRasterMat(ii,:));
                temp3 = c(1:halfL) + c((end + 1 - halfL):end);
                tempCorr(ii,:) = temp3;
            end
            meanXCorr = mean(tempCorr);
            T = step/1000; %s
            Fs = 1/T; %Hz
            Y = fft(meanXCorr);
            signalLength = halfL; % length of the signal (phase corr)
            P2 = abs(Y/signalLength);    % two-sided spectrum P2
            P1 = P2(1:signalLength/2+1); % single-sided spectrum P1
            P1(2:end-1) = 2*P1(2:end-1);
            f = Fs*(0:(signalLength/2))/signalLength;
            
%             figure;
%             lagTimes = step*(1:halfL);
%             plot(lagTimes,meanXCorr);
%             xlabel('lag - ms')
%             title('phase correlation')
%             
%             figure;
            plot(f(2:end),P1(2:end),'LineWidth',2)
            xlabel('f (Hz)')
            ylabel('|P1(f)|')

            % output
%             baseline = mean(P1(2:end));
            [argVal, argMax] = max(P1(f>2.5 & f<20));
            self.thetaPeak = argVal;
            temp = f(f>2.5 & f<20);
            self.thetaFreq = temp(argMax);
            [argVal, argMax] = max(P1(f>25 & f<100));
            self.gammaPeak = argVal;
            temp = f(f>25 & f<100);
            self.gammaFreq = temp(argMax);
            self.relativeStrength = (self.thetaPeak) / (self.gammaPeak);
        end
        
        function run2DNet(self)
            % random initilization %%%%%%
            self.states = rand(4,self.neuronNum);
            self.states(4,:) = -70+40*self.states(4,:);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            startT = self.startT;
            stopT = self.stopT;
            deltaT = self.deltaT;
            tPoints = self.tPoints;
            neuronNum = self.neuronNum;
            neuronIdVec = (1:neuronNum)';
            gKs = self.gKs;
            adjMatrix = self.adjacencyMatrix;
            channelZ = self.states(1,:)';
            channelH = self.states(2,:)';
            channelN = self.states(3,:)';
            memV = self.states(4,:)';
            tau = self.tau;
            synE = self.synE;
            DC = self.driveI(:);
            THRESHOLD_AP = -20;
            C = 1; %uf/cm2
            v_Na = 55.0; %mV
            v_K = -90; %mV
            v_L = -60; %mV
            g_Na = 24; %mS/cm2
            g_Kdr = 3.0; %mS/cm2
            g_L = 0.02; %mS/cm2
            tau_z = self.tau_z;
            
            spikeTimes = zeros(neuronNum,stopT-startT);
           
            spikeCounts = zeros(neuronNum,1);
            
%             vPoints = zeros(neuronNum, length(tPoints));
            
            %%%%%%%%%%%%%%%%%%%% add poisson kicks %%%%%%%%%%%%%%%%%%%%%%%
            poissonRate = 1/150; 
            poissonKickAmp = 0;%6;
            poissonKickDuration = 1; %ms
            noiseI=0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for i = 1:length(tPoints)
                t = tPoints(i);
                if ~mod(t,poissonKickDuration)
                    noiseI = poissonKickAmp*(rand(neuronNum,1)<poissonRate);                    
                end
                vPoints(:,i) = memV;
                % determine synI vector
                [synEMat, memVMat] = meshgrid(synE,memV);
                ithLatestSpike = 1;
                tempIdx = neuronIdVec + ...
                    neuronNum*(spikeCounts - ithLatestSpike);
                tempIdx(tempIdx<=0) = 1;   
                synI=0;
                while any((t-spikeTimes(tempIdx)) < 20) && any(ithLatestSpike<(spikeCounts))
                synI = synI + adjMatrix .* (memVMat - synEMat) * (logical((spikeCounts+1-ithLatestSpike)>0) .* ...
                       exp((spikeTimes(tempIdx)-t)./ tau));
                ithLatestSpike = ithLatestSpike + 1;
                tempIdx = neuronIdVec + ...
                    neuronNum*(spikeCounts - ithLatestSpike);
                tempIdx(tempIdx<=0) = 1;                
                end               
                   
                isFiring = (memV < THRESHOLD_AP);
                totalI = noiseI + DC - synI;
                
                % RK4 method starts
                kV = [memV,memV,memV,memV];
                kZ = [channelZ,channelZ,channelZ,channelZ];
                kH = [channelH,channelH,channelH,channelH];
                kN = [channelN,channelN,channelN,channelN];
                for colInd = 1:4
                    mInf = 1./(1 + exp((-kV(:,colInd)-30.0)/9.5));
                    hInf = 1./(1 + exp((kV(:,colInd)+53.0)/7.0));
                    nInf = 1./(1 + exp((-kV(:,colInd)-30.0)/10));
                    zInf = 1./(1 + exp((-kV(:,colInd)-39.0)/5.0));
                    hTau = 0.37 + 2.78./(1 + exp((kV(:,colInd)+40.5)/6));
                    nTau = 0.37 + 1.85./(1 + exp((kV(:,colInd)+27.0)/15));
                    fh = (hInf - kH(:,colInd))./hTau;
                    fn = (nInf - kN(:,colInd))./nTau;
                    fz = (zInf - kZ(:,colInd))./tau_z;
                    fv = (1/C)*(g_Na*(mInf.^3).*kH(:,colInd).* ...
                        (v_Na-kV(:,colInd))+ ...
                        g_Kdr*(kN(:,colInd).^4).*(v_K-kV(:,colInd))+ ...
                        gKs.*kZ(:,colInd).*(v_K - kV(:,colInd))+ ...
                        g_L*(v_L-kV(:,colInd))+ totalI);
                    kH(:,colInd) = deltaT*fh;
                    kN(:,colInd) = deltaT*fn;
                    kZ(:,colInd) = deltaT*fz;
                    kV(:,colInd) = deltaT*fv;
                    if (colInd == 1 || colInd == 2)
                        kH(:,colInd+1) = kH(:,colInd+1)+0.5*kH(:,colInd);
                        kN(:,colInd+1) = kN(:,colInd+1)+0.5*kN(:,colInd);
                        kZ(:,colInd+1) = kZ(:,colInd+1)+0.5*kZ(:,colInd);
                        kV(:,colInd+1) = kV(:,colInd+1)+0.5*kV(:,colInd);
                    elseif (colInd == 3)
                        kH(:,colInd+1) = kH(:,colInd+1)+kH(:,colInd);
                        kN(:,colInd+1) = kN(:,colInd+1)+kN(:,colInd);
                        kZ(:,colInd+1) = kZ(:,colInd+1)+kZ(:,colInd);
                        kV(:,colInd+1) = kV(:,colInd+1)+kV(:,colInd);
                    end
                end
                memV =     memV +     (kV(:,1) + 2 * kV(:,2) + 2 * kV(:,3) + kV(:,4))/6.0;
                channelH = channelH + (kH(:,1) + 2 * kH(:,2) + 2 * kH(:,3) + kH(:,4))/6.0;
                channelN = channelN + (kN(:,1) + 2 * kN(:,2) + 2 * kN(:,3) + kN(:,4))/6.0;
                channelZ = channelZ + (kZ(:,1) + 2 * kZ(:,2) + 2 * kZ(:,3) + kZ(:,4))/6.0;
                % RK4 ends
                isFiring = isFiring & (memV >= THRESHOLD_AP);
                spikeCounts = spikeCounts + isFiring;
                tempIdx = neuronIdVec(isFiring) + ...
                    neuronNum*(spikeCounts(isFiring) - 1);                
                spikeTimes(tempIdx) = t;
            end
            self.spikeTimes = [];
            for i = 1:neuronNum
                self.spikeTimes = [self.spikeTimes, ...
                    [spikeTimes(i,1:spikeCounts(i)); ...
                    i*ones(1,spikeCounts(i))]];
            end
            self.vPoints = vPoints;
        end
        
        function setBasics(self,inSideNumI,gKsCorrCutoff)
            self.sideNumI = inSideNumI;
            self.sideNumE = self.sqrtE2Iratio * self.sideNumI;
            self.sideNumG = self.sqrtG2Iratio * self.sideNumI;
            self.numI = self.sideNumI^2;
            self.numE = self.sideNumE^2;
            self.numG = self.sideNumG^2;
            self.neuronNum = self.numI + self.numE;
            
            % mV   E_syn = 0 for excitatory synapses;
            % -75mV for inhibitory synapses
            self.synE = [zeros(self.numE,1); -75*ones(self.numI,1)];
            % tauE = 3 ms tauI = 4 ms
            self.tau = [3*ones(self.numE,1); 3*ones(self.numI,1)];
            
            % build gKs concentration and assign gKs data to neurons
            self.map2DgKsWithReleaseAndDecay(gKsCorrCutoff);
                        
            % placeholder for adjacency matrix
            self.adjacencyMatrix = zeros(500);
            
            DC = 3;      
            DCspread = 0;  
            ACamp = 0;
            ACfreq = 0;
            self.states = rand(4,self.neuronNum);
            self.states(4,:) = -70+40*self.states(4,:);
            self.driveI = DC + DCspread*(rand(1,self.neuronNum)-0.5);
            self.noiseI = zeros(self.neuronNum,1);
            self.ACamp = ACamp;
            self.ACfreq = ACfreq;
            self.rearrangedIndex = 1:self.neuronNum;
            self.sortedIndex = 1:self.neuronNum;
            self.tPoints = self.startT:self.deltaT:self.stopT;
        end
        
        function map2DgKsWithReleaseAndDecay(self,gKsCorrCutoff)
            gKsTemp = 0.2 + 1.3*rand(self.numE,1);
            self.gKsMatrix = build2DAdjMatrix( ...
                self.sideNumE, ...
                1, ...
                self.dffsConst ...
                );
%             gKsTemp = ones(self.numG,1);
%             figure;
%             heatmapMatrix = ...
%                 reshape(gKsTemp,[self.sideNumE,self.sideNumE]);
%             pcolor(heatmapMatrix);pcolor(heatmapMatrix);caxis([0,1.5]);colorbar;
%             title('t = 0');
%             axis tight manual
%             set(gca,'nextplot','replacechildren');
%             v = VideoWriter('gksMap.avi');
%             open(v);
            
%             releaseProbTwist = (1:self.numG)';
%             tempX = mod(releaseProbTwist,self.sideNumG);
%             tempY = ceil(releaseProbTwist/self.sideNumG);
%             tempR = sqrt((tempX-self.sideNumG/2).^2 + (tempY-self.sideNumG/2).^2);
%             [~,tempSortedId] = sort(tempR);
%             releaseProbTwist(tempSortedId) = 10.^linspace(0,0,self.numG);
            
            %gKsLimits = [self.gKsMin self.gKsMax];
            for t = 1:gKsCorrCutoff
%                 heatmapMatrix = reshape(gKsTemp, ...
%                    [self.sideNumE,self.sideNumE]);
%                 figure(6666);pcolor(heatmapMatrix);%caxis(gKsLimits);
%                 colorbar;
%                 title(num2str(t));
                %gKsTemp(releaseProbTwist.*rand(self.numG,1) < self.releaseNum/self.numG) = 0.2;
                if (1 == mod (t,self.releaseDuration))
                    releaseSourceId = randsample(self.numE, self.releaseNum);
                    tempMap = zeros(self.numE,1);
                    if ~isempty(self.halfDistance)
                        releaseSourceId = [sub2ind([20 20],(11-self.halfDistance),(11-self.halfDistance)); ...
                                               sub2ind([20 20],(10+self.halfDistance),(10+self.halfDistance))];
                    end
                    if ~isempty(releaseSourceId)
                        tempMap(releaseSourceId) = self.releaseAmp;
                        releaseDffsConst = 0.12;
                        for tempT = 1:self.releaseRadius
                            tempMap = tempMap + releaseDffsConst.*logical(self.gKsMatrix) * tempMap - ...
                                (releaseDffsConst*sum(logical(self.gKsMatrix)).*tempMap')';
                        end
                    end
                end
                gKsTemp = gKsTemp - tempMap + ...
                    self.decayConst*(self.gKsMax - gKsTemp) + ...
                                    self.gKsMatrix * gKsTemp - ...
                                    (sum(self.gKsMatrix).*gKsTemp')';
                
%                 heatmapMatrix = ...
%                     reshape(gKsTemp,[self.sideNumE,self.sideNumE]);
%                 
%                 pcolor(heatmapMatrix);caxis([0.2,1.5]);colorbar;
%                 title(['t = ',num2str(t)]);
%                 frame = getframe(gcf);
%                 writeVideo(v,frame);
            gKsTemp(gKsTemp<self.gKsMin) = self.gKsMin;
            gKsTemp(gKsTemp>self.gKsMax) = self.gKsMax;
            end
%                         close(v);
            %             gKsTemp(gKsTemp<0.2) = 0.2;
            %             gKsTemp(gKsTemp>1.5) = 1.5;
            gKsTempMat = reshape(gKsTemp,[self.sideNumE,self.sideNumE]);
            gKsTempForICells = 0.5* ...
                (gKsTempMat(1:2:end,:) + gKsTempMat(2:2:end,:));
            gKsTempForICells = 0.5* ...
                (gKsTempForICells(:,1:2:end) + gKsTempForICells(:,2:2:end));         
            self.gKs = [gKsTemp;gKsTempForICells(:)];
        end
        
    end
end

