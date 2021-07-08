classdef NeuralNet < handle
    % Neural Network Object created by Yihao Yang Feb 2019
    methods
        function clusterBy(self,byWhat,cutoffVal)
            if ~isempty(self.clusterIndex)
                disp("Warning - you have already clustered the index!")
            end
            if byWhat == "MPC"
                mpcCutoff = cutoffVal;
                distByMpc = 1 - self.pairWiseMpc; 
                %distByMpc = tril(distByMpc)%triu(distByMpc);
                distByMpc = (distByMpc + distByMpc')/2;
                distByMpc = squareform(distByMpc);
                linkageByMpc = linkage(distByMpc,'complete'); %Farthest distance 'complete'               
                self.clusterIndex = cluster(linkageByMpc, ...
                    'cutoff',1-mpcCutoff,'criterion','distance');
                figure;
                D = dendrogram(linkageByMpc,'ColorThreshold',1-mpcCutoff);
                set(D,'LineWidth',2);
                ylim([0,1]);
            end
            self.rearrangeBy("cluster");            
        end
        
        function rearrangeBy(self,byWhat)            
%             if self.rearrangedIndex ~= 1:self.neuronNum
%                 disp("Warning - you have already rearranged the index!");
%             end
            if byWhat == "DC"
                [~,sortedID] = sort(self.driveI);                
            elseif byWhat == "gKs"
                [~,sortedID] = sort(self.gKs);
            elseif byWhat == "cluster"
                [~,sortedID] = sort(self.clusterIndex);
            end
            self.sortedIndex = sortedID;
            self.rearrangedIndex = sortedID; %consistent num
            for i = 1:length(sortedID)
                self.rearrangedIndex(i) = find(sortedID == i);
            end            
        end

        function [traces, traces_all] = computeNormSpikeTraces(self,inStrictness)
            strictness = inStrictness*50;            
            spiketimes = self.spikeTimes';
            % ignore spikes in initial transient for time
            spiketimes = spiketimes(spiketimes(:,1)>self.stopT/2,:);
            % get the traces
            x = (self.stopT/2):self.deltaT:self.stopT;%temp x points
            traces = zeros(self.neuronNum,length(x));
            for i = 1:self.neuronNum
                mu = spiketimes(spiketimes(:,2)==i,1);
                ISIs = diff(mu);
                if ~isempty(ISIs)
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
            traces_all = mean(traces);
%             subplot(2,1,2)
%             plot(x,traces_all)
%             figure(666); hold off;
        end
       
        function showColorCodedRasterPlot(self)           
%             figure('Name','spiking raster plot','NumberTitle','off');
            hold all
            if isempty(self.clusterIndex)  
                for neuronId = 1:self.neuronNum
                    xPoints = self.spikeTimes(1, ...
                        self.spikeTimes(2,:) == neuronId);
                    scatter(xPoints, ...
                        self.spikeTimes(2,self.spikeTimes(2,:) == neuronId), ...
                        10,(self.gKs(neuronId)).*ones(size(xPoints')),'filled')
                end
            else  
                groupID = self.clusterIndex(self.spikeTimes(2,:));
                for i = min(self.clusterIndex):max(self.clusterIndex)
                    neuronIndexPoints = self.rearrangedIndex( ...
                        self.spikeTimes(2,(groupID == i)));
                    plot(self.spikeTimes(1,(groupID == i)), ...
                        neuronIndexPoints,'.');
                end
            end
%             title(['Syn Weight = ', num2str(self.synapticStrength)]);

            ylim([1,self.neuronNum]);
            xlim([4500 5000]);
        end
        
        function setProperties(self, ...
                neuronNum, ...
                synapticStrength, ...
                gKs, ...
                DC, ...
                DCspread, ...
                ACamp, ...
                ACfreq, ...
                neighborNum, ...
                rewiringProb)
            
            self.neuronNum = neuronNum;
            self.rearrangedIndex = 1:neuronNum;
            self.sortedIndex = 1:neuronNum;
            %default orderd the same as original

            self.synapticStrength = synapticStrength;
            self.gKs = gKs;
            self.tPoints = self.startT:self.deltaT:self.stopT;
            
            if nargin == 10
                self.neighborNum = neighborNum;
                self.rewiringProb = rewiringProb;
                self.generateAdjacencyMatrix;
            elseif nargin == 8
                self.generateAdjacencyMatrix("2D");
            end
            %self.matrix2List;
            
            self.states = rand(4,neuronNum);
            self.states(4,:) = -70+40*self.states(4,:);
            self.driveI = DC + DCspread*(rand(1,neuronNum)-0.5);
            %self.driveI = linspace(DC-DCspread/2,DC+DCspread/2,neuronNum);
            self.noiseI = zeros(neuronNum,1);
            self.ACamp = ACamp;
            self.ACfreq = ACfreq;            
        end

        function uncluster(self)
            self.clusterIndex = [];
        end
        
        function unrearrange(self)
            self.rearrangedIndex = 1:self.neuronNum;
            self.sortedIndex = 1:self.neuronNum;
        end        
    end
    
    properties%properties
        neuronNum
        neighborNum
        rewiringProb
        adjacencyMatrix
        synapticStrength = 0.01
        
        startT = 0
        stopT = 5000%ms
        deltaT = 0.05
        tPoints
        
        rearrangedIndex
        sortedIndex
        clusterIndex        
        ACamp
        ACfreq
        MPC
        pairWiseMpc
        GBM
    end    
    properties%neuron properties
        gKs
        gKsMax = 1.5
        gKsMin = 0.2
        states
        driveI
        noiseI
        
        tau = 0.5%ms
        synE = 0%mV   E_syn = 0 for excitatory synapses; -75mV for inhibitory synapses
        spikeTimes
        vPoints
    end  
    properties%visualization properties
        fontSize = 16
    end
end