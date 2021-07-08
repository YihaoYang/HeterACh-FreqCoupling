#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 11:00:07 2020

@author: yihao yang
"""

import matplotlib.pyplot as plt
from scipy.spatial import distance
from scipy import signal
import numpy as np

    
# Constants
DEFAULT_NEURONUM = 500
DEFAULT_TEND = 7000
DEFAULT_IDRIVE = 3

DEFAULT_XNUME = 20
DEFAULT_YNUME = 20
DEFAULT_XNUMI = 10
DEFAULT_YNUMI = 10
DEFAULT_DEGREE_EE = 40
DEFAULT_DEGREE_EI = 10
DEFAULT_DEGREE_IE = 400
DEFAULT_DEGREE_II = 100
DEFAULT_WEIGHT_EE = 0.01
DEFAULT_WEIGHT_EI = 0.05
DEFAULT_WEIGHT_IE = 0.04
DEFAULT_WEIGHT_II = 0.04

DEFAULT_TAU_SYN = 3
DEFAULT_GKS_MIN = 0.2
DEFAULT_GKS_MAX = 1.5

# Class
class NeuroNet():
    
    def __init__(self,
                 neuroNum = DEFAULT_NEURONUM,
                 tEnd = DEFAULT_TEND,
                 Idrive = DEFAULT_IDRIVE,
                 tauSyn = DEFAULT_TAU_SYN,
                 gKsMin = DEFAULT_GKS_MIN,
                 gKsMax = DEFAULT_GKS_MAX):
        '''
        

        Parameters
        ----------
        neuroNum : TYPE, optional
            DESCRIPTION. The default is DEFAULT_NEURONUM.
        tEnd : TYPE, optional
            DESCRIPTION. The default is DEFAULT_TEND.
        Idrive : TYPE, optional
            DESCRIPTION. The default is DEFAULT_IDRIVE.
        tauSyn : TYPE, optional
            DESCRIPTION. The default is DEFAULT_TAU_SYN.

        Returns
        -------
        None.

        '''
        # simulation properties
        self.tEnd = tEnd # ms
        self.tStep = 0.05 # ms
        self.tPoints = np.arange(0,self.tEnd,self.tStep)
        # ensemble properties
        self.neuroNum = neuroNum
        self.Idrive = Idrive*np.ones(shape=(self.neuroNum,1))
        # neuronal properties
        self.gKsMin = gKsMin
        self.gKsMax = gKsMax
        self.randomInitialStates()
        self.gKs = self.gKsMax
        # initial adjMat
        self.adjMat = np.zeros(shape=(self.neuroNum,self.neuroNum))
        self.Esyn = np.zeros((self.neuroNum,1)) 
        # 0 mV for excitatory synapses;
        # -75mV for inhibitory synapses
        self.tauSyn = DEFAULT_TAU_SYN*np.ones((self.neuroNum,1)) # ms
        
    def randomInitialStates(self):
        self.states = np.random.rand(self.neuroNum,4)
        self.states[:,3] = -70 + 40 * self.states[:,3]
        return self
    
    def zerolikeInitialStates(self,logV=False):
        originalDC = self.Idrive.copy()
        originalT = self.tEnd
        self.Idrive[:] = -1
        self.tEnd = 500
        self.tPoints = np.arange(0,self.tEnd,self.tStep) - self.tEnd        
        self.runSimulation(isNet = False,logV=logV)
        if logV: self.tPoints_before = self.tPoints.copy()
        self.Idrive = originalDC
        self.tEnd = originalT
        self.tPoints = np.arange(0,self.tEnd,self.tStep)
        return self
        
    def mexicanHat(self,
                 xNumE = DEFAULT_XNUME,
                 yNumE = DEFAULT_YNUME,
                 xNumI = DEFAULT_XNUMI,
                 yNumI = DEFAULT_YNUMI,
                 degreeEE = DEFAULT_DEGREE_EE,
                 degreeEI = DEFAULT_DEGREE_EI,
                 degreeIE = DEFAULT_DEGREE_IE,
                 degreeII = DEFAULT_DEGREE_II,
                 weightEE = DEFAULT_WEIGHT_EE,
                 weightEI = DEFAULT_WEIGHT_EI,
                 weightIE = DEFAULT_WEIGHT_IE,
                 weightII = DEFAULT_WEIGHT_II):
        '''
        

        Parameters
        ----------
        xNumE : TYPE, optional
            DESCRIPTION. The default is DEFAULT_XNUME.
        yNumE : TYPE, optional
            DESCRIPTION. The default is DEFAULT_YNUME.
        xNumI : TYPE, optional
            DESCRIPTION. The default is DEFAULT_XNUMI.
        yNumI : TYPE, optional
            DESCRIPTION. The default is DEFAULT_YNUMI.
        degreeEE : TYPE, optional
            DESCRIPTION. The default is DEFAULT_DEGREE_EE.
        degreeEI : TYPE, optional
            DESCRIPTION. The default is DEFAULT_DEGREE_EI.
        weightEE : TYPE, optional
            DESCRIPTION. The default is DEFAULT_WEIGHT_EE.
        weightEI : TYPE, optional
            DESCRIPTION. The default is DEFAULT_WEIGHT_EI.
        weightIE : TYPE, optional
            DESCRIPTION. The default is DEFAULT_WEIGHT_IE.
        weightII : TYPE, optional
            DESCRIPTION. The default is DEFAULT_WEIGHT_II.

        Returns
        -------
        None.

        '''
        
        self.numE = xNumE * yNumE
        self.xNumE,self.yNumE = xNumE,yNumE
        
        self.numI = self.neuroNum - self.numE
        self.xNumI,self.yNumI = xNumI,yNumI
        
        if self.numI != xNumI * yNumI:
            print('ERROR!!')
        self.Esyn[-self.numI:,:] = -75 # mV
        
        # assign x, y coordinates
        xLocE = np.arange(xNumE) + 0.5 # + 0.5 for periodic condition
        yLocE = np.arange(yNumE) + 0.5
        xLocE,yLocE = np.meshgrid(xLocE,yLocE)
        
        self.coordsE = np.stack((xLocE.reshape(-1),yLocE.reshape(-1))).T
        
        xLocI = (np.arange(xNumI) + 0.5) * (xNumE / xNumI)
        yLocI = (np.arange(yNumI) + 0.5) * (yNumE / yNumI)
        xLocI,yLocI = np.meshgrid(xLocI,yLocI)
        
        self.coordsI = np.stack((xLocI.reshape(-1),yLocI.reshape(-1))).T
        
        # compute mexican-hat adjacency matrix
        # compute distance matrices       
        distEE = distance.cdist(self.coordsE,self.coordsE,
                                lambda a,b: self.computeDist(a,b))
        distEI = distance.cdist(self.coordsI,self.coordsE,
                                lambda a,b: self.computeDist(a,b))
        self.distEE = distEE
        self.distEI = distEI
        
        # compute adjEE and adjEI 
        if degreeEE >= self.numE:
            adjMatEE = weightEE * np.ones(shape = (self.numE,self.numE))
        else:
            adjMatEE = np.zeros(shape = (self.numE,self.numE))
    
            adjMatEE[
                np.argsort(distEE,axis = 0)[1:degreeEE+1,:].T.reshape(-1),
                np.concatenate(
                    [i*np.ones(degreeEE,dtype=int) for i in np.arange(self.numE)])
                ] = weightEE
            
        if degreeEI >= self.numI:
            adjMatEI = weightEI * np.ones(shape = (self.numI,self.numE))
        else:
            adjMatEI = np.zeros(shape = (self.numI,self.numE))
            adjMatEI[
                np.argsort(distEI,axis = 0)[:degreeEI,:].T.reshape(-1),
                np.concatenate(
                    [i*np.ones(degreeEI,dtype=int) for i in np.arange(self.numE)])
                ] = weightEI
        
        # compute adjIE and adjII: all to all connection if degree < # of cells
        if degreeIE >= self.numE:
            adjMatIE = weightIE * np.ones(shape = (self.numE,self.numI))
        else:
            distIE = distance.cdist(self.coordsE,self.coordsI,
                                    lambda a,b: self.computeDist(a, b))
            adjMatIE = np.zeros(shape = (self.numE,self.numI))
            adjMatIE[
                np.argsort(distIE,axis=0)[:degreeIE,:].T.reshape(-1),
                np.concatenate(
                    [i*np.ones(degreeIE,dtype=int) for i in np.arange(self.numI)])
                ] = weightIE
            
        if degreeII >= self.numI:
            adjMatII = weightII * np.ones(shape = (self.numI,self.numI))
        else:
            distII = distance.cdist(self.coordsI,self.coordsI,
                                    lambda a,b: self.computeDist(a,b))
            adjMatII = np.zeros(shape = (self.numI,self.numI))

            adjMatII[
                np.argsort(distII,axis = 0)[1:degreeII+1,:].T.reshape(-1),
                np.concatenate(
                    [i*np.ones(degreeII,dtype=int) for i in np.arange(self.numI)])
                ] = weightII
        
        # finally get the adjMat
        self.adjMat = np.vstack((np.hstack((adjMatEE,adjMatIE)),
                       np.hstack((adjMatEI,adjMatII))))
        return self
    

    # compute the euclidean distance with periodic boundary conditions
    def computeDist(self,a,b):
        bounds = np.array([self.xNumE,self.yNumE])
        delta = np.abs(a-b)
        delta = np.where(delta > 0.5 * bounds,delta - bounds,delta)
        return np.sqrt((delta ** 2).sum(axis = -1)) 
    
    def mapGks(self, 
               r, 
               releaseLocs = np.array([[0.25,0.25],[0.75,0.75]]),
               sharpness = 2):
        '''
        Parameters
        ----------
        releaseLocs : TYPE, optional
            DESCRIPTION. The default is np.array([]). Normalized by x,y ranges.

        Returns
        -------
        None.

        '''        
        
        if releaseLocs.size>0:
            self.releaseR = r
            self.coordsRelease = np.array([self.xNumE,self.yNumE]) * releaseLocs      
        
            distER = (distance.cdist(self.coordsRelease,self.coordsE,
                                         lambda a,b: self.computeDist(a,b))
                                         .min(axis=0).reshape(-1,1))
            distIR = (distance.cdist(self.coordsRelease,self.coordsI,
                                         lambda a,b: self.computeDist(a,b))
                                         .min(axis=0).reshape(-1,1))  
            distToR = np.vstack((distER,distIR))
            self.distToR = distToR
            sigmoid = lambda x: 1/(1 + np.exp(-x))
            # self.sigmoidDistToR = sigmoidDistToR
            # sigmoidDistToR -= sigmoidDistToR.min()
            self.gKs = self.gKsMin + sigmoid(sharpness*(distToR - r)) * (
                self.gKsMax - self.gKsMin) 
            return self
        
    def runSimulation(self, 
                      isNet = True, 
                      isSTDP = False, 
                      silentSynapse = False,
                      externalInput = False,
                      ex_drive_strength = 0.1,
                      poisson_noise = False,
                      poisson_rate = 1/200,
                      poisson_amp = 6,
                      logV = False):        
        
        THRESHOLD_AP = -20 # mV
        C = 1 # uf/cm2
        v_Na = 55.0 # mV
        v_K = -90 # mV
        v_L = -60 # mV
        g_Na = 24 # mS/cm2
        g_Kdr = 3.0 # mS/cm2
        g_L = 0.02 # mS/cm2

        spikeTimes = np.zeros((self.neuroNum,self.tEnd))          
        spikeCounts = np.zeros((self.neuroNum,1),dtype=int)     
        # vPoints = np.zeros(size(tPoints));

        channelZ = self.states[:,[0]]
        channelH = self.states[:,[1]]
        channelN = self.states[:,[2]]
        memV = self.states[:,[3]]
        
        if logV: 
            logCounter = 0
            self.vPoints = np.zeros(shape=(self.neuroNum,self.tPoints.size))
            # temp current logger
            self.iPoints = np.zeros(shape=(self.neuroNum,self.tPoints.size))
        
        colIdx = np.arange(4)
        neuroIdx = np.arange(self.neuroNum).reshape(-1,1)
        Itotal = self.Idrive
        STDPon = False
        STDPoff = False
        windowIsyn = 20 # ms
        
        ### external input ###
        if externalInput:
            distToRs = []
            for releaseId in range(self.num_external_input):
                distER = (distance.cdist(self.coordsRelease[[releaseId],:],self.coordsE,
                                         lambda a,b: self.computeDist(a,b))
                                         .reshape(-1,1))
                
                distIR = (distance.cdist(self.coordsRelease[[releaseId],:],self.coordsI,
                                         lambda a,b: self.computeDist(a,b))
                                         .reshape(-1,1))  
                distToRs.append(np.vstack((distER,
                                           100*np.ones(shape=distIR.shape))))
                # self.Idrive = DEFAULT_IDRIVE*np.ones(shape=(self.neuroNum,1))
                self.Idrive[distToRs[releaseId]<self.releaseR] = (1+ex_drive_strength) * self.Idrive.min()
            
        ### poisson noise ###        
        if poisson_noise:
            poissonRate = poisson_rate #s-1
            poissonKickAmp = poisson_amp
            poissonKickDur = 1    
            Ipoisson = 0
            
        # ### temp current logger
        # self.meanItotal = 0

        for t in self.tPoints:     
            if logV: 
                self.vPoints[:,[logCounter]] = memV
                self.iPoints[:,[logCounter]] = Itotal
                logCounter += 1
            
            # determine synI vector (for sub class NeuroNet) 
            # and record spike times
            isFiring = (memV < THRESHOLD_AP)
            if isNet:
                EsynMat,memVMat = np.meshgrid(self.Esyn,memV)
                expTerm = np.zeros(shape = (self.neuroNum,1))
                ithLatestSpike = 1
                deltaTs = t - spikeTimes[neuroIdx,spikeCounts-ithLatestSpike]
                while ((deltaTs<windowIsyn) & (spikeCounts>ithLatestSpike)).any():                            
                    expTerm += ((deltaTs < windowIsyn) & 
                                (spikeCounts>ithLatestSpike)) * np.exp(
                                -deltaTs /self.tauSyn)
                    ithLatestSpike += 1 
                    deltaTs = t-spikeTimes[neuroIdx,spikeCounts-ithLatestSpike]
                
                Isyn =self.adjMat * (memVMat - EsynMat) @ expTerm
                Itotal = self.Idrive - Isyn
                # ### temp current logger
                # self.meanItotal += Itotal
        ### poisson noise ###    
            if poisson_noise:
                if not t%poissonKickDur:
                    Ipoisson = poissonKickAmp * (np.random.rand(self.neuroNum,1)<poissonRate)
                Itotal += Ipoisson 
                
            # RK4 method
            kV = np.tile(memV,4)
            kZ = np.tile(channelZ,4)
            kH = np.tile(channelH,4)
            kN = np.tile(channelN,4)
            for colInd in colIdx:
                mInf = 1 / (1 + np.exp((-kV[:,[colInd]]-30.0)/9.5))
                hInf = 1 / (1 + np.exp((kV[:,[colInd]]+53.0)/7.0))
                nInf = 1 / (1 + np.exp((-kV[:,[colInd]]-30.0)/10))
                zInf = 1 / (1 + np.exp((-kV[:,[colInd]]-39.0)/5.0))
                hTau = 0.37 + 2.78 / (1 + np.exp((kV[:,[colInd]]+40.5)/6))
                nTau = 0.37 + 1.85 / (1 + np.exp((kV[:,[colInd]]+27.0)/15))
                fh = (hInf - kH[:,[colInd]]) / hTau
                fn = (nInf - kN[:,[colInd]]) / nTau
                fz = (zInf - kZ[:,[colInd]]) / 75.0
                fv = (1/C)*(g_Na*(mInf**3) * kH[:,[colInd]] *  
                    (v_Na-kV[:,[colInd]]) + 
                    g_Kdr*(kN[:,[colInd]]**4) * (v_K - kV[:,[colInd]])+ 
                    self.gKs * kZ[:,[colInd]] * (v_K - kV[:,[colInd]])+ 
                    g_L*(v_L-kV[:,[colInd]]) + Itotal)
                kH[:,[colInd]] = self.tStep*fh
                kN[:,[colInd]] = self.tStep*fn
                kZ[:,[colInd]] = self.tStep*fz
                kV[:,[colInd]] = self.tStep*fv
                if colInd == 0 or colInd == 1:
                    kH[:,[colInd+1]] = kH[:,[colInd+1]] + 0.5*kH[:,[colInd]]
                    kN[:,[colInd+1]] = kN[:,[colInd+1]] + 0.5*kN[:,[colInd]]
                    kZ[:,[colInd+1]] = kZ[:,[colInd+1]] + 0.5*kZ[:,[colInd]]
                    kV[:,[colInd+1]] = kV[:,[colInd+1]] + 0.5*kV[:,[colInd]]
                elif colInd == 2:
                    kH[:,[colInd+1]] = kH[:,[colInd+1]] + kH[:,[colInd]]
                    kN[:,[colInd+1]] = kN[:,[colInd+1]] + kN[:,[colInd]]
                    kZ[:,[colInd+1]] = kZ[:,[colInd+1]] + kZ[:,[colInd]]
                    kV[:,[colInd+1]] = kV[:,[colInd+1]] + kV[:,[colInd]]
            memV =     memV +     (kV[:,[0]] + 2 * kV[:,[1]] + 
                                   2 * kV[:,[2]] + kV[:,[3]])/6.0
            channelH = channelH + (kH[:,[0]] + 2 * kH[:,[1]] + 
                                   2 * kH[:,[2]] + kH[:,[3]])/6.0
            channelN = channelN + (kN[:,[0]] + 2 * kN[:,[1]] + 
                                   2 * kN[:,[2]] + kN[:,[3]])/6.0
            channelZ = channelZ + (kZ[:,[0]] + 2 * kZ[:,[1]] + 
                                   2 * kZ[:,[2]] + kZ[:,[3]])/6.0
            # RK4 ends
            isFiring &= (memV > THRESHOLD_AP)   
        
        ### STDP part ###       
            # when STDP turned on, initialize adjMat_max,A+, A-, tau+,tau- etc.
            if STDPon: # if STDP rule is taking place
                if not STDPoff: 
                # if STDP has already been turned off, nothing should be done
                    # STDP rule taking effect here! 
                    if isFiring.any(): 
                    # only change weights when at least one cell is firing
                    # This if statement can not combine with above one 
                    # to make sure keep track of time to turn off STDP
                        # iteration for get all the terms 
                        # within cutoff STDP time window    
                        
                       
                        
                        
                        ithLatestSpike = 1
                        deltaTs = t-spikeTimes[neuroIdx,spikeCounts-ithLatestSpike]
                        # if spikeCounts is zeros then -1 index leads to time at 0
                        # deltaWeights = 0
                        deltaWeightsPlus,deltaWeightsMinus = 0,0
                        
                        ### nearest spike
                        deltaWeightsPlus += (deltaTs < windowSTDP) * np.exp(
                                                            -deltaTs / tauPlus)
                        deltaWeightsMinus += (deltaTs < windowSTDP) * np.exp(
                                                            -deltaTs / tauMinus) * 0.5
                        
                        
                        # STDPAdjMat[idxPostSyn[(isFiring&depressions)[:numPostSyn]],:] -= 
                        
                        STDPAdjMat[idxPostSyn[isFiring[:numPostSyn]],:] += (
                            deltaWeightConst[idxPostSyn[isFiring[:numPostSyn]],:]
                            * deltaWeightsPlus[:numPreSyn].T)
                        
                        STDPAdjMat[:,idxPreSyn[isFiring[:numPreSyn]]] -= (
                            deltaWeightConst[:,idxPreSyn[isFiring[:numPreSyn]]]
                            * deltaWeightsMinus[:numPostSyn]) 

                       
                        # make sure weights in [0,weightmax]
                        STDPAdjMat[STDPAdjMat>STDPAdjMatMax] = STDPAdjMatMax[
                            STDPAdjMat>STDPAdjMatMax]
                        STDPAdjMat[STDPAdjMat<STDPAdjMatMin] = STDPAdjMatMin[
                            STDPAdjMat<STDPAdjMatMin]                        
                        # STDP update done!            
                    if t>self.tSTDP_off: # time to turn off STDP rule
                        STDPoff = True
            elif isSTDP and t>self.tSTDP_on: # turn on STDP at the right time
                STDPon = True
                # initialize important STDP parameters
                numPreSyn = self.numE # considering all excitatory synapses
                numPostSyn = self.numE         
                idxPreSyn = np.arange(numPreSyn).reshape(-1,1)
                idxPostSyn = np.arange(numPostSyn).reshape(-1,1)
                STDPAdjMat = self.adjMat[:numPostSyn,:numPreSyn].copy()
                if silentSynapse:
                    STDPAdjMatEE = STDPAdjMat[:self.numE,:self.numE]
                    tempVec = STDPAdjMatEE[STDPAdjMatEE!=0]
                    synapNum = tempVec.size
                    weightEE = tempVec[0]
                    silentNum =  round(synapNum*silentSynapse)                    
                    # except the diagonal elements
                    tempVec = STDPAdjMatEE[(STDPAdjMatEE+np.eye(self.numE))==0]
                    tempVec[np.random.choice(tempVec.size,
                                             silentNum,
                                             replace=False)] = weightEE
                    STDPAdjMatEE[(STDPAdjMatEE+np.eye(self.numE))==0] = tempVec
                    
                    
                STDPAdjMatMax = STDPAdjMat * (1 + self.STDPlevel)
                STDPAdjMatMin = STDPAdjMat * (1 - self.STDPlevel)                
                deltaWeightConst = STDPAdjMat * self.STDPlevel/20.0       
                STDPAdjMat = self.adjMat[:numPostSyn,:numPreSyn]
                tauSTDP = 10 # ms
                # assymetrical STDP learning rule
                tauPlus = 14 # ms
                tauMinus = 34 #ms
                
                windowSTDP = 100 # ms
                
                
            spikeTimes[neuroIdx[isFiring],spikeCounts[isFiring]] = t
            spikeCounts += isFiring
        # main simulation over

        # compress spikeTimes to a 2D array        
        timingVec = np.concatenate(
            [spikeTimes[i,:spikeCounts[i,0]] for i in neuroIdx.reshape(-1)])
        idVec = np.concatenate(
            [i*np.ones(spikeCounts[i,0]) for i in neuroIdx.reshape(-1)])
        self.spikeTimes = np.stack((timingVec,idVec))
        self.spikeCounts = spikeCounts
        if not isNet: self.states = np.hstack((channelZ,channelH,channelN,memV))
        return self
        # return spikeCounts    

    def detectRhythm(self,tMin=DEFAULT_TEND-4000,tMax=DEFAULT_TEND):
        thresholdTheta = 2
        thresholdGamma = 2
        infThetaBand = 2.5
        supThetaBand = 20
        infGammaBand = 30
        supGammaBand = 200
        freqUpperLimit = 200
        timeWindow = 0.5 * 1000/freqUpperLimit # ms
        
        tempSpikeTimes = self.spikeTimes[:,
                    (self.spikeTimes[0,:]>tMin) & (self.spikeTimes[0,:]<tMax)]
        timePoints = np.arange(tMin,tMax,timeWindow)
        iterationTimes = timePoints.shape[0]
        logicalRaster = np.zeros(shape=(iterationTimes,self.neuroNum))
        
        for ithTimeWindow in range(iterationTimes):
            temp1 = np.stack((np.abs(tempSpikeTimes[0,:]-
                                     timePoints[ithTimeWindow]),
                              tempSpikeTimes[1,:]),axis = 0)
            temp2 = temp1[1,temp1[0,:] <= timeWindow/2]
            logicalRaster[[ithTimeWindow],temp2.astype(int)] = 1
        logicalRaster = logicalRaster.T
        fPoints,PxxDensity = signal.periodogram(logicalRaster,
                                        fs=1000/timeWindow)
        
        # Network Pxx Density and its normalization
        netPxxDensity = PxxDensity.mean(axis=0)
        netPxxDensity = netPxxDensity/netPxxDensity.mean() 
        self.fPoints = fPoints
        self.netPxxDensity = netPxxDensity
        # log peak power and corresponding freq in theta and gamma band
        # peak power
        self.Ptheta = netPxxDensity[
            (fPoints>infThetaBand)&(fPoints<supThetaBand)].max()
        self.Pgamma = netPxxDensity[
            (fPoints>infGammaBand)&(fPoints<supGammaBand)].max()
        # freq part
        if self.Ptheta > thresholdTheta:
            self.thetaFreq = fPoints[
                (fPoints>infThetaBand)&(fPoints<supThetaBand)][
                    netPxxDensity[(
                        fPoints>infThetaBand)&(fPoints<supThetaBand)].argmax()] 
        else:
            self.thetaFreq = np.nan
        
        if self.Pgamma > thresholdGamma:
            self.gammaFreq = fPoints[
                (fPoints>infGammaBand)&(fPoints<supGammaBand)][
                    netPxxDensity[(
                        fPoints>infGammaBand)&(fPoints<supGammaBand)].argmax()] 
        else:
            self.gammaFreq = np.nan
        
        # Use neuronal Pxx Density to map neuronal type (rhythmic)
        self.neuroRhythmType = np.zeros(self.neuroNum) # 0 for null type
        neuroPtheta = PxxDensity[:,(fPoints>infThetaBand)&(fPoints<supThetaBand)].max(axis=1)
        self.neuroRhythmType[neuroPtheta>thresholdTheta*PxxDensity.mean(axis=1)] = 1 # 1 for theta type
        neuroPgamma = PxxDensity[:,(fPoints>infGammaBand)&(fPoints<supGammaBand)].max(axis=1)
        self.neuroRhythmType[neuroPgamma>thresholdGamma*PxxDensity.mean(axis=1)] = 2 # 2 for gamma type
        self.neuroRhythmType[
            (neuroPtheta>thresholdTheta*PxxDensity.mean(axis=1)) & 
            (neuroPgamma>thresholdGamma*PxxDensity.mean(axis=1))] = 3 # 3 for mixed     
        return self
            
    def showRaster(self):
        # preliminary method needs improvement
        plt.figure()
        plt.plot(self.spikeTimes[0,:],self.spikeTimes[1,:],'o',markersize = 2)
        plt.xlim(self.tEnd - 500,self.tEnd)
        plt.ylim(0,500)
        plt.show()
        
        
    def rewireEE(self,rewiringProb=0.2): # assuming const weight
        adjMatEE = self.adjMat[:self.numE,:self.numE]
        tempVec = adjMatEE[adjMatEE!=0]
        synapNum = tempVec.shape[0]
        weightEE = tempVec[0]
        rewiringNum =  round(synapNum*rewiringProb)
        breakId = np.random.choice(synapNum,rewiringNum,replace=False)
        tempVec[breakId] = 0
        adjMatEE[adjMatEE!=0] = tempVec
        # except the diagonal elements
        tempVec = adjMatEE[(adjMatEE+np.eye(self.numE))==0]
        tempVec[np.random.choice(tempVec.shape[0],rewiringNum,replace=False)] = weightEE
        adjMatEE[(adjMatEE+np.eye(self.numE))==0] = tempVec
        return self
    
    def sparsenEE(self,sparsity=0.5): # assuming const weight
        adjMatEE = self.adjMat[:self.numE,:self.numE]
        tempVec = adjMatEE[adjMatEE!=0]
        synapNum = tempVec.shape[0]

        breakId = np.random.choice(synapNum,
                                   round(synapNum*sparsity),
                                   replace=False)
        tempVec[breakId] = 0
        adjMatEE[adjMatEE!=0] = tempVec
        
        return self        
