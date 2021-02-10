#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 12:19:18 2019

@author: christianaguilarfuster
"""
import HarmonySearch as HS
import Genetic as GA
import InitializationFunctions as IF
import numpy as np
from scipy import spatial

class storedSolution:
    def __init__(self, metrics):
        self.metrics=metrics
        self.permutation=None
    
    def normalizeMetrics(self,maxs, mins):
        self.normalizedMetrics=[]
        for i in range(len(maxs)):
            self.normalizedMetrics.append(float(self.metrics[i]-mins[i])/float((maxs[i]-mins[i])+np.spacing(1)))
    
    def __str__(self):
        return str(self.metrics)

class VNEAlgorithm:
    def __init__(self,algorithm, version,amountSolutions,maxIterations,param1,param2):
        self.version=version
        self.amountSolutions=amountSolutions
        if algorithm==0:
            #print "VNE Algorithm: Genetic"
            self.algorithm=GA.GeneticAlgorithm(self.version,self.amountSolutions,maxIterations,param1,param2)
        else:
            #print "VNE Algorithm: Harmony Search"
            self.algorithm=HS.HarmonySearch(self.version,self.amountSolutions,maxIterations,param1,param2)
        if self.version==2 or self.version==8 or self.version>=11:
            self.globalMemory = dict()
            self.maxMetrics=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
            self.minMetrics=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    
    def globalMemoryOperations(self,graphMetrics,physicalNetwork,virtualNetwork):
        numberNodes=virtualNetwork.graph.vcount()
        if numberNodes in self.globalMemory:
            #_____________________Actulizando maximos y minimos_______________________
            #print "__________cantidad de metricas: "+str(len(metricas))+"_--------------"
            for position in range(len(graphMetrics.metrics)):
                if graphMetrics.metrics[position]>self.maxMetrics[numberNodes][position]:
                    self.maxMetrics[numberNodes][position]=graphMetrics.metrics[position]
                if graphMetrics.metrics[position]<self.minMetrics[numberNodes][position]:
                    self.minMetrics[numberNodes][position]=graphMetrics.metrics[position]
            graphMetrics.normalizeMetrics(self.maxMetrics[numberNodes], self.minMetrics[numberNodes])
            elements=self.globalMemory.get(numberNodes)
            minDist=1000
            posMinDist=0
            for i in range(len(elements)):
                elements[i].normalizeMetrics(self.maxMetrics[numberNodes], self.minMetrics[numberNodes])
                cosineDist=spatial.distance.cosine(graphMetrics.normalizedMetrics,elements[i].normalizedMetrics)
                if cosineDist<minDist:
                    minDist=cosineDist
                    posMinDist=i
            return elements[posMinDist].permutation.copy()
        else:
            return []

    def embeddVirtualNetwork(self,physicalNetwork, virtualNetwork,objectiveFunction,revenue):
        v_nodes=virtualNetwork.graph.vcount()
        phy_nodes=physicalNetwork.graph.vcount()
        if self.version==1 or self.version==7:
            #print "versions with L2S2"
            initialSolutions=IF.L2S2intialization(self.amountSolutions,v_nodes,phy_nodes,physicalNetwork.graph,virtualNetwork.graph)
        elif self.version==2 or self.version==8:
            #print "versions with memory initialization"
            metrics=virtualNetwork.calculateMetricsGraph(1)
            graphMetrics=storedSolution(metrics)
            metricSolution=self.globalMemoryOperations(graphMetrics,physicalNetwork,virtualNetwork)
            #print metricSolution
            initialSolutions=IF.memoryInitialization(self.amountSolutions,v_nodes, phy_nodes,metricSolution,0)
        elif self.version==3 or self.version==9:
            #print "versions with comunity initialization"
            weight_tmp = np.array(physicalNetwork.graph.es["bandwidth"]) + 1 
            physicalNetwork.graph.es["bandwidth"]=list(weight_tmp.copy())
            comunities=physicalNetwork.graph.community_edge_betweenness(directed=False, weights="bandwidth").as_clustering()
            weight_tmp = np.array(physicalNetwork.graph.es["bandwidth"]) - 1 
            physicalNetwork.graph.es["bandwidth"]=list(weight_tmp.copy())
            comunitiesDict = dict()
            for i in range(len(comunities)):
                #print cluster[i]
                comunitySize=len(comunities[i])
                if comunitySize in comunitiesDict:
                    comunitiesDict.get(comunitySize).append(comunities[i])
                else:
                    comunitiesDict[comunitySize]=[comunities[i]]
            initialSolutions=IF.comunitiesIntialization(self.amountSolutions,v_nodes,phy_nodes,comunitiesDict)
        elif self.version==4 or self.version==10:
            #print "versions with candidate nodes"
            initialSolutions,candidateNodes,order=IF.selectCandidateNodes(self.amountSolutions,physicalNetwork,virtualNetwork)
            if len(initialSolutions)==1:
                return initialSolutions[0]
            self.algorithm.setCandidateNodes(candidateNodes,order)
        elif self.version==5 or self.version==11:
            print "versions with DFS"
        else: 
            #print "versions with basic initialization"
            initialSolutions=IF.initializeSolutions(self.amountSolutions,v_nodes,phy_nodes)
        
        finalSolution=self.algorithm.startAlgorithm(physicalNetwork, virtualNetwork,initialSolutions,objectiveFunction,revenue)
        
        if finalSolution.penaltyValue==0 and (self.version==2 or self.version==8):
            graphMetrics.permutation=finalSolution.permutation.copy()
            if v_nodes in self.globalMemory:
                self.globalMemory.get(v_nodes).append(graphMetrics)
            else:
                #print "---------------------------Maximos y Minimos---------------------------"
                for i in range(len(graphMetrics.metrics)):
                    self.maxMetrics[v_nodes].append(graphMetrics.metrics[i])
                    self.minMetrics[v_nodes].append(graphMetrics.metrics[i])
                self.globalMemory[v_nodes]=[graphMetrics]
        return finalSolution