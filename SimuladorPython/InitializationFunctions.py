#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 10:24:18 2019

@author: christianaguilarfuster
"""
import numpy as np
import SimulationTools as ST

def initializeSolutions(amountSolutions,vectorSize,amountNodes):
    solutions=[]
    for i in range(amountSolutions):
        vector=np.arange(amountNodes)
        np.random.shuffle(vector)
        vector=vector[0:vectorSize]
        solutions.append(ST.Solution(vector.copy()))
    return solutions
        

def L2S2intialization(amountSolutions,vectorSize,amountNodes,physicalNetwork,virtualNetwork):
    solutions=[]
    for j in range(amountSolutions):
        vector=[]
        for i in range(vectorSize):
            flag2=True
            selectedNodes=physicalNetwork.vs.select(CPU_ge=virtualNetwork.vs[i]['CPU'])
            candidateNodesList=[]
            for nodo in selectedNodes:
                candidateNodesList.append(nodo.index)
            #print len(candidateNodesList)
            if len(candidateNodesList)!=0:
                flag2=False
                cpu=np.asarray(physicalNetwork.vs[candidateNodesList]['CPU'])
                sumbandwidth=np.asarray(physicalNetwork.strength(candidateNodesList,weights='bandwidth'))
                physicalNRU=cpu*sumbandwidth
                sumNRU=np.sum(physicalNRU)
                probabilitiesNRU=physicalNRU/float(sumNRU)
                flag=True
                count=0
                while flag==True:
                    count+=1
                    cumulativeSum=0
                    coin=np.random.uniform()
                    for k in range(len(probabilitiesNRU)):
                        cumulativeSum=cumulativeSum+probabilitiesNRU[k]
                        if cumulativeSum>=coin and candidateNodesList[k] not in vector:
                            vector.append(candidateNodesList[k])
                            flag=False
                            break
                    if count==1000:
                        flag=False
                        flag2=True
                        
            if flag2==True:
                number=np.random.randint(0,amountNodes)
                while number in vector:
                    number=np.random.randint(0,amountNodes)                   
                vector.append(number)
        vector=np.asarray(vector)
        solutions.append(ST.Solution(vector.copy()))
    return solutions