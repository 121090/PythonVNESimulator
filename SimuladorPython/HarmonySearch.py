#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 09:46:08 2019

@author: christianaguilarfuster
"""
import numpy as np
import SimulationTools as ST


class HarmonySearch:
    
    def __init__(self,version,HMSize,maxIterations,hcmr,par):
        self.HMSize=HMSize
        self.maxIterations=maxIterations
        self.hcmr=hcmr
        self.par=par
        self.version=version
        self.candidateNodes=None
        
    def setCandidateNodes(self, candidateNodes,order):
        self.candidateNodes=candidateNodes
        self.order=order
        #print "Order: "+str(self.order)
        #for i in range(len(candidateNodes)):
         #   print candidateNodes[i]

    def improviseNewHarmony(self,solutions,amountNodes):
        amountNotes=len(solutions[0].permutation)
        newVector=np.ones(amountNotes,dtype=int)*amountNodes
        for i in range(amountNotes):
            coin=np.random.uniform()
            if coin<=self.hcmr:
                memoryPosition=np.random.randint(0,len(solutions))
                number=solutions[memoryPosition].permutation[i]
                count=0
                while number in newVector:
                    count+=1
                    memoryPosition=np.random.randint(0,len(solutions))
                    number=solutions[memoryPosition].permutation[i]     
                    if count>1000:
                        number=np.random.randint(0,amountNodes) 
                newVector[i]=number
                coin=np.random.uniform()
                if coin<=self.par/float(2):
                    number=newVector[i]-1
                    if number<0:
                        number=amountNodes-1
                    while number in newVector:
                        number=number-1
                        if number<0:
                            number=amountNodes-1
                    newVector[i]=number
                elif coin<=self.par:
                    number=newVector[i]+1
                    #print numero
                    if number>=amountNodes:
                        number=0
                    while number in newVector:
                        number=number+1
                        if number>=amountNodes:
                            number=0
                    #print numero
                    newVector[i]=number
                    
            else:
                number=np.random.randint(0,amountNodes)
                while number in newVector:
                    number=np.random.randint(0,amountNodes)                   
                newVector[i]=number
        #print newVector
        return ST.Solution(newVector.copy())
    
    def improviseNewHarmony2(self,solutions,amountNodes):
        amountNotes=len(solutions[0].permutation)
        newVector=np.ones(amountNotes,dtype=int)*amountNodes
        for i in range(amountNotes):
            vectorPosition=self.order[i]
            coin=np.random.uniform()
            size=len(self.candidateNodes[vectorPosition])
            if coin<=self.hcmr:
                memoryPosition=np.random.randint(0,len(solutions))
                number=solutions[memoryPosition].permutation[vectorPosition]
                count=0
                while number in newVector:
                    count+=1
                    memoryPosition=np.random.randint(0,len(solutions))
                    number=solutions[memoryPosition].permutation[vectorPosition]     
                    if count>1000:
                        number=np.random.randint(0,len(self.candidateNodes[vectorPosition]))
                        number=self.candidateNodes[vectorPosition][number]
                newVector[vectorPosition]=number
                coin=np.random.uniform()
                #print "number: "+str(number)
                if coin<=self.par/float(2) and size!=1:
                    #print "***************************"+str(vectorPosition)
                    #print self.candidateNodes[vectorPosition]
                    indexList=self.candidateNodes[vectorPosition].index(newVector[vectorPosition])
                    #print indexList
                    newPosition=indexList-1
                    if newPosition<0:
                        newPosition=len(self.candidateNodes[vectorPosition])-1
                    number=self.candidateNodes[vectorPosition][newPosition]
                    while number in newVector:
                        newPosition=newPosition-1
                        if newPosition<0:
                            newPosition=len(self.candidateNodes[vectorPosition])-1
                        number=self.candidateNodes[vectorPosition][newPosition]
                    newVector[vectorPosition]=number
                elif coin<=self.par and size!=1:
                    #print "----------------------------"+str(vectorPosition)
                    #print self.candidateNodes[vectorPosition]
                    indexList=self.candidateNodes[vectorPosition].index(newVector[vectorPosition])
                    #print indexList
                    newPosition=indexList+1
                    if newPosition>=len(self.candidateNodes[vectorPosition]):
                        newPosition=0
                    #print "position: "+str(newPosition)+" size: "+str(size)
                    number=self.candidateNodes[vectorPosition][newPosition]
                    while number in newVector:
                        newPosition=newPosition+1
                        if newPosition>=len(self.candidateNodes[vectorPosition]):
                            newPosition=0
                        #print "position: "+str(newPosition)
                        number=self.candidateNodes[vectorPosition][newPosition]
                    newVector[vectorPosition]=number                    
            else:
                if size!=1:
                    position=np.random.randint(0,size)
                    number=self.candidateNodes[vectorPosition][position]
                    while number in newVector:
                        position=np.random.randint(0,size)
                        number=self.candidateNodes[vectorPosition][position]     
                    newVector[vectorPosition]=number
                else:
                    newVector[vectorPosition]=self.candidateNodes[vectorPosition][0]
            #print newVector
        #print newVector
        return ST.Solution(newVector.copy())

    def compareSolutions(self,solution1,solution2,optFX):
        flag=0
        if solution1.aptitude<solution2.aptitude:
            flag=1
        elif solution1.aptitude==solution2.aptitude:
            if solution1.aptitude==float('inf') and self.version>5:
                if solution1.penaltyValue<solution2.penaltyValue:
                    flag=1
            elif solution1.FXs[1]<solution2.FXs[1] and optFX==1:
                flag=1
        return flag
    
    def startAlgorithm(self,physicalNetwork, virtualNetwork, initialSolutions, objectiveFunction,revenue):
        iterationsCount=0
        physicalNodes=physicalNetwork.graph.vcount()
        hm=initialSolutions
        phy_degree=np.array(physicalNetwork.graph.degree())
        v_degree=np.array(virtualNetwork.graph.degree())
        objectiveFunction.initializeMax()
        #objectiveFunction.printSolutionsFX(hm)
        objectiveFunction.evaluateSolutions(hm,physicalNetwork,virtualNetwork,phy_degree,v_degree)
        while iterationsCount<self.maxIterations:
            #print "---------------------- Generacion: "+str(iterationsCount)+"----------------------"
            #objectiveFunction.printSolutionsFX(hm)
            maxPosition=objectiveFunction.maxPositionFX(hm,self.version)
            if self.candidateNodes==None:
                newHarmony=self.improviseNewHarmony(hm,physicalNodes)
            else:
                newHarmony=self.improviseNewHarmony2(hm,physicalNodes)
            #print newHarmony
            objectiveFunction.evaluateOneSolution(newHarmony,physicalNetwork,virtualNetwork,phy_degree,v_degree)
            #objectiveFunction.printSolutionsFX([newHarmony])
            flag=self.compareSolutions(newHarmony,hm[maxPosition],objectiveFunction.optFX)
            if flag==1:
                hm[maxPosition]=newHarmony
            iterationsCount=iterationsCount+1
        minPosition=objectiveFunction.minPositionFX(hm,self.version)
        return hm[minPosition]