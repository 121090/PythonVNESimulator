# -*- coding: utf-8 -*-
"""
Created on Sun Aug 04 22:52:57 2019

@author: LTI
"""
import numpy as np
import SimulationTools as ST

class objectiveFunction:
    def __init__(self,optFX, flagForwarding, inputDirectory,weights,representation):
        self.optFX=optFX
        self.flagForwarding=flagForwarding
        if self.flagForwarding==1:
            self.vectorForwarding=np.loadtxt(str(inputDirectory)+'/forwarding.csv',delimiter=',')
        else:
             self.vectorForwarding=[]
        self.weights=weights
        self.representation=representation
        
    def initializeMax(self):
        self.maxG=[0,0]
        self.maxFX=[0,0,0,0,0,0]
        self.minFX=[float('inf'),float('inf'),float('inf'),float('inf'),float('inf'),float('inf')]
    
    def printSolutionsFX(self,solutions):
        for i in range(len(solutions)):
            if self.representation==0:
                cad= "Solution: "+str(solutions[i])
            else:
                cad= "Solution: "+str(solutions[i].factoradic)
            print cad+", cost: "+str(solutions[i].vneCost)+", FXs: "+str(solutions[i].FXs)+", constraints: "+str(solutions[i].constraints)+", penaltyValue: "+str(solutions[i].penaltyValue)                
    
    def maxPositionFX(self,solutions,version):
        maxValue=0
        maxPenaltyValue=0
        maxPosition=0
        maxDegreeDiference=0
        for i in range(len(solutions)):
            solutions[i].getAptitude(self.optFX,self.weights,self.maxFX,self.minFX)
            if solutions[i].aptitude>maxValue:
                maxValue=solutions[i].aptitude
                maxPenaltyValue=solutions[i].penaltyValue
                maxDegreeDiference=solutions[i].FXs[1]
                maxPosition=i
            elif solutions[i].aptitude==maxValue:
                if maxValue==float('inf') and version>3:
                    if solutions[i].penaltyValue>maxPenaltyValue:
                        maxPenaltyValue=solutions[i].penaltyValue
                        maxPosition=i
                else:
                    if solutions[i].FXs[1]>maxDegreeDiference and self.optFX==1:
                        maxDegreeDiference=solutions[i].FXs[1]
                        maxPosition=i
        return maxPosition
    
    def minPositionFX(self,solutions,version,alg=None):
        minvalue=float("inf")
        minPenaltyValue=float("inf")
        minPosition=0
        minDegreeDiference=float("inf")
        for i in range(len(solutions)):
            if alg!=None:
                solutions[i].getAptitude(self.optFX,self.weights,self.maxFX,self.minFX)
            if solutions[i].aptitude<minvalue:
                minvalue=solutions[i].aptitude
                minPenaltyValue=solutions[i].penaltyValue
                minDegreeDiference=solutions[i].FXs[1]
                minPosition=i
            elif solutions[i].aptitude==minvalue:
                if minvalue==float('inf') and version>3:
                    if solutions[i].penaltyValue<minPenaltyValue:
                        minPenaltyValue=solutions[i].penaltyValue
                        minPosition=i
                else:
                    if solutions[i].FXs[1]<minDegreeDiference and self.optFX==1:
                        minDegreeDiference=solutions[i].FXs[1]
                        minPosition=i
        return minPosition
    
    def mapLinks(self,vector,PhysicalNetwork,virtualNetwork):
        paths=[]
        infeasiblePaths=0
        pathsCost=0
        adjacencyList=virtualNetwork.graph.get_adjlist(mode='ALL')
        c=0
        bandwidthNormalization=0
        #bandera=0
        #print "__________Solucion: "+str(vector)+"__________"
        for i in range(virtualNetwork.graph.vcount()):    
            for j in range(len(adjacencyList[i])):
                if adjacencyList[i][j]>i:
                    #print "******************************************************"
                    #print "Enlace virtual: "+str(i)+" , "+str(listaAdyacencia[i][j])
                    bandwitdhRequired=virtualNetwork.graph.es[virtualNetwork.graph.get_eid(adjacencyList[i][j],i)]["bandwidth"]
                    graphCopy=ST.prune(PhysicalNetwork.graph,bandwitdhRequired,self.vectorForwarding)
                    #if copia.ecount()==66 and bandera==0:
                        #print copia
                     #   bandera=2
                    #print copia
                    source=vector[i]
                    dest=vector[adjacencyList[i][j]]
                    #print "Enlace Fisico: "+ str(origen)+" , "+str(destino)
                    #camino=copia.shortest_paths_dijkstra(source=origen, target=destino, weights=None, mode='ALL')
                    #camino=copia.shortest_paths(source=origen, target=destino, weights=None, mode='OUT')
                    path=graphCopy.get_shortest_paths(v=source, to=dest, weights=None, mode='ALL',output="vpath")
                    #print "sdjdhkdhkdjhdkjhdkjhdjkhdjkdhfjkhfkjhfjkfhkjfh"
                    c=c+1
                    #print camino
                    #if virtualNetwork.vcount()==8 and virtualNetwork.ecount()==16:
                    #print "camino: "+str(c)+" tamano: "+str(len(camino[0]))+" nodos: "+str(camino)+ " BW: "+str(anchoBanda)
                    if(len(path[0])==0):
                        infeasiblePaths=infeasiblePaths+1
                        pathCost=0
                        #print "no se encontro camino"
                    else:
                        if len(self.vectorForwarding)==0:
                            pathCost=(bandwitdhRequired*(len(path[0])-1))
                        else:
                            pathCost=(bandwitdhRequired*(len(path[0])-1))+(self.vectorForwarding[int(bandwitdhRequired)-1]*(len(path[0])-2))
                    #print (vectorForwarding[int(anchoBanda)-1]*(len(camino[0])-2))
                    #print path[0]
                    paths.append(ST.Path(path[0],bandwitdhRequired,pathCost))
                    pathsCost=pathsCost+pathCost
                    bandwidthNormalization+=PhysicalNetwork.reserveFreeResourcesPath(paths[len(paths)-1],0,self.vectorForwarding,self.optFX)
        PhysicalNetwork.reserveFreeResourcesPaths(paths,1,self.vectorForwarding)
        #print PhysicalNetwork.es["bandwidth"]
        #print costoCaminos
        return paths,pathsCost,infeasiblePaths,bandwidthNormalization           
    
    def evaluateSolutions(self,solutions,physicalNetwork,virtualNetwork,phy_degree,v_degree):
        infeasibleSolution=[]
        for i in range(len(solutions)):
            if solutions[i].vneCost==0:
                infeasibleNodes=0
                sumCPU=0
                degreeDiference=0
                cpuDiference=0
                cpuNormalization1=0
                cpuNormalization2=0
                cpuNormalization3=0
                #print "Antes caminos"
                paths,linksCost,infeasibleLinks,bandwidthNormalization=self.mapLinks(solutions[i].permutation,physicalNetwork,virtualNetwork)
                solutions[i].pathsList=paths
                #print "enlacesNoFactibles: "+str(enlacesNoFactibles)
                for j in range(len(solutions[i].permutation)):
                    #print poblacion[i].permutacion[j]
                    #print redVirtual.vs[j]["CPU"]
                    #print redFisica.vs[poblacion[i].permutacion[j]]["CPU"]
                    physicalCPU=physicalNetwork.graph.vs[solutions[i].permutation[j]]["CPU"]
                    virtualCPU=virtualNetwork.graph.vs[j]["CPU"]
                    if physicalCPU<virtualCPU:
                        #print "hola"
                        infeasibleNodes=infeasibleNodes+1
                    sumCPU=sumCPU+virtualCPU#virtualNetwork.graph.vs[j]["CPU"]
                    if self.optFX==1:
                        degreeDiference=degreeDiference+(phy_degree[solutions[i].permutation[j]]-v_degree[j])
                    if self.optFX==2:
                        cpuDiference+=physicalCPU-virtualCPU
                    if self.optFX==3:
                        cpuNormalization1+=(physicalCPU+virtualCPU+0.0)
                    if self.optFX==4:
                        cpuNormalization2+=(1.0/(physicalCPU+np.spacing(1)))
                    if self.optFX==5:
                        cpuNormalization3+=virtualCPU*(1.0/(physicalCPU+np.spacing(1)))
                solutions[i].vneCost=sumCPU+linksCost
                solutions[i].FXs[0]=linksCost
                solutions[i].FXs[1]=degreeDiference
                solutions[i].FXs[2]=cpuDiference
                solutions[i].FXs[3]=cpuNormalization1
                solutions[i].FXs[4]=cpuNormalization2+bandwidthNormalization
                solutions[i].FXs[5]=cpuNormalization3+bandwidthNormalization
                if infeasibleNodes>0 or infeasibleLinks>0:
                    infeasibleSolution.append(i)
                    solutions[i].constraints[0]=infeasibleNodes
                    solutions[i].constraints[1]=infeasibleLinks
                    if infeasibleNodes>self.maxG[0]:
                        self.maxG[0]=infeasibleNodes
                    if infeasibleLinks>self.maxG[1]:
                        self.maxG[1]=infeasibleLinks
                else:
                    if solutions[i].FXs[0]>self.maxFX[0]:
                        self.maxFX[0]=solutions[i].FXs[0]
                    if solutions[i].FXs[1]>self.maxFX[1]:
                        self.maxFX[1]=solutions[i].FXs[1]
                    if solutions[i].FXs[0]<self.minFX[0]:
                        self.minFX[0]=solutions[i].FXs[0]
                    if solutions[i].FXs[1]<self.minFX[1]:
                        self.minFX[1]=solutions[i].FXs[1]
        for i in range(len(infeasibleSolution)):
            penaltyValue=0
            if self.maxG[0]!=0:
                penaltyValue=penaltyValue+(float(solutions[infeasibleSolution[i]].constraints[0])/self.maxG[0])
            if self.maxG[1]!=0:
                penaltyValue=penaltyValue+(float(solutions[infeasibleSolution[i]].constraints[1])/self.maxG[1])
            solutions[infeasibleSolution[i]].penaltyValue=penaltyValue
            
    
    def evaluateOneSolution(self,solution,physicalNetwork,virtualNetwork,phy_degree,v_degree):
        infeasibleNodes=0
        sumCPU=0
        degreeDiference=0
        cpuDiference=0
        cpuNormalization1=0
        cpuNormalization2=0
        cpuNormalization3=0
        paths,linksCost,infeasibleLinks,bandwidthNormalization=self.mapLinks(solution.permutation,physicalNetwork,virtualNetwork)
        solution.pathsList=paths
        for j in range(len(solution.permutation)):
            physicalCPU=physicalNetwork.graph.vs[solution.permutation[j]]["CPU"]
            virtualCPU=virtualNetwork.graph.vs[j]["CPU"]
            if physicalCPU<virtualCPU:
                infeasibleNodes=infeasibleNodes+1
            sumCPU=sumCPU+virtualNetwork.graph.vs[j]["CPU"]
            if self.optFX==1:
                degreeDiference=degreeDiference+(phy_degree[solution.permutation[j]]-v_degree[j])
            if self.optFX==2:
                cpuDiference+=physicalCPU-virtualCPU
            if self.optFX==3:
                cpuNormalization1+=(physicalCPU+virtualCPU+0.0)
            if self.optFX==4:
                cpuNormalization2+=(1.0/(physicalCPU+np.spacing(1)))
            if self.optFX==5:
                cpuNormalization3+=virtualCPU*(1.0/(physicalCPU+np.spacing(1)))
        solution.vneCost=sumCPU+linksCost
        solution.FXs[0]=linksCost
        solution.FXs[1]=degreeDiference
        solution.FXs[2]=cpuDiference
        solution.FXs[3]=cpuNormalization1
        solution.FXs[4]=cpuNormalization2+bandwidthNormalization
        solution.FXs[5]=cpuNormalization3+bandwidthNormalization
        if infeasibleNodes>0 or infeasibleLinks>0:
            solution.constraints[0]=infeasibleNodes
            solution.constraints[1]=infeasibleLinks
            if infeasibleNodes>self.maxG[0]:
                self.maxG[0]=infeasibleNodes
            if infeasibleLinks>self.maxG[1]:
                self.maxG[1]=infeasibleLinks
            penaltyValue=0
            if self.maxG[0]!=0:
                penaltyValue=penaltyValue+(float(solution.constraints[0])/self.maxG[0])
            if self.maxG[1]!=0:
                penaltyValue=penaltyValue+(float(solution.constraints[1])/self.maxG[1])
            solution.penaltyValue=penaltyValue
        else:
            if solution.FXs[0]>self.maxFX[0]:
                self.maxFX[0]=solution.FXs[0]
            if solution.FXs[1]>self.maxFX[1]:
                self.maxFX[1]=solution.FXs[1]
            if solution.FXs[0]<self.minFX[0]:
                self.minFX[0]=solution.FXs[0]
            if solution.FXs[1]<self.minFX[1]:
                self.minFX[1]=solution.FXs[1]
        solution.getAptitude(self.optFX,self.weights,self.maxFX,self.minFX)

