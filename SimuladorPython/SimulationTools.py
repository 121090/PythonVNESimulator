import Queue as queue
import numpy as np
import igraph as ig

ARRIVAL=0
DEPARTURE=1
PHYSICAL=0
VIRTUAL=1

class Path:
    def __init__(self, path,bandwidth,cost):
        self.path=path
        self.bandwidth=bandwidth
        self.cost=cost
            
    def __str__(self):
        return "["+str(self.path)+", "+str(self.bandwidth)+", "+str(self.cost)+"]"

class Solution:
    def __init__(self, permutation, factoradic=None,vneCost=None,FXs=None,constraints=None,penaltyValue=None,pathList=None):
        self.factoradic=factoradic
        self.permutation=permutation            
        if vneCost==None:
            self.vneCost=0
            self.FXs=np.zeros(6) 
            self.aptidud=0
            self.constraints=np.zeros(2) 
            self.pathsList=None
            self.penaltyValue=0
        else:
            self.vneCost=vneCost
            self.FXs=FXs
            self.constraints=constraints
            self.pathsList=pathList
            self.penaltyValue=penaltyValue
            
    def getAptitude(self,positionFX,weights,maxFXs,minFXs):
        if self.penaltyValue>0:
            self.aptitude=float('inf')
        else:
            if weights==0:
                self.aptitude=self.vneCost
            elif weights[0]!=1:
                Fx0=float(self.FXs[0]-minFXs[0])/(maxFXs[0]-minFXs[0]+np.spacing(1))
                Fx1=float(self.FXs[positionFX]-minFXs[positionFX])/(maxFXs[positionFX]-minFXs[positionFX]+np.spacing(1))
                self.aptitude=(weights[0]*Fx0)+(weights[1]*Fx1)
            else:
                if positionFX<3:
                    self.aptitude=(self.FXs[0])+(self.FXs[positionFX])
                else:
                    self.aptitude=self.FXs[positionFX]
            
    def __str__(self):
        return str(self.permutation)

class Request:
    def __init__(self, idRequest, typeRequest,lifetime):
        self.idRequest = idRequest
        self.typeRequest = typeRequest
        self.lifetime=lifetime
    
    def __str__(self):
        if self.typeRequest==0:
            return "ARRIVAL, "+str(self.idRequest)+" , "+str(self.lifetime)
        else:
            return "DEPARTURE, "+str(self.idRequest)+" , "+str(self.lifetime)

class Network:
    def __init__(self, inputDirectory, networkType,networkId=None):
        self.networkType=networkType
        if networkType==PHYSICAL:
            self.graph=ig.Graph.Read_GraphML(inputDirectory+"/copiaRed0.graphml")
        else:
            self.networkId=networkId
            self.graph=ig.Graph.Read_GraphML(inputDirectory+"/copiaRed"+str(networkId)+".graphml")
            
    def reserveFreeResourcesPath(self,path,option,vectorForwarding,optFX=None):
        size=len(path.path)-1
        bandwidthNormalization=0
        for i in range(size):
            source=path.path[i]
            dest=path.path[i+1]
            current=self.graph.es[self.graph.get_eid(source,dest)]["bandwidth"]
            if optFX>3:
                bandwidthNormalization+=1.0/current
            if option==0:
                if len(vectorForwarding)!=0 and i!=(size-1):
                    #print vectorForwarding[int(camino.bandwidth)-1]
                    self.graph.vs[dest]["CPU"]=self.graph.vs[dest]["CPU"]-vectorForwarding[int(path.bandwidth)-1]
                self.graph.es[self.graph.get_eid(source,dest)]["bandwidth"]=current-path.bandwidth
            else:
                if len(vectorForwarding)!=0 and i!=(size-1):
                    #print vectorForwarding[int(camino.bandwidth)-1]
                    self.graph.vs[dest]["CPU"]=self.graph.vs[dest]["CPU"]+vectorForwarding[int(path.bandwidth)-1]
                self.graph.es[self.graph.get_eid(source,dest)]["bandwidth"]=current+path.bandwidth
        if optFX==5:
            bandwidthNormalization=(path.bandwidth)*bandwidthNormalization
        return bandwidthNormalization 
            #print self.graph.es[self.graph.get_eid(source,dest)]["bandwidth"]
              
    def reserveFreeResourcesPaths(self,paths,option,vectorForwarding):
        for i in range(len(paths)):
            self.reserveFreeResourcesPath(paths[i],option,vectorForwarding)
    
    def reserveFreeResourcesNodes(self, virtualNetwork, asignatedNodes,option):
        for i in range(len(asignatedNodes)):
            #print i
            #print redVirtual.vs[i]["CPU"]
            #print redFisica.vs[asignacionNodos[i]]["CPU"]
            if option==0:
                self.graph.vs[asignatedNodes[i]]["CPU"]=self.graph.vs[asignatedNodes[i]]["CPU"]-virtualNetwork.vs[i]["CPU"]
            else:
                self.graph.vs[asignatedNodes[i]]["CPU"]=self.graph.vs[asignatedNodes[i]]["CPU"]+virtualNetwork.vs[i]["CPU"]
            
    def calculateRevenue(self):
        vector=self.graph.vs["CPU"]
        sumCPU=0
        minCPU=200
        sumBandwith=0
        minBandwith=200
        for i in vector:
            sumCPU=sumCPU+i
            if i<minCPU:
                minCPU=i
        averageCPU=sumCPU/len(vector)
        vector=self.graph.es["bandwidth"]
        for i in vector:
            sumBandwith=sumBandwith+i
            if i<minBandwith:
                minBandwith=i
        averageBandwith=sumBandwith/len(vector)
        return sumCPU+sumBandwith,minCPU,averageCPU,minBandwith,averageBandwith
    
    def __str__(self):
        if self.networkType==PHYSICAL:
            cad= "Physical Network:\n"
        else:
            cad= "Virtual Network :"+str(self.networkId)+"\n"
        return cad+str(self.graph)

def prune(physicalNetwork,bandwidth,vectorForwarding):
        selectedLinks=physicalNetwork.es.select(bandwidth_ge=bandwidth)
        graphCopy = ig.Graph()
        graphCopy.add_vertices(physicalNetwork.vcount())
        if len(vectorForwarding)==0:
            graphCopy.add_edges([(e.source, e.target) for e in selectedLinks])
        else:
            threshold=vectorForwarding[int(bandwidth)-1]
            for e in selectedLinks:  
                if physicalNetwork.vs[e.source]["CPU"]>=threshold and physicalNetwork.vs[e.target]["CPU"]>=threshold:
                    graphCopy.add_edge(e.source, e.target)
        return graphCopy     

def readSimulationFile(inputDirectory,simulationFile,amountRequests):
    priorityQueue = queue.PriorityQueue()
    file = open(str(inputDirectory)+"/simulacion"+str(simulationFile)+".txt", 'r')
    #print('>>> Reading the simulation file')
    for i in range(amountRequests):
        line=file.readline()
        line=line.split(",")
        idRequest=int(line[0])
        arrivalTime=int(line[1])
        lifetime=int(line[2])
        event=Request(idRequest,ARRIVAL,lifetime)
        priorityQueue.put((arrivalTime,event))
    file.close()
    return priorityQueue

def readRequests(inputDirectory,amountRequests):
    requestlist=[]
    for i in range(1,amountRequests+1):
         request=Network(inputDirectory,VIRTUAL,i)
         requestlist.append(request)
    return requestlist