import argparse #Parser for command line options
import SimulationTools as ST 
import AlgorithmVNE as VNE
import ObjectiveFunction as OF
import numpy as np #Numpy Library
from time import time
import sys
import os,errno
import pickle

#Class that contains the simulation procces
class Simulator:
    #Method to initialize the variables of simulation
    #InputDirectory=Directory tha contains the physical and virtual networks and simulation File.
    #OutputDirectory=Directory where the results of the simulation will be written
    #SimulationFile=Name of the file that contains the order, arrival time, and lifetime of the requests
    #Number of request to Simulate.
    def __init__(self, inputDirectory, outputDirectory,simulationFile,amountRequest):
        self.inputDirectory=inputDirectory
        self.outputDirectory=outputDirectory
        self.simulationFile=simulationFile
        self.amountRequest=amountRequest
        self.priorityQueue=ST.readSimulationFile(inputDirectory,simulationFile,amountRequest)
        self.countRequestAccepted=0
        self.countRequestRejected=0
        self.solutions=[]
        self.totalCost=0
        self.totalRevenue1=0
        self.totalRevenue2=0
        self.optimalSolutions=0
        self.executionTime=0
        self.minCPU=[]
        self.minAverageCPU=[]
        self.minBandwith=[]
        self.minAverageBandwith=[]
        self.minPhysicalRevenue=[]
        self.rejectedRequest=[]
        self.links=[]
        self.penaltyValue=[]
        self.physicalNetwork=ST.Network(inputDirectory,ST.PHYSICAL)
        self.initialNetwork=ST.Network(inputDirectory,ST.PHYSICAL)
        self.requestList=ST.readRequests(inputDirectory,amountRequest)    

    #Method to write the results of the simulation.
    #Name of the file that will contain the results of the simulation.
    #Parameters that the VNE algorithm used during the simulation.
    def writeOutput(self,outputFileName,parameters):
        try:
            os.mkdir(outputDirectory)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        f = open(str(self.outputDirectory)+"/"+str(outputFileName),"a")
        f.write(str(self.countRequestAccepted)+","+str(float(self.countRequestAccepted)/self.amountRequest)+","+str(self.totalRevenue1)+","+str(self.totalRevenue2))
        f.write(","+str(self.totalCost)+","+str(float(self.totalRevenue2)/self.totalCost)+","+str(self.optimalSolutions)+","+str(self.executionTime))
        f.write(","+str(float(self.totalRevenue2)/self.countRequestAccepted)+","+str(float(self.totalCost)/self.countRequestAccepted)+",")
        f.write(str((float(self.totalRevenue2)/self.countRequestAccepted)/(float(self.totalCost)/self.countRequestAccepted))+",")
        f.write(str(parameters)+"\n")
        #f.write(str(estadisticas[4])+","+str(estadisticas[0])+","+str(estadisticas[1])+","+str(estadisticas[2])+",")
        #f.write(str(estadisticas[3])+","+str(estadisticas[9]/amountRequest)+","+str(estadisticas[5]/amountRequest)+",")
        #f.write(str(estadisticas[6]/amountRequest)+","+str(estadisticas[7]/amountRequest)+","+str(estadisticas[8]/amountRequest)+"\n")
        f.close()
        #print "Best "+str(1-float(self.countRequestAccepted)/self.amountRequest)     
        #print "Best "+str(1-((.9*(float(self.countRequestAccepted)/self.amountRequest))+(.1*(float(self.totalRevenue2)/self.totalCost))))
#        print self.countRequestAccepted
#        print self.amountRequest
#        print self.totalRevenue1
#        print self.totalCost
#        print self.totalRevenue2
#        print 2-(float(self.countRequestAccepted)/self.amountRequest+float(self.totalRevenue2)/self.totalCost)
#        print abs((float(self.countRequestAccepted)/self.amountRequest)-(float(self.totalRevenue2)/self.totalCost))
#        print 1-((.9*(float(self.countRequestAccepted)/self.amountRequest))-(.1*(float(self.totalRevenue2)/self.totalCost)))
    
    #Method that simulates the arrival of the requests
    def simulate(self,vneAlgorithm,objFunction):
        c=0
        while self.priorityQueue.empty()!=True:
            c=c+1
            #Get the first event in the Queue
            event=self.priorityQueue.get()
            #print "_--------------------------------------------------------------------------------------"
            #print "Time: "+str(event[0])+" event type: "+str(event[1])
            #"If" to determine if the event is the arrival of a request
            if event[1].typeRequest==ST.ARRIVAL:
                #We calculate the revenue of the request
                revenue=self.requestList[event[1].idRequest].calculateRevenue()
                self.totalRevenue1+=revenue[0]
                #print "try to embedd the request with revnue"+str(revenue)
#                bw_tmp = np.array(self.physicalNetwork.graph.es["bandwidth"])
#                averagebw=np.average(bw_tmp)
#                cpu_tmp=np.array(self.physicalNetwork.graph.vs["CPU"])
#                averageCPU=np.average(cpu_tmp)
#                if averagebw>40 and averagebw<50 and averageCPU<=64:
#                    print str(averagebw)+"-------------- "+str(averageCPU)
#                    self.physicalNetwork.graph.write_graphml("/Users/CAF/Documents/DoctoralInstances/medios/red"+str(c)+".graphml")
#                if averageCPU<=65:
#                    print "cpu"
#                    print averageCPU
#                if averagebw<=50:
#                    print "bw"
#                    print averagebw
                #We use the VNE Algorithm to find a solution.
                #self.physicalNetwork=object with the current state of physical network.
                #self.requestList[event[1].idRequest] we get the current request form the virtual networks stored in the memory
                #objFunction we pass tothe algorithm the object with the objective function
                #revenue we pass the revenue to the algorithm
                intialTime = time()     
                vneSolution=vneAlgorithm.embeddVirtualNetwork(self.physicalNetwork,self.requestList[event[1].idRequest],objFunction,revenue)
                finalTime = time()
                #print vneSolution.vneCost
                self.executionTime=self.executionTime+(finalTime - intialTime)
                self.solutions.append(vneSolution)
                if vneSolution.penaltyValue==0:
                    #print vneSolution
                    self.countRequestAccepted=self.countRequestAccepted+1
                    self.totalCost=self.totalCost+vneSolution.vneCost
                    self.totalRevenue2=self.totalRevenue2+revenue[0]
                    if vneSolution.vneCost==revenue[0]:
                        self.optimalSolutions=self.optimalSolutions+1
                    self.physicalNetwork.reserveFreeResourcesPaths(vneSolution.pathsList,0,objFunction.vectorForwarding)
                    self.physicalNetwork.reserveFreeResourcesNodes(self.requestList[event[1].idRequest].graph,vneSolution.permutation,0)
                    departureEvent=ST.Request(event[1].idRequest,ST.DEPARTURE,0)
                    departureTime=event[1].lifetime+event[0]
                    self.priorityQueue.put((departureTime,departureEvent))
                else:
                    self.rejectedRequest.append(self.requestList[event[1].idRequest].graph.vcount())
                    self.links.append(self.requestList[event[1].idRequest].graph.ecount())
                    self.penaltyValue.append(vneSolution.constraints)
                    self.countRequestRejected=self.countRequestRejected+1
                pyshicalRevenue=self.physicalNetwork.calculateRevenue()
                self.minPhysicalRevenue.append(float(pyshicalRevenue[0]))
                self.minCPU.append(pyshicalRevenue[1])
                self.minAverageCPU.append(pyshicalRevenue[2])
                self.minBandwith.append(pyshicalRevenue[3])
                self.minAverageBandwith.append(pyshicalRevenue[4])
            #We need to add here an elif if exist an event of elasticity.
            #If the event is a departure free the resources.
            else:
            #print "Evento Salida :" +str(evento[1])
                solution=self.solutions[event[1].idRequest]
                self.physicalNetwork.reserveFreeResourcesPaths(solution.pathsList,1,objFunction.vectorForwarding)
                self.physicalNetwork.reserveFreeResourcesNodes(self.requestList[event[1].idRequest].graph,solution.permutation,1)
            #print "Physical network resources: "+str(self.physicalNetwork.calculateRevenue())
#        print "Amount of processed requests: "+str(self.amountRequest)
#        print "Amount of accepted requests: "+str(self.countRequestAccepted)
#        print "Amount of rejected requests: "+str(self.countRequestRejected)
#        print "Acceptance rate: "+str(float(self.countRequestAccepted)/self.amountRequest)
#        print "Total revenue of requests: "+str(self.totalRevenue1)
#        print "Total revenue of accepted requests: "+str(self.totalRevenue2)
#        print "Total cost of accepted requests: "+str(self.totalCost)
#        print "Revenue/Cost rate of accepted requests: "+str(float(self.totalRevenue2)/self.totalCost)
#        print "Amount of accepted requests with optimal cost: "+str(self.optimalSolutions)
#    print 'El tiempo de ejecucion fue:',tiempo_ejecucion #En segundos
#    print "Tasa de revenue/solicitudes: "+str(float(totalRevenue2)/countRequestAccepted) 
#    print "Tasa de cost/solicitudes: "+str(float(totalCost)/countRequestAccepted) 
#    print "Tasa de (revenue/solicitudes)/(cost/solicitudes): "+str((float(totalRevenue2)/countRequestAccepted)/(float(totalCost)/countRequestAccepted)) 
    #estadisticas=ST.calculateFreeResources(minCPU,minAverageCPU,minBandwith,minAverageBandwith,minPhysicalRevenue)
        
    def __str__(self):
        return "Number of request loaded: "+str(self.priorityQueue.qsize())
    

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="VNE Simulator.")
    parser.add_argument('-i', '--inputDirectory', help="Directory that contains the graphs.", required=True)
    parser.add_argument('-o', '--outputDirectory', help="Directory to save the results.", required=True)
    parser.add_argument('-r', '--Number_of_Request', help="Amount of request to simulate.",type=int, required=True)
    parser.add_argument('-sf', '--Simulation_file', help="File that contains the arrival order.",type=int, required=True)
    parser.add_argument('-alg', '--VNE_algorithm', help="Algorithms: 0.-Genetic, 1.- Harmony Search.",type=int, choices=range(0,2), required=True)
    parser.add_argument('-v', '--Version_algorithm', help="Algorithm version: 0.-basic, 1.-L2S2, 2.-Memory, 3.-Cluster, 4.-Candidate nodes, 5.-DFS, 6.-Penalty function, 7.-B+L2S2+PF, 8.-B+Mem+PF, 9.-B+Clus+PF, 10.-CD+PF, 11.-DFS+PF.", choices=range(0,12),type=int, required=True)
    parser.add_argument('-p', '--paramaters', help="Parameters for the algorithm.", required=True)
    parser.add_argument('-sd', '--Seed', help="Seed for the random number genertion.",type=int, default=0, required=False)
    parser.add_argument('-forw', '--forwarding_version', help="Forwarding is 0.-Activated, 1.-Disable.", choices=range(0,2),type=int, default=0, required=False)
    parser.add_argument('-fx', '--FX_version', help="Objective function 0.-Normal, 1.-AlternDegree, 2.-AlternCPU.", choices=range(0,6), type=int, default=0, required=False)
    parser.add_argument('-we', '--weights', help="Weighst for the altern obejective function", default=0, required=False)
    parser.add_argument('-rep', '--representation', help="Representation is 0.-Permutation, 1.-Factoradic.", choices=range(0,2), default=0, required=False)
    
    args = parser.parse_args()
    #Parameters required for Simulation
    inputDirectory=args.inputDirectory
    outputDirectory=args.outputDirectory
    amountRequest=args.Number_of_Request
    simulationFile=args.Simulation_file

    #parameters Required for VNE algorithms
    vneAlgorithm=args.VNE_algorithm
    algVersion=args.Version_algorithm
    parameters=args.paramaters
    params=parameters
    parameters=parameters.split()
    if len(parameters)==4:
        amountSolutions=int(parameters[0])
        maxIterations=int(parameters[1])
        param1=float(parameters[2])
        param2=float(parameters[3])
    else:
        print argparse.ArgumentTypeError("Parameters are incorrect.")
        sys.exit()
    simulationSeed=args.Seed
    #Parameters for objective function
    flagRep=args.representation
    flagForw=args.forwarding_version
    flagFX=args.FX_version
    weights=args.weights
    if weights!=0:
        #if flagFX==1:
        weights=weights.split()
        flag=False
        if len(weights)==2:
            weights[0]=float(weights[0])
            weights[1]=float(weights[1])
            if weights[0]+weights[1]==1 or weights[0]+weights[1]==2:
                flag=True
        if flag==False:
            print argparse.ArgumentTypeError("Weights are incorrect.")
            sys.exit()
#        else:
#            print argparse.ArgumentTypeError("Flag FX is not activated.")
#            sys.exit()
        
    #print np.random.get_state()
    #np.random.seed(4000)
#    random.seed(simulationSeed)
    np.random.RandomState(seed=simulationSeed)
    simulator=Simulator(inputDirectory,outputDirectory,simulationFile,amountRequest)
    objectiveFunction=OF.objectiveFunction(flagFX,flagForw,inputDirectory,weights,flagRep)
    algorithmVNE=VNE.VNEAlgorithm(vneAlgorithm, algVersion,amountSolutions, maxIterations, param1, param2)
    simulator.simulate(algorithmVNE,objectiveFunction)
    outputFileName="Output-Alg"+str(vneAlgorithm)+"-"+str(algVersion)+"-R"+str(amountRequest)
    outputFileName+="-SF"+str(simulationFile)+"-Forw"+str(flagForw)+"-FX"+str(flagFX)+".csv"
    simulator.writeOutput(outputFileName,params)
    pickle_file = file('globalMemory.pickle', 'w')
    pickle_file2 = file('maxMetrics.pickle', 'w')
    pickle_file3 = file('minMetrics.pickle', 'w')
    pickle.dump(algorithmVNE.globalMemory, pickle_file)
    pickle.dump(algorithmVNE.maxMetrics, pickle_file)
    pickle.dump(algorithmVNE.minMetrics, pickle_file)
   # -i /Users/christianaguilarfuster/Documents/VNE/SolicitudesWaxman2/50 -o /Users/christianaguilarfuster/Documents/VNE/PruebaCluster -r 1000 -sf 500 -a 1 -sd 10 -v 0 -p "6 20 .88 .1" 
   
   #-forw 0 -optFX 0    
   
   #-i /Users/christianaguilarfuster/Documents/VNE/pruebac++ -o /Users/christianaguilarfuster/Documents/VNE/PruebaCluster -r 1 -sf 500 -a 1 -sd 10 -v 0 -p "6 300 .88 .1" -fx 1
    