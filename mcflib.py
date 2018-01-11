
#!/usr/bin/env python
from __future__ import print_function
import csv
import sys, getopt
import numpy


from ortools.graph import pywrapgraph
from apgl.graph import *
from utils import *




def build_graph(G, min_cost_flow, normalPatterns, tumorPatterns, start , end):
    #source supply for mcf
    min_cost_flow.SetNodeSupply(0, 100)
    
    normalCounts = len(normalPatterns)
    tumorCounts = len(tumorPatterns)
    
    normalSum = 0
    for normal in normalPatterns:
        normalSum += normal.abundance
    
    
    tumorSum = 0
    for tumor in tumorPatterns:
        tumorSum += tumor.abundance
    
    #print("tumorSum = ", tumorSum)
    
    nRes = 100
    for normal in normalPatterns[:-1]:
        #print("normalCounts = ", normalCounts)
        #add edge btw source and normal patterns
        
        G.addEdge(0 , normal.getId(), (int(normal.abundance*100/normalSum), 0))
        nRes -= int(normal.abundance*100/normalSum)
        #normal supplies for mcf
        min_cost_flow.SetNodeSupply(normal.getId(), 0)
        
        #print("Normalabundance = ", int(normal.abundance*100/normalSum))
        #print("normalIDs = ", normal.getId())
        #tumorCounts = 0
        for tumor in tumorPatterns:
            #tumorCounts+=1
            #print("tumorCounts = ",tumorCounts)
            #print("normal mpat = ", normal)
            #print("tumor mpat = ", tumor)
            #add edge btw normal and tumor patterns
            #print("Tumor Abundance = ", int(tumor.abundance*100/tumorSum))
            #print("tumorIDs = ", tumor.getId())
            G.addEdge(normal.getId(), tumor.getId(), (100, int(float(getJaccardDistance(normal, tumor, start, end))*100)))



    G.addEdge(0 , normalPatterns[-1].getId(), (nRes, 0))
    #normal supplies for mcf
    min_cost_flow.SetNodeSupply(normalPatterns[-1].getId(), 0)
    for tumor in tumorPatterns:
        #print("tumorCounts = ", tumorCounts)
        #add edge btw normal and tumor patterns
        G.addEdge(normalPatterns[-1].getId(), tumor.getId(), (100, int(float(getJaccardDistance(normalPatterns[-1], tumor, start, end))*100)))


    tRes = 100
    for tumor in tumorPatterns[:-1]:
        #tumor supply nodes
        min_cost_flow.SetNodeSupply(tumor.getId(), 0)
        #add edge btw tumor patterns and sink
        tRes -= int(tumor.abundance*100/tumorSum)
        G.addEdge( tumor.getId(), normalCounts + tumorCounts + 1, (int(tumor.abundance*100/tumorSum), 0))
        #print("Tumorabundance = ", int(tumor.abundance*100/tumorSum))
        #print("tumorIDs = " , tumor.getId())

    G.addEdge( tumorPatterns[-1].getId(), normalCounts + tumorCounts + 1, (tRes, 0))
    
    
    
    #print("here 1")
    #print(G.getAllEdges())
    if G.getAllEdges() == []:
        print("There is no Pattern in Selected Region")
    else:
        for x in G.getAllEdges():
            edge = G.getEdge(x[0],x[1])
            #print(x[0], " ----", x[1])
            #print( "capacity = ", edge[0], ", cost = ", edge[1])
            min_cost_flow.AddArcWithCapacityAndUnitCost(x[0], x[1], edge[0], edge[1])


    #print("here 2")
    
    #print("normalCounts = ", normalCounts)
    #print("tumorCounts = ", tumorCounts)
    min_cost_flow.SetNodeSupply( normalCounts + tumorCounts + 1, -100)





def run_mcf(normalSelectedPat, tumorSelectedPat, start , end):
    G = DictGraph()
    
    min_cost_flow = pywrapgraph.SimpleMinCostFlow()
    
    count = 0
    for normal in normalSelectedPat:
        count += 1
        normal.setId(count)
    
    for tumor in tumorSelectedPat:
        count += 1
        tumor.setId(count)


    build_graph(G, min_cost_flow, normalSelectedPat, tumorSelectedPat, start , end)
    
    
    if min_cost_flow.Solve() == min_cost_flow.OPTIMAL:
        #print('winStart = ', start, '  winEnd = ', end)
        #print()
        #print('Total cost = ', min_cost_flow.OptimalCost())
        #print()
        #for i in range(min_cost_flow.NumArcs()):
        #    print(min_cost_flow.Tail(i), min_cost_flow.Head(i), min_cost_flow.Flow(i))
        return min_cost_flow.OptimalCost()
    else:
        print('winStart = ', start, '  winEnd = ', end)
        print()
        print('There was an issue with the min cost flow input.')
        return -1







def find_overal_intersect(inputFileDirNormal, inputFileDirTumor, start, end):
    print("compute overal Intersection start = ", start)
    print("compute overal Intersection end = ", end)
    
    overalIntersect = []
    overalIntersect.append(Component(0, start, end))
    
    with open(inputFileDirNormal,'rb') as fileN:
        for Nline in fileN:
            patterns = read_methyl_patterns(str(Nline[0:]).strip())
            print(" patterns length = ", len(patterns))
            selectedPatterns = within_region_patterns(patterns, start, end)
            print("selected pattern len = ", len(selectedPatterns))
            componentList = getComponentList(selectedPatterns)
            print(" comp list length = ", len(componentList))

            overalIntersect = IntersectofTwoComponentLists(overalIntersect, componentList)

    
    with open(inputFileDirTumor,'rb') as fileT:
        for Tline in fileT:
            patterns = read_methyl_patterns(str(Tline[0:]).strip())
            selectedPatterns = within_region_patterns(patterns, start, end)
            componentList = getComponentList(selectedPatterns)
            overalIntersect = IntersectofTwoComponentLists(overalIntersect, componentList)

    return overalIntersect


def mcf_for_two(file1, file2, windowStart, windowEnd):
    windowSize = windowEnd - windowStart

    normalPatterns = read_methyl_patterns(file1)
    tumorPatterns = read_methyl_patterns(file2)
    
    normalSelected = within_region_patterns(normalPatterns, windowStart, windowEnd)
    tumorSelected = within_region_patterns(tumorPatterns, windowStart, windowEnd)
   

    return run_mcf(normalSelected, tumorSelected, windowStart, windowEnd)


def compute_distance_matrix_sim(normal_rep_pat, tumor_rep_pat, windowStart, windowEnd):

    normal_rep_size = len(normal_rep_pat)
    tumor_rep_size = len(tumor_rep_pat)
    distMat = numpy.zeros(shape=(normal_rep_size, tumor_rep_size))

    for  i, normal in zip(range(normal_rep_size), normal_rep_pat):
        for j, tumor in zip(range(tumor_rep_size), tumor_rep_pat):
            distMat[i][j] = run_mcf(normal, tumor, windowStart, windowEnd)

    return distMat


def compute_distance_matrix(file1, file2, start, end):
    
    # Creates a list containing 5 lists, each of 8 items, all set to 0
    with open(file1) as f:
        normal_rep_size =  sum(1 for _ in f)
    with open(file2) as f:
        tumor_rep_size =  sum(1 for _ in f)
    distMat = numpy.zeros(shape=(normal_rep_size, tumor_rep_size))

    with open(file1,'rb') as fileN:
        for i, Nline in zip(range(normal_rep_size), fileN):
            with open(file2,'rb') as fileT:
                for j, Tline in zip(range(tumor_rep_size), fileT):
                    normalFile = str(Nline[0:]).strip()
                    tumorFile = str(Tline[0:]).strip()
                    print("normal file befor runForTwo = ", normalFile)
                    print("tumor file befor runForTwo = ", tumorFile)
                        
                    distMat[i][j] = mcf_for_two(str(Nline[0:]).strip(), str(Tline[0:]).strip(), start, end)
    return distMat
