
#!/usr/bin/env python
from __future__ import print_function
import csv
import sys, getopt
import numpy as np
import pandas as pd
import random



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


def patterns_2_is_DMR_arrays(selectedPatterns, componentList, rep_size, is_DMR, DMR_change_param, replicate_change_param, comp):
    
    for comp in componentList:
        normalSelected = within_region_patterns(selectedPatterns, comp.start, comp.end)

        #rand = random.randint(0,101)
        rand = random.random()

        if rand < 0.50:
            is_DMR = 1
        else:
            is_DMR = 0

        normal_rep_pats, tumor_rep_pats, case = simulate(normalSelected, rep_size, is_DMR, thr, DMR_change_param, replicate_change_param)

        distMat = compute_distance_matrix_sim(normal_rep_pats, tumor_rep_pats, comp.start, comp.end)

        dist_mean = distMat.mean()
         
        dist_mean_array_all.append(dist_mean)
        is_DMR_array_all.append(is_DMR)
        case_array_all.append(case)

        if case == 'simple':
            dist_mean_array_simple.append(dist_mean)
            is_DMR_array_simple.append(is_DMR)
        else: 
            if case == 'moderate':
                dist_mean_array_moderate.append(dist_mean)
                is_DMR_array_moderate.append(is_DMR)
            else:
                dist_mean_array_hard.append(dist_mean)
                is_DMR_array_hard.append(is_DMR)

    return is_DMR_array_all, is_DMR_array_simple, is_DMR_array_moderate, is_DMR_array_hard


def compute_distance_matrix_sim(normal_rep_pat, tumor_rep_pat, windowStart, windowEnd):

    normal_rep_size = len(normal_rep_pat)
    tumor_rep_size = len(tumor_rep_pat)
    distMat = np.zeros(shape=(normal_rep_size, tumor_rep_size))

    for  i, normal in zip(range(normal_rep_size), normal_rep_pat):
        for j, tumor in zip(range(tumor_rep_size), tumor_rep_pat):
            distMat[i][j] = run_mcf(normal, tumor, windowStart, windowEnd)

    return distMat

def extended_dist_mat_2_distributions(extended_dist_mat, normal_size, tumor_size):
    normal_tumor_array = []
    normal_normal_array = []
    tumor_tumor_array = []

    normal_id_array = []
    tumor_id_array = []
    for i in range(normal_size):
        normal_id_array.append( 'n' + str(i))

    for i in range(tumor_size):
        tumor_id_array.append('t' +str(i))


    for i in range(normal_size):
        for j in range(i+1):
            normal_normal_array.append(extended_dist_mat[normal_id_array[i]][normal_id_array[j]] + random.randint(0,5))

    for i in range(tumor_size):
        for j in range(i+1):
            tumor_tumor_array.append(extended_dist_mat[tumor_id_array[i]][tumor_id_array[j]] + random.randint(0,5))

    for i in range(normal_size):
        for j in range(tumor_size):
            normal_tumor_array.append(extended_dist_mat[normal_id_array[i]][tumor_id_array[j]] + random.randint(0,5))

    return normal_normal_array, tumor_tumor_array, normal_tumor_array


def compute_extended_distance_matrix_sim(normal_rep_pat, tumor_rep_pat, windowStart, windowEnd):

    normal_rep_size = len(normal_rep_pat)
    tumor_rep_size = len(tumor_rep_pat)

    #print("normal_rep_size = " +str(normal_rep_size))
    #print("tumor_rep_size = " +str(tumor_rep_size))

    distMat = np.zeros(shape=(normal_rep_size + tumor_rep_size, normal_rep_size + tumor_rep_size))
    normal_id_array = []
    normal_id_array_b = []
    tumor_id_array = []
    tumor_id_array_b = []

    for i in range(normal_rep_size):
        normal_id_array.append( 'n' + str(i))
        normal_id_array_b.append(0)

    for i in range(tumor_rep_size):
        tumor_id_array.append('t' +str(i))
        tumor_id_array_b.append(1
            )
    id_array =  normal_id_array + tumor_id_array
    id_array_b = normal_id_array_b + tumor_id_array_b


    df_distMat = pd.DataFrame(np.zeros(shape=(normal_rep_size + tumor_rep_size, normal_rep_size + tumor_rep_size)), index = id_array, columns = id_array)


    for  i, normal in zip(range(normal_rep_size), normal_rep_pat):
        for j, tumor in zip(range(tumor_rep_size), tumor_rep_pat):
            df_distMat[normal_id_array[i]][tumor_id_array[j]] = run_mcf(normal, tumor, windowStart, windowEnd)
            df_distMat[tumor_id_array[j]][normal_id_array[i]] = run_mcf(normal, tumor, windowStart, windowEnd)

    for  i, normal_1 in zip(range(normal_rep_size), normal_rep_pat):
        for j, normal_2 in zip(range(normal_rep_size), normal_rep_pat):
            df_distMat[normal_id_array[i]][normal_id_array[j]] = run_mcf(normal_1, normal_2, windowStart, windowEnd)
     
    for  i, tumor_1 in zip(range(tumor_rep_size), tumor_rep_pat):
        for j, tumor_2 in zip(range(tumor_rep_size), tumor_rep_pat):
            df_distMat[tumor_id_array[i]][tumor_id_array[j]] = run_mcf(tumor_1, tumor_2, windowStart, windowEnd)
                    

    return df_distMat, id_array_b

def write_extended_dis_mat(normal_rep_pats, tumor_rep_pats, distMat, comp, is_DMR, case,  output_file):
     with open(output_file, 'a+') as outfile:
        outfile.write('component ' + str(comp.cid)+ ' = (' + str(comp.start)+ ', ' + str(comp.end) + ')\n')
        if is_DMR == 1:
            outfile.write('DMR component\n')
            outfile.write(case + '\n')
        else:
            outfile.write('none-DMR component\n')
            outfile.write(str(case) + '\n')

        np.savetxt(outfile, distMat, fmt='%1.0f')
        for norm_rep in normal_rep_pats:
            for pat in norm_rep:
                outfile.write(str(pat)+'\n')
            outfile.write('-----------------\n')


        outfile.write('-----------------\n')


        for tumor_rep in tumor_rep_pats:
            for pat in tumor_rep:
                outfile.write(str(pat)+'\n')
            outfile.write('-----------------\n')

                                        
        outfile.write('-----------------\n')
        outfile.write('-----------------\n')

def write_single_dist_mat(df_distMat, output_file):
    df_distMat.to_csv(output_file, sep='\t')

def write_health_status(id_array_b, output_file):
    np.savetxt(output_file, id_array_b)

def write_DMR_status(is_DMR, case, output_file):
    
     f = open(output_file,'w') 
     f.write(str(is_DMR) + '\t' + str(case))
     f.close()

def write_auc_pvalue_MiRKAT(auc_pvalue, type_I_error, type_II_error, sensitivity, specificity, output_file):
    f = open(output_file,'a+') 
    f.write(str(auc_pvalue) + '\t' + str(type_I_error)+ '\t' +  str(type_II_error) + '\t')
    f.write(str(sensitivity) + '\t' + str(specificity)+ '\n' )
    
    f.close()

def compute_extended_distance_matrix_from_file(file1, file2, start, end):

    
    with open(file1) as f:
        normal_rep_size =  sum(1 for _ in f)
    with open(file2) as f:
        tumor_rep_size =  sum(1 for _ in f)

    normal_id_array = []
    normal_id_array_b = []
    tumor_id_array = []
    tumor_id_array_b = []


    for i in range(normal_rep_size):
        normal_id_array.append( 'n' + str(i))
        normal_id_array_b.append(0)

    for i in range(tumor_rep_size):
        tumor_id_array.append('t' +str(i))
        tumor_id_array_b.append(1)

    id_array =  normal_id_array + tumor_id_array
    id_array_b = normal_id_array_b + tumor_id_array_b
    


    

    extended_dist_mat_size = normal_rep_size + tumor_rep_size
    df_distMat = pd.DataFrame(np.zeros(shape=(normal_rep_size + tumor_rep_size, normal_rep_size + tumor_rep_size)), index = id_array, columns = id_array)

    with open(file1,'rb') as fileN:
        for  i, Nline in zip(range(normal_rep_size), fileN):
             with open(file2,'rb') as fileT:
                for j, Tline in zip(range(tumor_rep_size), fileT):
                    normalFile = str(Nline[0:]).strip()
                    tumorFile = str(Tline[0:]).strip()
                    mcf_val = mcf_for_two(str(Nline[0:]).strip(), str(Tline[0:]).strip(), start, end) 
                    df_distMat[normal_id_array[i]][tumor_id_array[j]] = mcf_val
                    df_distMat[tumor_id_array[j]][normal_id_array[i]] = mcf_val

    with open(file1,'rb') as fileN1:
        for  i, N1line in zip(range(normal_rep_size), fileN1):
            with open(file1,'rb') as fileN2:
                for j, N2line in zip(range(normal_rep_size), fileN2):
                    normalFile = str(N1line[0:]).strip()
                    tumorFile = str(N2line[0:]).strip()
                    mcf_val = mcf_for_two(str(N1line[0:]).strip(), str(N2line[0:]).strip(), start, end) 
                    df_distMat[normal_id_array[i]][normal_id_array[j]] = mcf_val
             
    with open(file2,'rb') as fileT1:
        for  i, T1line in zip(range(tumor_rep_size), fileT1):
            with open(file2,'rb') as fileT2:
                for j, T2line in zip(range(tumor_rep_size), fileT2):
                    normalFile = str(T1line[0:]).strip()
                    tumorFile = str(T2line[0:]).strip()
                    mcf_val = mcf_for_two(str(T1line[0:]).strip(), str(T2line[0:]).strip(), start, end) 
                    df_distMat[tumor_id_array[i]][tumor_id_array[j]] = mcf_val

    return df_distMat, id_array_b

    #distMat = np.zeros(shape=(extended_dist_mat_size, extended_dist_mat_size))

    #with open(file1,'rb') as fileN:
    #    for i, Nline in zip(range(normal_rep_size), fileN):
     #       with open(file2,'rb') as fileT:
      #          for j, Tline in zip(range(tumor_rep_size), fileT):
       #             normalFile = str(Nline[0:]).strip()
        #            tumorFile = str(Tline[0:]).strip()
                    #print("normal file befor runForTwo = ", normalFile)
                    #print("tumor file befor runForTwo = ", tumorFile)
                    
         #           mcf_val = mcf_for_two(str(Nline[0:]).strip(), str(Tline[0:]).strip(), start, end)    
          #          distMat[i][j + normal_rep_size] = mcf_val
           #         distMat[j + normal_rep_size][i] = mcf_val

    #with open(file1,'rb') as fileN1:
     #   for i, N1line in zip(range(normal_rep_size), fileN1):
      #      with open(file1,'rb') as fileN2:
       #         for j, N2line in zip(range(normal_rep_size), fileN2):
        #            normal_1_File = str(N1line[0:]).strip()
         #           normal_2_File = str(N2line[0:]).strip()
                    #print("normal file befor runForTwo = ", normalFile)
                    #print("tumor file befor runForTwo = ", tumorFile)
                        
             #       distMat[i][j] = mcf_for_two(str(N1line[0:]).strip(), str(N2line[0:]).strip(), start, end)
   
   # with open(file2,'rb') as fileT1:
    #    for i, T1line in zip(range(tumor_rep_size), fileT1):
     #       with open(file2,'rb') as fileT2:
      #          for j, T2line in zip(range(tumor_rep_size), fileT2):
       #             tumor_1_File = str(T1line[0:]).strip()
        #            tumor_2_File = str(T2line[0:]).strip()
                    #print("normal file befor runForTwo = ", normalFile)
                    #print("tumor file befor runForTwo = ", tumorFile)
                        
          #          distMat[i + normal_rep_size][j + normal_rep_size] = mcf_for_two(str(T1line[0:]).strip(), str(T2line[0:]).strip(), start, end)


    #return distMat, id_array_b


def compute_distance_matrix(file1, file2, start, end):
    
    # Creates a list containing 5 lists, each of 8 items, all set to 0
    with open(file1) as f:
        normal_rep_size =  sum(1 for _ in f)
    with open(file2) as f:
        tumor_rep_size =  sum(1 for _ in f)
    distMat = np.zeros(shape=(normal_rep_size, tumor_rep_size))

    with open(file1,'rb') as fileN:
        for i, Nline in zip(range(normal_rep_size), fileN):
            with open(file2,'rb') as fileT:
                for j, Tline in zip(range(tumor_rep_size), fileT):
                    normalFile = str(Nline[0:]).strip()
                    tumorFile = str(Tline[0:]).strip()
                    #print("normal file befor runForTwo = ", normalFile)
                    #print("tumor file befor runForTwo = ", tumorFile)
                        
                    distMat[i][j] = mcf_for_two(str(Nline[0:]).strip(), str(Tline[0:]).strip(), start, end)
    return distMat
