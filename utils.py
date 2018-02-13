
#!/usr/bin/env python
from __future__ import print_function
import csv
import sys, getopt
import numpy


from ortools.graph import pywrapgraph
from apgl.graph import *



class methylPattern(object):
    def __init__(self,chr,start,end,cid,pid,abundance,mPat,regions,group, patId):
        self.chr=chr
        self.start=start
        self.end=end
        self.cid=cid
        self.pid=pid
        self.abundance=abundance
        self.mPat=mPat
        self.regions=regions
        self.group=group
        self.patId = patId
    
    
    def __str__(self):
        return 'patId = ' + str(self.patId) + ' start =' + str(self.start) + ' mPat = ' + str(self.mPat) + ' abd = ' + str(self.abundance)
    
    def __eq__(self, other):
        return (self.start == other.start) and (self.end == other.end) and (self.mPat == other.mPat)

    def setId(self, vertId):
        self.patId = vertId
    
    def getId(self):
        return self.patId






class Range:
    start = 0
    end = 0
    def __init__(self, start, end):
        self.start = start
        self.end = end
    
    def __str__(self):
        return '(' + str(self.start) + ', ' + str(self.end) + ')'



class Methyl:
    def __init__(self, pos, methyl):
        self.pos = pos
        self.methyl = methyl
    
    def __str__(self):
        return '(' + str(self.pos) + ', ' + str(self.methyl) + ')'



class Component:
    def __init__(self, cid, start , end):
        self.cid = cid
        self.start = start
        self.end = end
    
    def __str__(self):
        return '(' + str(self.cid) + ', ' + str(self.start) + ',' + str(self.end) + ')'


def findOverlapRange(metPat1, metPat2):
    overlapRange = Range(max(metPat1.start, metPat2.start), min(metPat1.end, metPat2.end))
    return overlapRange
# return range(overlapRange.start, overlapRange.end + 1)

def getSubPattern(metPat, overlapRange):
    subPat = []
    if str(metPat.mPat) == "*":
        return subPat
    
    tokens = str(metPat.mPat).split(",")
    
    for token in tokens:
        
        #print("token = ", token)
        element = token.split(":")
        #print("element = ", element)
        pos = int(element[0]) + metPat.start
        methyl = element[1]
        mpat = Methyl(pos, methyl)
        if mpat.pos >= overlapRange.start and mpat.pos <= overlapRange.end:
            subPat.append(mpat)
    
    #for subPat in selectedPat:
    #    metString.append(sub)

    return subPat;



def getJaccardDistance(metPat1, metPat2, start, end):
    alpha = 0.5
    overlapRange = findOverlapRange(metPat1, metPat2)
    if overlapRange == []:
        return 1
    
    finalRange = Range(max(start, overlapRange.start),min(end, overlapRange.end))


    subPat1 = getSubPattern(metPat1, finalRange)
    subPat2 = getSubPattern(metPat2, finalRange)
    
    if len(subPat1) == 0 or len(subPat2) == 0:
        print("Jaccard: a mission value (*) is detected")
        return 1
    
    if len(subPat1) != len(subPat2):
        print("Jaccard: there is some issue, missing values exist")
        return 1

    #print("subPat1 = ", subPat1[0])
    #print("subPat2 = ", subPat2[0])
    diffCount = sum(1 for a, b in zip(subPat1, subPat2) if a.pos == b.pos and a.methyl != b.methyl)
    missCount = sum(1 for a, b in zip(subPat1, subPat2) if a.pos != b.pos)
    #print("diff = ", diffCount)


    cost = (diffCount + alpha * missCount)/len(subPat1)
    return cost




def read_methyl_patterns(fileName):
    methylPatterns = []
    patId =0
    # 0 = chr, 1 = start, 2 = end, 3 = cid, 4 = pid, 5 = abundance, 6 = methylpat, 7 = regions , 8 = group
    with open(fileName,'r') as tsvin:
        tsvin=csv.reader(tsvin,delimiter='\t')
        next(tsvin)
        for row in tsvin:
            patId +=1
            #print("row6 = " , row[6])
            methylPatterns.append(methylPattern(row[0],int(row[1]),int(row[2]),row[3],row[4],int(float(row[5])*100),row[6],row[7], "normal", patId))

    return methylPatterns


def findComponentRange(patterns, cid):
    selectedPat = []
    for pat in patterns:
        if pat.cid == cid:
            selectedPat.append(pat)
    
    mm = min(pat.start for pat in selectedPat)
    MM = max(pat.end for pat in selectedPat)
    componentRange = Range(mm, MM)
    return componentRange

def within_region_patterns(patterns, start, end):
    patternList = []
    for pat in patterns:
        if pat.start > end or pat.end < start:
            continue
        patternList.append(pat)
    return patternList

def isInList(cid, compList):
    for comp in compList:
        if cid == comp.cid:
            print("Cid = ", cid , " is already in list" )
            return True
        else:
            continue
    return False

def getComponentList(patterns):
    unique_comp_list = []
    for pat in patterns:
        if not isInList(pat.cid, unique_comp_list):
            compRange = findComponentRange(patterns, pat.cid)
            cidComp = Component(pat.cid, compRange.start,compRange.end )
            unique_comp_list.append(cidComp)
            print( "cid = ", cidComp.cid, cidComp.start, cidComp.end, " is added to unique list" )
    return unique_comp_list




def IntersectofTwoComponentLists(componentList1, componentList2):
#assume both inputs are sorted based on their start
    cid = 0
    intersectList = []
    i = 0
    j = 0
    while i < len(componentList1) and j < len(componentList2):
        a = componentList1[i]
        b = componentList2[j]
        if a.start > b.end:
            j += 1
            continue

        if a.end < b.start:
            i += 1
            continue

        intersectList.append(Component(cid, max(a.start, b.start), min(a.end, b.end)))
        cid += 1
        if a.end < b.end:
            i += 1
        else:
            j += 1
    return intersectList







