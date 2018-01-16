
#!/usr/bin/env python
from __future__ import print_function
import csv
import sys, getopt
import numpy
import random
import sklearn

from sklearn.metrics import roc_auc_score


from ortools.graph import pywrapgraph
from apgl.graph import *


from mcflib import find_overal_intersect
from mcflib import compute_distance_matrix, compute_distance_matrix_sim
from simulation import simulate
from utils import *


def is_valid_file(parser, arg):
    """
        Check if arg is a valid file that already exists on the file system.
        
        Parameters
        ----------
        parser : argparse object
        arg : str
        
        Returns
        -------
        arg
        """
    arg = os.path.abspath(arg)
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg


def get_parser():
    import argparse


    parser = argparse.ArgumentParser(prog='MCFDiff')
    
    parser.add_argument('--version', action='version', version='MCFDiff-0.1.0')
    
    #subparsers = parser.add_subparsers()
    
    parser.add_argument("-run", "--run",
                        dest="run",
                        default=0,
                        type=int,
                        help="define type of process: simulation(0) or realdata(1) (sim/real)")
    parser.add_argument("-rep_size", "--rep_size", "-repSize",
                        dest="rep_size",
                        default=3,
                        type=int,
                        help="Number of replicate within each group for simulation run")
    parser.add_argument("-start", "--start",
                        dest="start",
                        default=10,
                        type=int,
                        help="start position of process")
    parser.add_argument("-end", "--end",
                        dest="end",
                        default=10,
                        type=int,
                        help="end position of process")
    parser.add_argument("-step", "--step",
                        dest="step",
                        default=100,
                        type=int,
                        help="moving window step size")
    parser.add_argument("-w", "--w",
                        dest="windowSize",
                        default=500,
                        type=int,
                        help="window size")
    parser.add_argument("-i1", "--i1",
                        dest="inputfileDir1",
                        #  type=lambda x: is_valid_file(parser, x),
                        help="A text file including list of input files name (Normal) in same directory")

    parser.add_argument("-i2", "--i2",
                        dest="inputfileDir2",
                        #  type=lambda x: is_valid_file(parser, x),
                        help="A text file including list of input files name (Tumor) in same directory")
    parser.add_argument("-o", "--o", "-output", "--output",
                        dest="outputfileDir",
                        #  type=lambda x: is_valid_file(parser, x),
                        help="outputfile directory to write outputs")


#parser.add_argument("-q", "--quiet",action="store_false",dest="verbose",default=True,help="don't print status messages to stdout")
    return parser


if __name__ == "__main__":
    
    parser = get_parser()
    result_args = parser.parse_args()
    
    inputFileDirNormal = result_args.inputfileDir1
    inputFileDirTumor = result_args.inputfileDir2
    outputDir = result_args.outputfileDir+'/distance.txt'
    open(outputDir, 'w').close()
    
    print(result_args.start, result_args.end, result_args.step, inputFileDirNormal, inputFileDirTumor, outputDir)
    run_type = result_args.run

    if run_type == 0:
        print("simulation start")
        fileN = open(inputFileDirNormal,'rb')
        line = fileN.readline()
        patterns = read_methyl_patterns(str(line).strip())
        print(" patterns length = ", len(patterns))
        selectedPatterns = within_region_patterns(patterns, result_args.start, result_args.end)
        print("selected pattern len = ", len(selectedPatterns))
        componentList = getComponentList(selectedPatterns)
        print(" comp list length = ", len(componentList))

        is_DMR_array = []
        dist_mean_array = []

        for comp in componentList:
            normalSelected = within_region_patterns(selectedPatterns, comp.start, comp.end)

            #rand = random.randint(0,101)
            rand = random.random()

            if rand < 0.50:
                is_DMR = 1
            else:
                is_DMR = 0

            normal_rep_pats, tumor_rep_pats = simulate(normalSelected, result_args.rep_size, is_DMR)


          #  with open(inputFileDirTumor,'rb') as fileT:
           #     for file_name, tumor_list in zip(fileT, tumor_rep_pat):
            #         with open(file_name, 'w') as f:
             #           for _string in tumor_list:
              #              f.write(str(_string) + '\n')

            distMat = compute_distance_matrix_sim(normal_rep_pats, tumor_rep_pats, comp.start, comp.end)


            dist_mean = distMat.mean()
            print("is_DMAR = " + str(is_DMR) + " , " + str(dist_mean))
            dist_mean_array.append(dist_mean)
            is_DMR_array.append(is_DMR)
                         
            with file(outputDir, 'a+') as outfile:
                outfile.write('component ' + str(comp.cid)+ ' = (' + str(comp.start)+ ', ' + str(comp.end) + ')\n')
                for pat in normalSelected:
                    outfile.write(str(pat)+'\n')
                if is_DMR == 1:
                    outfile.write('DMR component\n')
                else:
                    outfile.write('none-DMR component\n')
                numpy.savetxt(outfile, distMat, fmt='%1.0f')
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

        auc = sklearn.metrics.roc_auc_score(is_DMR_array, dist_mean_array)
        print(auc)

                    
    else:
        overalCompList = find_overal_intersect(inputFileDirNormal, inputFileDirTumor, result_args.start, result_args.end)

        for comp in overalCompList:
            distMat = compute_distance_matrix(inputFileDirNormal,  inputFileDirTumor, comp.start, comp.end)
            with file(outputDir, 'a+') as outfile:
                outfile.write('component ' + str(comp.cid)+ ' = (' + str(comp.start)+ ', ' + str(comp.end) + ')\n')
                numpy.savetxt(outfile, distMat, fmt='%1.0f')
                outfile.write('-----------------\n')



        print("comp List size = ", len(overalCompList))

