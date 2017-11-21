
#!/usr/bin/env python
from __future__ import print_function
import csv
import sys, getopt
import numpy


from ortools.graph import pywrapgraph
from apgl.graph import *

from mcflib import find_overal_intersect
from mcflib import compute_distance_matrix



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

    overalCompList = find_overal_intersect(inputFileDirNormal, inputFileDirTumor, result_args.start, result_args.end)

    for comp in overalCompList:
        distMat = compute_distance_matrix(inputFileDirNormal,  inputFileDirTumor, comp.start, comp.end)
        with file(outputDir, 'a+') as outfile:
            outfile.write('component ' + str(comp.cid)+ ' = (' + str(comp.start)+ ', ' + str(comp.end) + ')\n')
            numpy.savetxt(outfile, distMat, fmt='%1.0f')
            outfile.write('-----------------\n')



    print("comp List size = ", len(overalCompList))

