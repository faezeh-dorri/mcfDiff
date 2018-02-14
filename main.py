
#!/usr/bin/env python
from __future__ import print_function
import csv, os
import sys, getopt
import numpy 
import random
import sklearn
import subprocess as sp
from scipy import stats




from sklearn.metrics import roc_auc_score


from ortools.graph import pywrapgraph
from apgl.graph import *


from mcflib import find_overal_intersect, compute_extended_distance_matrix_from_file, compute_extended_distance_matrix_sim, write_extended_dis_mat, write_single_dist_mat, write_health_status, write_DMR_status, write_auc_pvalue_MiRKAT, extended_dist_mat_2_distributions
from mcflib import compute_distance_matrix, compute_distance_matrix_sim
from simulation import simulate
from utils import *
from pandas import DataFrame


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
                        help="define type of process: simulation(0) or realdata(1)  (sim/real)")
    parser.add_argument("-eval", "--eval",
                        dest="eval",
                        default=False,
                        type=bool,
                        help="define type of process: evaluate_param_effect (true/False) ")
    parser.add_argument("-method", "--method",
                        dest="method",
                        default='MiRKAT',
                        type=str,
                        help="define method of process: MiRKAT or  ttest or others ")
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
def patterns_2_p_values(selectedPatterns, componentList,full_filename_sim, output_file_info, output_file_health_status, output_file_p_value, output_file_extended_dist_mat, output_file_auc, output_file_MiRKAT):
    auc_all_list = []
    auc_simple_list =[]
    auc_moderate_list = []
    auc_hard_list = []

    thr = 0.15
    DMR_change_param = 0.8
    replicate_change_param = 0.03
    rep_size = result_args.rep_size

    is_DMR_array_all = []
    dist_mean_array_all = []
    case_array_all = []

    is_DMR_array_simple = []
    dist_mean_array_simple = []

    is_DMR_array_moderate = []
    dist_mean_array_moderate = []

    is_DMR_array_hard = []
    dist_mean_array_hard = []

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



        extended_dist_mat, id_array_b = compute_extended_distance_matrix_sim(normal_rep_pats, tumor_rep_pats, comp.start, comp.end)
        write_extended_dis_mat(normal_rep_pats, tumor_rep_pats, extended_dist_mat, comp, is_DMR, case, output_file_extended_dist_mat)

        open(output_file_extended_dist_mat_single, 'w').close()
        open(output_file_health_status, 'w').close()
        open(output_file_info, 'w').close()

        #print(extended_dist_mat)
        write_single_dist_mat(extended_dist_mat, output_file_extended_dist_mat_single)
        write_health_status(id_array_b, output_file_health_status)
        write_DMR_status(is_DMR, case, output_file_info)


        args =[output_file_extended_dist_mat_single, output_file_health_status, output_file_info, output_file_p_value]
        cmd = ["Rscript", full_filename_sim] + args
        sp.call(cmd)
        

    if len(is_DMR_array_all) >= 2:
        auc_all = sklearn.metrics.roc_auc_score(is_DMR_array_all, dist_mean_array_all)
        auc_all_list.append(auc_all)
    else:
        auc_all = 'NA'

    if len(is_DMR_array_simple) >= 2:
        auc_simple = sklearn.metrics.roc_auc_score(is_DMR_array_simple, dist_mean_array_simple)
        auc_simple_list.append(auc_simple)
    else:
        auc_simple = 'NA'

    if len(is_DMR_array_moderate) >= 2:
        auc_moderate = sklearn.metrics.roc_auc_score(is_DMR_array_moderate, dist_mean_array_moderate)
        auc_moderate_list.append(auc_moderate)
    else:
        auc_moderate = 'NA'

    if len(is_DMR_array_hard) >= 2:
        auc_hard = sklearn.metrics.roc_auc_score(is_DMR_array_hard, dist_mean_array_hard)
        auc_hard_list.append(auc_hard)
    else:
        auc_hard = 'NA'

    with open(output_file_auc, 'a+') as outfile:
        outfile.write(str(rep_size) + '\t'+ str(thr) + '\t' + str(DMR_change_param) + '\t' + str(replicate_change_param) + '\tALL\t')
        if len(auc_all_list) > 0:
            outfile.write(str(sum(auc_all_list)/len(auc_all_list)) + '\n')
        else:
            outfile.write('NA\n')

        outfile.write(str(rep_size) + '\t' + str(thr) + '\t' + str(DMR_change_param) + '\t' + str(replicate_change_param) + '\tSIM\t')
        if len(auc_simple_list) > 0:
            outfile.write(str(sum(auc_simple_list)/len(auc_simple_list)) + '\n')
        else:
            outfile.write('NA\n')

        outfile.write(str(rep_size) + '\t' + str(thr) + '\t' + str(DMR_change_param) + '\t' + str(replicate_change_param) + '\tMOD\t')
        if len(auc_moderate_list) > 0:            
            outfile.write(str(sum(auc_moderate_list)/len(auc_moderate_list)) + '\n')
        else:
            outfile.write('NA\n')

        outfile.write(str(rep_size) + '\t' + str(thr) + '\t' + str(DMR_change_param) + '\t' + str(replicate_change_param) + '\tHAR\t') 
        if len(auc_hard_list) > 0:              
            outfile.write(str(sum(auc_hard_list)/len(auc_hard_list)) + '\n')
        else:
            outfile.write('NA\n')


    pvalue_data = numpy.genfromtxt(output_file_p_value, delimiter="\t")
    pvalue_MiRKAT = pvalue_data[:,0]
    #pvalue_MMiRKAT = pvalue_data[:,1]
    #labels = 1 - pvalue_data[:,2]

    labels = 1 - pvalue_data[:,1]
    significance_level = 0.02
    auc_pvalue_MiRKAT = sklearn.metrics.roc_auc_score(labels, pvalue_MiRKAT)
    #auc_pvalue_MMiRKAT = sklearn.metrics.roc_auc_score(labels, pvalue_MMiRKAT)

    total_number_of_regions = len(pvalue_MiRKAT)
    type_I = 0
    type_II = 0
    for pval, lab  in zip(pvalue_MiRKAT, labels):
        if pval <= significance_level and lab == 0:
            type_I += 1
        if pval > significance_level and lab == 1:
            type_II +=1
    type_I_error = float(type_I) / float(total_number_of_regions)
    type_II_error = float(type_II) / float(total_number_of_regions)

    write_auc_pvalue_MiRKAT(auc_pvalue_MiRKAT, type_I_error, type_II_error, output_file_MiRKAT)
    


def evaluate_param_effect(selectedPatterns, componentList,full_filename_sim, method, output_file_auc_eval):
    random_iteration = 20

    thr_list = [0.15]
    DMR_change_param_list = [0.6, 0.8]
    replicate_change_param_list = [0.02, 0.04, 0.06, 0.08, 0.1]
    #DMR_change_param_list = [0.8]
    #replicate_change_param_list = [0.02]

    rep_size_list = [3, 5, 10, 20]

    for thr in thr_list:
        for DMR_change_param in DMR_change_param_list:
            for replicate_change_param in replicate_change_param_list:
                for rep_size  in rep_size_list:
                    auc_all_list = []
                    auc_simple_list =[]
                    auc_moderate_list = []
                    auc_hard_list = []

                    for iter in range(random_iteration):
                        is_DMR_array_all = []
                        dist_mean_array_all = []
                        case_array_all = []

                        is_DMR_array_simple = []
                        dist_mean_array_simple = []

                        is_DMR_array_moderate = []
                        dist_mean_array_moderate = []

                        is_DMR_array_hard = []
                        dist_mean_array_hard = []
                        
                      ###  is_DMR_array_all, is_DMR_array_simple, is_DMR_array_moderate, is_DMR_array_hard = patterns_2_is_DMR_arrays(selectedPatterns, componentList, rep_size, is_DMR, DMR_change_param, replicate_change_param, comp)

                        for comp in componentList:
                            normalSelected = within_region_patterns(selectedPatterns, comp.start, comp.end)

                            #rand = random.randint(0,101)
                            rand = random.random()

                            if rand < 0.50:
                                is_DMR = 1
                            else:
                                is_DMR = 0

                            normal_rep_pats, tumor_rep_pats, case = simulate(normalSelected, rep_size, is_DMR, thr, DMR_change_param, replicate_change_param)

                            

                            if method == 'MiRKAT':
                                extended_dist_mat, id_array_b = compute_extended_distance_matrix_sim(normal_rep_pats, tumor_rep_pats, comp.start, comp.end)
                                open(output_file_extended_dist_mat_single, 'w').close()
                                open(output_file_health_status, 'w').close()
                                open(output_file_info, 'w').close()

                                write_extended_dis_mat(normal_rep_pats, tumor_rep_pats, extended_dist_mat, comp, is_DMR, case, output_file_extended_dist_mat)
                                
                                write_single_dist_mat(extended_dist_mat, output_file_extended_dist_mat_single)
                                write_health_status(id_array_b, output_file_health_status)
                                write_DMR_status(is_DMR, case, output_file_info)



                                args =[output_file_extended_dist_mat_single, output_file_health_status, output_file_info, output_file_p_value, output_file_single_p_value]
                                cmd = ["Rscript", full_filename_sim] + args
                                sp.call(cmd)

                                single_p_value = open(output_file_single_p_value, 'r').readline().rstrip()

                                print("MiRKAT_p_ value = " + str(single_p_value))
                                score = 1 - float(single_p_value)

                            elif method == 'ttest':
                                extended_dist_mat, id_array_b = compute_extended_distance_matrix_sim(normal_rep_pats, tumor_rep_pats, comp.start, comp.end)
                                n_n_array, t_t_array, n_t_array = extended_dist_mat_2_distributions(extended_dist_mat, len(normal_rep_pats), len(tumor_rep_pats))

                                ttest_tscore, ttest_pvalue= stats.ttest_ind(n_n_array, n_t_array)
                                #print("ttest_tscore = " + str(ttest_tscore))
                                print("ttest_pvalue = " + str(2*ttest_pvalue))

                                #score = 1 - 2*ttest_pvalue
                                #score = abs(ttest_tscore)
                                score = ttest_pvalue


                            else:
                                distMat = compute_distance_matrix_sim(normal_rep_pats, tumor_rep_pats, comp.start, comp.end)
                                dist_mean = distMat.mean()
                                score = dist_mean
                     

                            dist_mean_array_all.append(score)
                            is_DMR_array_all.append(is_DMR)
                            case_array_all.append(case)

                            if case == 'simple':
                                dist_mean_array_simple.append(score)
                                is_DMR_array_simple.append(is_DMR)
                            else: 
                                if case == 'moderate':
                                    dist_mean_array_moderate.append(score)
                                    is_DMR_array_moderate.append(is_DMR)
                                else:
                                    dist_mean_array_hard.append(score)
                                    is_DMR_array_hard.append(is_DMR)

                            '''
                            with file(output_file_distance, 'a+') as outfile:
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
                                outfile.write('-----------------\n')'''
                            
                        print("len_all = " + str(len(is_DMR_array_all)))
                        print("len_simple = " + str(len(is_DMR_array_simple)))
                        print("len_moderate = " + str(len(is_DMR_array_moderate)))
                        print("len_hard = " + str(len(is_DMR_array_hard)))

                        if len(is_DMR_array_all) >= 2:
                            auc_all = sklearn.metrics.roc_auc_score(is_DMR_array_all, dist_mean_array_all)
                            auc_all_list.append(auc_all)
                        else:
                            auc_all = 'NA'

                        if len(is_DMR_array_simple) >= 2:
                            auc_simple = sklearn.metrics.roc_auc_score(is_DMR_array_simple, dist_mean_array_simple)
                            auc_simple_list.append(auc_simple)
                        else:
                            auc_simple = 'NA'

                        if len(is_DMR_array_moderate) >= 2:
                            auc_moderate = sklearn.metrics.roc_auc_score(is_DMR_array_moderate, dist_mean_array_moderate)
                            auc_moderate_list.append(auc_moderate)
                        else:
                            auc_moderate = 'NA'

                        if len(is_DMR_array_hard) >= 2:
                            auc_hard = sklearn.metrics.roc_auc_score(is_DMR_array_hard, dist_mean_array_hard)
                            auc_hard_list.append(auc_hard)
                        else:
                            auc_hard = 'NA'

                    with open(output_file_auc_eval, 'a+') as outfile:
                        outfile.write(str(rep_size) + '\t'+ str(thr) + '\t' + str(DMR_change_param) + '\t' + str(replicate_change_param) + '\tALL\t')
                        if len(auc_all_list) > 0:
                            outfile.write(str(sum(auc_all_list)/len(auc_all_list)) + '\n')
                        else:
                            outfile.write('NA\n')

                        outfile.write(str(rep_size) + '\t' + str(thr) + '\t' + str(DMR_change_param) + '\t' + str(replicate_change_param) + '\tSIM\t')
                        if len(auc_simple_list) > 0:
                            outfile.write(str(sum(auc_simple_list)/len(auc_simple_list)) + '\n')
                        else:
                            outfile.write('NA\n')

                        outfile.write(str(rep_size) + '\t' + str(thr) + '\t' + str(DMR_change_param) + '\t' + str(replicate_change_param) + '\tMOD\t')
                        if len(auc_moderate_list) > 0:            
                            outfile.write(str(sum(auc_moderate_list)/len(auc_moderate_list)) + '\n')
                        else:
                            outfile.write('NA\n')

                        outfile.write(str(rep_size) + '\t' + str(thr) + '\t' + str(DMR_change_param) + '\t' + str(replicate_change_param) + '\tHAR\t') 
                        if len(auc_hard_list) > 0:              
                            outfile.write(str(sum(auc_hard_list)/len(auc_hard_list)) + '\n')
                        else:
                            outfile.write('NA\n')




if __name__ == "__main__":
    
    parser = get_parser()
    result_args = parser.parse_args()
    
    inputFileDirNormal = result_args.inputfileDir1
    inputFileDirTumor = result_args.inputfileDir2
    method = result_args.method
    output_file_auc_eval = result_args.outputfileDir+'/auc_eval.txt'
    output_file_auc = result_args.outputfileDir+'/auc.txt'
    output_file = result_args.outputfileDir+'/final.txt'

    output_file_extended_dist_mat = result_args.outputfileDir + '/extended_dist_mat.txt'
    output_file_extended_dist_mat_single = result_args.outputfileDir + '/extended_dist_mat_single.txt'
    output_file_health_status = result_args.outputfileDir + '/health_status.txt'
    output_file_info = result_args.outputfileDir + '/info.txt'
    output_file_p_value = result_args.outputfileDir + '/p_value.txt'
    output_file_single_p_value = result_args.outputfileDir + '/single_p_value.txt'
    output_file_auc_values = result_args.outputfileDir + '/AUC_values.txt'
    output_file_MiRKAT = result_args.outputfileDir + '/MiRKAT_p_value.txt'


    fileDir = os.path.dirname(os.path.realpath('__file__'))
    print(fileDir)

    #For accessing the file in the same folder
    filename_MiRKAT_sim = "MiRKAT_p_value.r"
    filename_MiRKAT_real = "MiRKAT_p_value_real.r"

    #readFile(filename)

    #For accessing the file in a folder contained in the current folder
    full_filename_sim = os.path.join(fileDir, filename_MiRKAT_sim)
    full_filename_real = os.path.join(fileDir, filename_MiRKAT_real)


    open(output_file_auc, 'w').close()
    open(output_file, 'w').close()
    output_file_distance = result_args.outputfileDir+'/distance.txt'
    open(output_file_distance, 'w').close()
    open(output_file_p_value, 'w').close()


    
    print(result_args.start, result_args.end, result_args.step, inputFileDirNormal, inputFileDirTumor, output_file_distance, output_file_auc)
    run_type = result_args.run
    eval_type = result_args.eval

    if run_type == 0:
        print("simulation start")
        fileN = open(inputFileDirNormal,'r')
        line = fileN.readline()
        print(str(line).strip())
        patterns = read_methyl_patterns(str(line).strip())
        print(" patterns length = ", len(patterns))
        selectedPatterns = within_region_patterns(patterns, result_args.start, result_args.end)
        print("selected pattern len = ", len(selectedPatterns))
        componentList = getComponentList(selectedPatterns)
        print(" comp list length = ", len(componentList))

        if eval_type == False:
            patterns_2_p_values(selectedPatterns, componentList, full_filename_sim, output_file_info, output_file_health_status, output_file_p_value, output_file_extended_dist_mat, output_file_auc, output_file_MiRKAT)

        if eval_type == True:
            evaluate_param_effect(selectedPatterns, componentList, full_filename_sim, method, output_file_auc_eval)
    
    else:
        overalCompList = find_overal_intersect(inputFileDirNormal, inputFileDirTumor, result_args.start, result_args.end)
        print("length of compList = " + str(len(overalCompList)))
        score_list = []
        for comp in overalCompList:
            is_DMR = 'unknown'

            extended_dist_mat, id_array_b = compute_extended_distance_matrix_from_file(inputFileDirNormal,  inputFileDirTumor, comp.start, comp.end)
            #extended_dist_mat, id_array_b = compute_extended_distance_matrix_sim(normal_rep_pats, tumor_rep_pats, comp.start, comp.end)
            #write_extended_dis_mat(normal_rep_pats, tumor_rep_pats, extended_dist_mat, comp, is_DMR, case, output_file_extended_dist_mat)

            with open(inputFileDirNormal) as f:
                normal_rep_size =  sum(1 for _ in f)
            with open(inputFileDirTumor) as f:
                tumor_rep_size =  sum(1 for _ in f)

            open(output_file_extended_dist_mat_single, 'w').close()
            open(output_file_health_status, 'w').close()
            open(output_file_info, 'w').close()

            #extended_dist_mat.to_csv(output_file, sep='\t')
            write_single_dist_mat(extended_dist_mat, output_file_extended_dist_mat_single)
            write_health_status(id_array_b, output_file_health_status)
            #write_DMR_status(is_DMR, case, output_file_info)


            args =[output_file_extended_dist_mat_single, output_file_health_status, output_file_p_value, output_file_single_p_value]
            cmd = ["Rscript", full_filename_real] + args
            sp.call(cmd)


            n_n_array, t_t_array, n_t_array = extended_dist_mat_2_distributions(extended_dist_mat, normal_rep_size, tumor_rep_size)

            ttest_tscore, ttest_pvalue= stats.ttest_ind(n_n_array, n_t_array)
            #print("ttest_tscore = " + str(ttest_tscore))
            print("ttest_pvalue = " + str(2*ttest_pvalue))

            #score = 1 - 2*ttest_pvalue
            score = ttest_pvalue
            score_list.append(score)
            
        all = []
        qbfile = open(output_file_p_value,"r")

        for score, comp, aline in zip(score_list, overalCompList, qbfile.readlines()):
            values = str(aline).strip() + '\t' + str(score) + '\t' + str(comp.start) + '\t' + str(comp.end)
            all.append(values)

        qbfile.close()
            


        with open(output_file, 'w') as outfile:
            outfile.write('pvalue_MiRKAT\tttest_pvalue\tcompStart\tcompEnd\n')
            for line in all:
                outfile.write(str(line)+'\n')
            


        print("comp List size = ", len(overalCompList))

