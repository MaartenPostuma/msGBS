import pandas as pd
import numpy as np
import csv
import os
import argparse
import datetime
import subprocess
import time
from warnings import simplefilter

start_time = time.time()

# THIS script is adjusted for the new msGBS Snakemake pipeline 2024

# It uses the rootmono read mapping to execute the cluster filtering

### Default filter parameters : ###
#   Filter 1 : 4  # effectively this means that a minimum of 4 mapped reads (total over all samples) are needed
#                         to ratain a cluster in the mapping database
#   Filter 2 : 15 # effectively this means that a cluster is filtered from the mapping data base if the  Total
#                      mapped reads (of all mono's combined) > (reads target mono + (reads target mono / 15))
#                      This taxonomic uninformative clusters.
#   Filter 3 : 1000 # min reads per samples to be retained (applied later in this script)

# names of mono samples         : rootmono1, rootmono2, rootmono3 ...
#
# names of reference contigs    : rootmono1_1, rootmono1_2, etc...

simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
pd.set_option('display.max_columns', 1000)
pd.options.mode.chained_assignment = None # getting rid of som error messages

def parse_args():
    """Pass command line arguments"""
    parser = argparse.ArgumentParser(description='morph msGBS stats into species percentages (after filtering)')
    #input files
    parser.add_argument('-i', '--input',
                        help='tab delimited .csv input file as produced by msGBS_STATS.py')
    parser.add_argument('-op', '--output_prefix',
                        help='general output location (folder) and first part of name that endswith _'
                             'example: /Users/NielsWagemaker/SGBS_NIELS/pilot_3/mapping_DINA_OKT_2024_')
    parser.add_argument('-f1', '--filter_1',
                        help='filtering parameter_1, minumum reads threshold',default=8)
    parser.add_argument('-f2', '--filter_2',
                        help='filtering parameter_2, sensitivity threshold',default=15)
    parser.add_argument('-f3', '--filter_3',
                        help='filtering parameter_2, minimum read count per sample',default=1000)
    parser.add_argument('-log', '--log',
                        help='log file, ')

    args = parser.parse_args()

    #create suffix argument from filter arguments
    parser.add_argument('-os', '--output_suffix',
                        help='final part of final output endswith .csv',default='%s_%s_%s'%(args.filter_1,args.filter_2,args.filter_3))
    args = parser.parse_args()

    return args

def run_subprocess(cmd,args,log_message):
    "Run subprocess under standardized settings"
    #force the cmds to be a string.
    if len(cmd) != 1:
        cmd = [" ".join(cmd)]
    with open(args.log,'a') as log:
        log.write("now starting:\t%s\n"%log_message)
        log.write('%s\n'%(' '.join(cmd)))
        log.flush()
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, executable='bash', universal_newlines=True)
        stdout, stderr = p.communicate()
        stdout = stdout.replace('\r','\n')
        stderr = stderr.replace('\r','\n')
        if stdout:
            log.write('stdout:\n%s\n'%stdout)
        if stderr:
            log.write('stderr:\n%s\n'%stderr)
        log.write('\n')
    return 0

def calculate_ESA(args):
    #ESA = estimated species abundance
    cmd = ["echo "" > %s"%args.log]
    log = "STARTLOG"
    run_subprocess(cmd,args,log)

    logstring = str("input file (created by msGBS_STATS.py) : %s \n"%args.input)
    logstring += str("\n")
    logstring += str("output_prefix : %s\n"%args.output_prefix)
    logstring += str("output_suffix : %s\n"%args.output_suffix)
    logstring += str("\n")
    logstring += str("filter_1 : %s \n"%args.filter_1)
    logstring += str("filter_2 : %s \n"%args.filter_2)
    logstring += str("filter_3 : %s \n"%args.filter_3)
    logstring += str("\n")
    logstring += str("reading input file \n")
    logstring += str("time started : %s\n"% datetime.datetime.now())

    cmd = [""]
    log = logstring
    run_subprocess(cmd,args,log)

    print("calculate Estimated Species Abundance")
    print("reading input file")

    print("start loading .csv file                  :  --- %s seconds ---" % round((time.time() - start_time),0))

    df = pd.read_csv(args.input, sep='\t| |;', engine='python')

    logstring = str("ready \n\n")
    logstring += str("time finished : %s\n"% datetime.datetime.now())
    logstring += str("Start: remove duplicate rootmono samples")

    cmd = [""]
    log = logstring
    run_subprocess(cmd,args,log)

    df_rootmono = df.filter(regex='^Locus|^rootmono', axis=1)

    df_rootmono.set_index('Locus', inplace=True)
    df_sum_sample_rootmono = df_rootmono.sum(axis=0)
    df_rootmono_div = df_rootmono.div(df_sum_sample_rootmono.squeeze())
    df_rootmono_div['max']=df_rootmono_div.idxmax(axis=1)
    getridof=[]
    df.set_index('Locus', inplace=True)

    '''FILTERING of CLUSTERS in the mapping database based on mapping data of monoculture samples (rootmono): 
        if too many reads of a mono sample map to other mono reference.
     This is called cluster filtering (See also: Wagemaker et al. 2020 Figure 3)'''

    f0low = 0
    f1low = 0
    f2low = 0
    number = 1000000

    cmd = [""]
    log = "Cluster filtering \n"
    log += "\nSTART cluster Filtering\n"
    log += "Number of Clusters before cluster filtering : %s\n"%(len(df.index))
    run_subprocess(cmd,args,log)

    ### The below file is created to be able to asses which clusters are discarded and why.... are there mistakes in
    ### the mono's (due to pollution of samples swapping) or are they causing problemns due to high homology to closely
    ### related species). The file is later processed to get an overview (deleted_clusters_max_%s_%s_%s_processed.csv)

    print("start_cluster_filtering                  :  --- %s seconds ---" % round((time.time() - start_time),0))

    resultFyle = open((args.output_prefix + "TMP_Clusters_Target_vs_Reason_to_remove_%s_%s_%s.txt" % (args.filter_1,args.filter_2,args.filter_3)),'w')
    wr = csv.writer(resultFyle, dialect='excel')

    for index, row in df_rootmono_div.iterrows():
        locus = index
        max = row[-1]
        maxi = row[-1]+'_'
        number += 1
        SUM =row.filter(regex='rootmono').sum(axis=0)
        high = row.get(max)
        ID = locus.split("_")
        item = str(ID[0]) + '_' + max
        if maxi not in locus:
            getridof.append(locus)
            f0low +=1
            #ID = locus.split("_")
            wr.writerow(str(item))
        else:
            mono = locus.split('_')
            mono = mono[0]
            summie = df_sum_sample_rootmono[mono]
            test =(int(args.filter_1)/float(summie))
            if high < test:
                getridof.append(locus)
                f1low +=1
                wr.writerow(str(item))
            elif SUM > (high + (high / int(args.filter_2))):
                getridof.append(locus)
                f2low +=1
                wr.writerow(str(item))
    print("-")
    print("prefilter step : remove cluster if target mono signal < highest mono signal")
    print("clusters removed by prefilter         : %s"%f0low)
    print("clusters removed by filter f1 (=%s)    : %s"%((args.filter_1),f1low))
    print("clusters removed by filter f2 (=%s)    : %s"%((args.filter_2),f2low))
    print("-")
    print("finished cluster_filtering               :  --- %s seconds ---" % round((time.time() - start_time),0))

    cmd = [""]
    log = "prefilter step : remove cluster if target mono signal < highest mono signal"
    log += "clusters removed by prefilter   : %s\n"%f0low
    log += "clusters removed by filter f1   : %s\n"%f1low
    log += "clusters removed by filter f2   : %s\n"%f2low
    run_subprocess(cmd,args,log)

    resultFyle.close()
    resultFyle_to_process = str(args.output_prefix + "TMP_Clusters_Target_vs_Reason_to_remove_%s.txt" %(args.output_suffix))
    resultFyle_processed = str(args.output_prefix + "Data_1_Clusters_Target_vs_Reason_to_remove_%s_summed_per_species.txt" %(args.output_suffix))

    cmd = ["cat %s |sed 's/,//g'|cut -f1 -d '_'|sort|uniq -c > %s"%(resultFyle_to_process, resultFyle_processed)]

    log = "Finished cluster filtering. \n"
    log += "Resultfile removed clusters processed!\n"
    run_subprocess(cmd,args,log)

    df_mono_clusters_removed = pd.read_csv(resultFyle_processed, header=None, sep='\t| |,', engine='python')

    # Make a list of all the rootmono samples in the input .csv file  #

    newTargetSpecieslist = []

    for index,row in df_mono_clusters_removed.iterrows():
        #print(row[1])
        newTargetSpecieslist.append(str(row.iloc[1]))

    resultFyle_to_process_2 = str(args.output_prefix + "TMP2_Clusters_Target_vs_Reason_to_remove_%s.txt" %(args.output_suffix))
    df_mono_clusters_removed_all = pd.DataFrame()

    # remove filtered clusters from data #

    print("remove filtered clusters from all data   :  --- %s seconds ---" % round((time.time() - start_time),0))

    for Species in newTargetSpecieslist:
        cmd = ["cat %s | sed 's/,//g' |grep '^%s_' |cut -f2 -d '_'|sort|uniq -c > %s"%(resultFyle_to_process, Species, resultFyle_to_process_2)]
        log = "%s is processed"%(Species)
        run_subprocess(cmd,args,log)
        df_clusters_removed_TMP = pd.read_csv(resultFyle_to_process_2, sep='\t| |,', header=None,engine='python')
        column = list(df_clusters_removed_TMP.columns)
        column = column[1]
        df_clusters_removed_TMP = df_clusters_removed_TMP.set_index(column)
        df_mono_clusters_removed_all = pd.concat([df_mono_clusters_removed_all,df_clusters_removed_TMP], axis=1, join="outer")
        df_mono_clusters_removed_all = df_mono_clusters_removed_all.rename(columns={0:Species})

    print("processed resultFyle_to_process          :  --- %s seconds ---" % round((time.time() - start_time),0))

    df_mono_clusters_removed_all = df_mono_clusters_removed_all.transpose()

    df_mono_clusters_removed_all = df_mono_clusters_removed_all.fillna(0)

    csv_name = (args.output_prefix + "Data_2_Clusters_filtered_due_to_homology_to_%s.txt"%(args.output_suffix))
    df_mono_clusters_removed_all = df_mono_clusters_removed_all.round(0)
    df_mono_clusters_removed_all.to_csv(csv_name,sep='\t' , decimal =',')

    os.remove(str(args.output_prefix + "TMP_Clusters_Target_vs_Reason_to_remove_%s.txt" %(args.output_suffix)))
    os.remove(str(args.output_prefix + "TMP2_Clusters_Target_vs_Reason_to_remove_%s.txt" %(args.output_suffix)))

    cmd = [""]
    log += "resultfile removed clusters processed!\n"
    run_subprocess(cmd,args,log)

    ''' Write removed loci to file ; This can be removed after implementation of the below two sections'''

    df_removed = df[df.index.isin(getridof)]
    print("-")

    ''' Here i remove the taxonomic lowering clusters from the mapping dataset'''

    df = df.drop(getridof,axis=0)

    print("clusters removed from all data           :  --- %s seconds ---" % round((time.time() - start_time),0))

    ''' HERE make an output with How many clusters remained and the total number of reads mapped to a 
    rootmono after filtering'''

    cmd = [""]
    log = "Number of clusters retained : %s\n"%(len(df.index))
    run_subprocess(cmd,args,log)

    print("-")

    df_rootmono2 = df.filter(regex='^Locus|^rootmono',axis=1)
    df_rootmono2.sort_values(by='rootmono1',ascending=False)

    ''' EXPORT number of reads mapped to removed clusters '''

    df_removed_sum = df_removed.groupby('Species').sum()
    df_removed_sum_T = df_removed_sum.transpose()
    df_removed_sum_T['Removed_Total'] =df_removed_sum_T.sum(axis=1)

    ''' Continue analysis on remaining clusters '''

    # Here the total read mapping per species are culculated #

    df_sum = df.groupby('Species').sum()
    df_sum = df_sum.transpose()
    df_sum['Total'] = df_sum.sum(axis=1)

    ''' Intermezzo to include the retained number of reads to the removed reads '''

    df_removed_sum_T['Retained_Total'] = df_sum['Total']
    csv_name = (args.output_prefix + 'Data_3_READ_COUNT_removed_CLUSTERS_SUM_%s.csv'%(args.output_suffix))
    df_removed_sum_T.to_csv(csv_name,sep='\t' , decimal =',')
    df_sum = df_sum.reset_index()
    csv_name = (args.output_prefix + 'Data_4_SUM_%s.csv'%(args.output_suffix))
    df_sum.to_csv(csv_name,sep='\t' , decimal =',')

    # Filter minimum number of reads per sample: #
    print("start application of filter f3 (minimum of %s reads per sample) :"%args.filter_3)
    df_sum = df_sum[df_sum['Total'] > int(args.filter_3)]
    csv_name = (args.output_prefix + 'Data_5_SUM_MINREAD_FILTER_%s.csv'%(args.output_suffix))
    df_sum.to_csv(csv_name,sep='\t' , decimal =',')

    print("finished save final output               :  --- %s seconds ---" % round((time.time() - start_time),0))

    return args

def main():
    """main function loop"""
    # get command line arguments
    args = parse_args()
    args = calculate_ESA(args)
    current_dateTime = datetime.datetime.now()
    print("analysis executed :")
    print("     Date     : %s-%s-%s"%(current_dateTime.day,current_dateTime.month,current_dateTime.year))

if __name__ == '__main__':
    main()


