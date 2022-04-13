#######################################################################################################################
#######################################################################################################################
#########    GET AVERAGE LFQ FROM INDIVIUAL LFQ AND OUTPUT MATRIX>                                         ############
#########
#########    
#########                                                                                                
#########   -by Travis Schrank                                                                  
#######################################################################################################################

#######################################################################################################################
######### *** Notice : Several Parameters must be set manully in the code for it to work             ********##########
######### *** These should be marked with >>>!!!ACTION ITEM!!!<<< flag in the comments above them    ********##########
#######################################################################################################################

#### Import all relavent packages ####

import sys
import os
#### for math.log(x) and math.exp(x)
import math
#### The search path includes the directory that the top level program is in ###
#import time 
####  time:   not sure I need this
import datetime
#### datetime:this module lets me subtract dates
#### form is datetime.date(YYYY,MM,DD)=OBJ; Diff = OBJ-OBJ2
import random
### random.shuffle(List) shuffles the top level of a list (only)
import zlib
import gzip
### I need this go read in database files
import string
import csv
import re
import pandas
import pandas as pd
import scipy
#from astropy import stats as aps
from scipy import stats as sps
import numpy as np
import pickle
import re


###### import pickled mutation data

manifest = open("/home/parke/ssdb2/UNCseq_oropharynx_fastqs/oropharynx_fastq_links/oropharynx_tumor_manifest.txt")

runs = []

for i in manifest.readlines():
    #print(i.split()[0])
    runs.append(i.split()[0])


sample_l = []
    
human_length_l = []
human_aligns_l = []
human_not_aligned_l = []
    
other_hpvs_length_l = []
other_hpvs_aligns_l = []
other_hpvs_not_aligned_l = []

hpv16_length_l = []
hpv16_aligns_l = []
hpv16_not_aligned_l = []

merged_out_lines = []

for i in runs:
    
    my_line_list = []
    
    human_length = 0
    human_aligns = 0
    human_not_aligned = 0
    
    other_hpvs_length = 0
    other_hpvs_aligns = 0
    other_hpvs_not_aligned = 0

    hpv16_length = 0
    hpv16_aligns = 0
    hpv16_not_aligned = 0
    
    data = open("/home/parke/ssdb2/UNCseq_oropharynx_fastqs/ViFi_Run/" + i + "/output.sorted.bam.idxstat.txt")

    for ii in data.readlines():
        if ii[0:1] == "*":
            test = 1
        else:
            merged_out_lines.append(i[0:-2] + "\t" + ii)
            
            my_line_list.append(ii)
            split_line = ii.split()

            if split_line[0][0:3] == "chr":
                    human_length = human_length + int(split_line[1])
                    human_aligns = human_aligns + int(split_line[2])
                    human_not_aligned = human_not_aligned + int(split_line[3])
            elif split_line[0][0:10] == "hpv16ref_1":
                    hpv16_length= hpv16_length + int(split_line[1])
                    hpv16_aligns = hpv16_aligns + int(split_line[2])
                    hpv16_not_aligned = hpv16_not_aligned + int(split_line[3])
            else:
                    other_hpvs_length = other_hpvs_length + int(split_line[1])
                    other_hpvs_aligns = other_hpvs_aligns + int(split_line[2])
                    hpv16_not_aligned = hpv16_not_aligned + int(split_line[3])
            
    human_length_l.append(human_length)
    human_aligns_l.append(human_aligns)
    human_not_aligned_l.append(human_not_aligned)
    
    other_hpvs_length_l.append(other_hpvs_length)
    other_hpvs_aligns_l.append(other_hpvs_aligns)
    other_hpvs_not_aligned_l.append(other_hpvs_not_aligned)

    hpv16_length_l.append(hpv16_length)
    hpv16_aligns_l.append(hpv16_aligns)
    hpv16_not_aligned_l.append(hpv16_not_aligned)

    sample_l.append(i)

    data.close()

        

            









out_data = pandas.DataFrame({
    
    'sample':sample_l,
    
    'human_length':human_length_l,
    'human_aligns':human_aligns_l,
    'human_not_aligned':human_not_aligned_l,
    
    'other_hpvs_length':other_hpvs_length_l,
    'other_hpvs_aligns':other_hpvs_aligns_l,
    'other_hpvs_not_aligned':other_hpvs_not_aligned_l,

    'hpv16_length':hpv16_length_l,
    'hpv16_aligns':hpv16_aligns_l,
    'hpv16_not_aligned':hpv16_not_aligned_l,

    })






                    
out_data.to_csv("/home/parke/Desktop/UNCseq_cervix_op/oropahrynx_cohort/4_batch_idx_stats/idx_stats_summary.csv", na_rep='NA', index = False)
            
with open("merged_idx_stats.txt", 'a') as the_file:
    for iii in merged_out_lines:
        the_file.write(iii)



