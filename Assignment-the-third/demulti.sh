#!/usr/bin/env python

import argparse
import itertools
import gzip
import re


def get_args():
    '''adds argparse arguments to run script from the command line'''
    parser = argparse.ArgumentParser(description="takes arguments for input files")
    parser.add_argument("-r1", "--read1", help="read1 filename")
    parser.add_argument("-r2", "--read2", help="read2 filename")
    parser.add_argument("-r3", "--read3", help="read3 filename")
    parser.add_argument("-r4", "--read4", help="read4 filename")
    parser.add_argument("-i", "--indexfile", help="index info filename")
    parser.add_argument("-c", "--cutoffscore", help="cutoff score for average index quality")
    parser.add_argument("-n", "--numberrecords", help="number of records in files", type=int)
    return parser.parse_args()

args = get_args()

read1 = args.read1
read2 = args.read4
index1 = args.read2
index2 = args.read3
index_file = args.indexfile
cutoff_score = args.cutoffscore
number_records = args.numberrecords


def convert_phred(letter):
    '''Takes a single letter phred score and coverts it to a number, phred-33 only'''
    qscore = ord(letter) - 33
    return qscore

def index_score(phred_score, cutoff_score):
    '''returns True if average qscore of index is above cutoff score, or False if not '''
    total = 0
    for letter in phred_score:
        score = convert_phred(letter) 
        total += score
    average = total/len(phred_score)
    if average >= float(cutoff_score):
        return True
    else:
        return False

def create_index_dicts(index_file):
    '''opens file that contains indexes, make keys from indexes with empty string as value'''
    index_dict = {}
    line_counter = 0
    with open(index_file) as fh:
        for line in fh:
            line = line.strip()
            if line_counter == 0:
                line_counter += 1
            else:
                split_line = line.split('\t')
                index_key = split_line[-1]
                index_value = split_line[0:4]
                index_dict.setdefault(index_key, index_value)
    return index_dict

#creates index dictionary
index_dict = create_index_dicts(index_file)

def rev_comp(index):
    '''take an index and returns the reverse complement, accounts for base calls of N'''
    rev_comp_index = ''
    ref_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    for i in index:
        rev_comp_index = ref_dict[i] + rev_comp_index
    return rev_comp_index

#creates dictionary of reverse complements of indexes in index_dict
rev_comp_dict = {}
for key in index_dict:
    rev_comp_dict.setdefault(rev_comp(key), '')

#creates a list of all indexes and reverse complements to use in creating permutations dictionary
all_indexes = []
for key in index_dict.keys():
    all_indexes.append(key)

def create_perm_dict(all_indexes):
    '''uses itertools to create dictionary of all possible permutations of the indexes in the index dictionary and complement dictionary'''
    perm_dict = {}
    for item in itertools.product(all_indexes, repeat=2):
        start, end = item
        perm_dict.setdefault(start + "-" + end, 0)
    return perm_dict

#creates permutations dictionary in format 'index-index'
perm_dict = create_perm_dict(all_indexes)

#creates output files names for index matched reads using index dictionary
split_name = index1.split('_')
file_names = {}
for index, value in index_dict.items():
    file_name_1 = value[0] + '_' + value[1] + '_' + value[2] + '_' + value[3] + '_' + split_name[2] + '_' + 'R1' + '_' + '001.fastq'
    file_name_2 = value[0] + '_' + value[1] + '_' + value[2] + '_' + value[3] + '_' + split_name[2] + '_' + 'R2' + '_' + '001.fastq'
    file_names.setdefault(index, [file_name_1, file_name_2])

#opens output files using file_names dictionary
for key in file_names:
    for i in range(0, 2):
        file_names[key][i] = open(file_names[key][i], "w")
      
#open hardcoded files here 
hop_r1 = open('index_hopped_r1.fastq', 'w')
hop_r2 = open('index_hopped_r2.fastq', 'w')
bad_egg_r1 = open('bad_eggs_r1.fastq', 'w')
bad_egg_r2 = open('bad_eggs_r2.fastq', 'w')

def sort_indexes(index_dict, rev_comp_dict, perm_dict, read1, read2, index1, index2):
    '''uses fastq files containing reads 1-4 and sorts the records into new files for each index, also sorts out index hopped reads
    and reads with poor quality/unknown indexes into their own files'''
    num_matched_reads = 0
    num_hopped_reads = 0
    num_unknown = 0
    num_bad_quality = 0
    r1 = gzip.open(read1, 'tr')
    r2 = gzip.open(read2, 'tr')
    i1 = gzip.open(index1, 'tr')
    i2 = gzip.open(index2, 'tr')
#creates chunks of files to read as records
    for i in range(number_records):
        record_r1 = list(itertools.islice(r1, 4))
        record_r2 = list(itertools.islice(r2, 4))
        record_i1 = list(itertools.islice(i1, 4))
        record_i2 = list(itertools.islice(i2, 4))
        #strips new line characters
        record_r1 = [str(i[:-1]) for i in record_r1]
        record_r2 = [str(i[:-1]) for i in record_r2]
        record_i1 = [str(i[:-1]) for i in record_i1]
        record_i2 = [str(i[:-1]) for i in record_i2]
        #creates index pair to add to header lines
        index_pair = record_i1[1] + '-' + rev_comp(record_i2[1])

        #check indexes and sort into files
        #first checks if average quality score of index is above cutoff
        if index_score(record_i1[3], cutoff_score) is True and index_score(record_i2[3], cutoff_score) is True:
            #checks if indexes match
            if record_i1[1] in index_dict.keys() and record_i2[1] == rev_comp(record_i1[1]):
                #adds count to permutations dictionary
                perm_dict[index_pair] += 1
                for key in file_names:
                    if key == record_i1[1]:
                        file_names[key][0].write(record_r1[0] + ' ' + str(index_pair) + '\n')
                        file_names[key][0].write(record_r1[1] + '\n')
                        file_names[key][0].write(record_r1[2] + '\n')
                        file_names[key][0].write(record_r1[3] + '\n')
                        file_names[key][1].write(record_r2[0] + ' ' + str(index_pair) + '\n')
                        file_names[key][1].write(record_r2[1] + '\n')
                        file_names[key][1].write(record_r2[2] + '\n')
                        file_names[key][1].write(record_r2[3] + '\n')
                num_matched_reads += 1
            #checks if indexes are hopped
            elif record_i1[1] in index_dict.keys() and record_i2[1] in rev_comp_dict:
                #adds count to permutations dictionary
                perm_dict[index_pair] += 1
                hop_r1.write(record_r1[0] + ' ' + str(index_pair) + '\n')
                hop_r1.write(record_r1[1] + '\n')
                hop_r1.write(record_r1[2] + '\n')
                hop_r1.write(record_r1[3] + '\n')
                hop_r2.write(record_r2[0] + ' ' + str(index_pair) + '\n')
                hop_r2.write(record_r2[1] + '\n')
                hop_r2.write(record_r2[2] + '\n')
                hop_r2.write(record_r2[3] + '\n')
                num_hopped_reads += 1
            #indexes here are of high quality but unknown
            else:
                bad_egg_r1.write(record_r1[0] + ' ' + str(index_pair) + '\n')
                bad_egg_r1.write(record_r1[1] + '\n')
                bad_egg_r1.write(record_r1[2] + '\n')
                bad_egg_r1.write(record_r1[3] + '\n')
                bad_egg_r2.write(record_r2[0] + ' ' + str(index_pair) + '\n')
                bad_egg_r2.write(record_r2[1] + '\n')
                bad_egg_r2.write(record_r2[2] + '\n')
                bad_egg_r2.write(record_r2[3] + '\n')
                num_unknown += 1
        #indexes here do not pass quality test but could also be unknown
        else:
            bad_egg_r1.write(record_r1[0] + ' ' + str(index_pair) + '\n')
            bad_egg_r1.write(record_r1[1] + '\n')
            bad_egg_r1.write(record_r1[2] + '\n')
            bad_egg_r1.write(record_r1[3] + '\n')
            bad_egg_r2.write(record_r2[0] + ' ' + str(index_pair) + '\n')
            bad_egg_r2.write(record_r2[1] + '\n')
            bad_egg_r2.write(record_r2[2] + '\n')
            bad_egg_r2.write(record_r2[3] + '\n')
            #checks if indexes are unknown or bad quality for stats purposes
            if record_i1[1] in index_dict.keys() and record_i2[1] in rev_comp_dict:
                num_bad_quality += 1
            else:
                num_unknown += 1
    r1.close()
    r2.close()
    i1.close()
    i2.close()
    return num_bad_quality, num_hopped_reads, num_matched_reads, num_unknown

num_bad_quality, num_hopped_reads, num_matched_reads, num_unknown = sort_indexes(index_dict, rev_comp_dict, perm_dict, read1, read2, index1, index2)

#calculate stats and output to files 
percent_match = (num_matched_reads/number_records) * 100
percent_hopped = (num_hopped_reads/number_records) * 100
percent_bad_qual = (num_bad_quality/number_records) * 100
percent_unknown = (num_unknown/number_records) * 100

#overall counts and percentages
stats_file = open('stats.txt', 'w') 
stats_file.write("Total number of reads: " + str(number_records) + '\n')
stats_file.write("Number matched read pairs: " + str(num_matched_reads) + ', ' + str(percent_match) + '%' + '\n')
stats_file.write("Number hopped read pairs: " + str(num_hopped_reads) + ', ' + str(percent_hopped) + '%' + '\n')
stats_file.write("Number bad quality reads: " + str(num_bad_quality) + ', ' + str(percent_bad_qual) + '%' + '\n')
stats_file.write("Number unknown reads: " + str(num_unknown) + ', ' + str(percent_unknown) + '%' + '\n\n')
stats_file.write("Sample" + '\t' + 'Group' + '\t' + 'Treatment' + '\t' + 'Index' + '\t' + 'Percentage of reads' + '\n')

#percentage of reads per sample
for key in index_dict:
    for i in range(0, 4):
        stats_file.write(str(index_dict[key][i] + '\t'))
    index_match = str(key) + '-' + str(key)
    for key in perm_dict:
        if index_match == key:
            percent_total = (perm_dict[key]/number_records) * 100
            stats_file.write(str(percent_total))
    stats_file.write('\n')
stats_file.close()

#ouputs permutations to tab separated file
permutations_file = open('permutations.tsv', 'w')
permutations_file.write('Index pair' + '\t' + 'Occurence' + '\n')
for key, value in perm_dict.items():
    permutations_file.write(str(key) + "\t" + str(value) + '\n') 
permutations_file.close()


#close remaining files
for key in file_names:
    for i in range(0, 2):
        file_names[key][i].close()

hop_r1.close()
hop_r2.close()
bad_egg_r1.close()
bad_egg_r2.close()
