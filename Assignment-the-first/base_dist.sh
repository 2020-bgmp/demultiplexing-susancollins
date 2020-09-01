#!/usr/bin/env python

import gzip
import argparse
import math
import matplotlib
import matplotlib.pyplot as plt

def get_args():
    '''Create arguments to run script in command line'''
    parser = argparse.ArgumentParser(description="Makes histogram of average quality score per base location")
    parser.add_argument("-f", "--file", help="filename")
    parser.add_argument("-r", "--readlength", help="length of reads, needed to make array", type=int)
    parser.add_argument("-o", "--outputfile", help="name of outputfile")
    return parser.parse_args()

args = get_args()

input_file = args.file
read_length = args.readlength
output_file = args.outputfile
read_length = read_length+1



def init_list(array, value=0.0):
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list with 101 values of 0.0.'''
    for i in range(0, read_length):
        array.append(value)
    return array

mean_scores = []
mean_scores = init_list(mean_scores)

def convert_phred(letter):
    """Converts a single character into a phred score"""
    phred = ord(letter) - 33
    return phred

def populate_list(input_file):
    """Populates a list of quality scores totals for each nucleotide location in a FASTQ file"""
    line_counter = 0
    mean_scores = []
    mean_scores = init_list(mean_scores)
    with gzip.open(input_file, 'tr') as fh:
        for line in fh:
            line_counter += 1
            if line_counter%4 == 0:
                letter_counter = 0
                for letter in line.rstrip():
                    score = convert_phred(letter)
                    mean_scores[letter_counter] += score
                    letter_counter += 1
    
    return mean_scores, line_counter

mean_scores, NR = populate_list(input_file)

def total_to_mean(list_of_totals, number_lines):
    """Takes a list of total quality score at each location and finds the mean"""
    real_mean_scores = []
    for index, value in enumerate(list_of_totals):
        real_mean_scores.append(list_of_totals[index]/(number_lines/4))
    return real_mean_scores

mean_scores = total_to_mean(mean_scores, NR)

index_list=[]

for index, value in enumerate(mean_scores):
    index_list.append(index)

plt.bar(index_list, mean_scores)
plt.xlim(left=-1, right=(read_length-1))
plt.ylabel('Mean Quality Score')
plt.xlabel('# Base')
plt.title('Average Quality Scores along Base Locations')
plt.savefig(output_file)
