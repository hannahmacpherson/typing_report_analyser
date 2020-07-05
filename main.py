
import sys
import openpyxl
from openpyxl import workbook
import json
import numpy as np
import pandas as pd
import xlrd
import csv
import os
from os import listdir
from os.path import isfile, join
from run_each_typing_report import *


# locating the folder of input typing reports
dir = os.getcwd()
folder = "import_files"
local_files_path = dir + "/" + folder

# making files in the folder into a list
input_files_list_raw = listdir(local_files_path)
typing_files_list_checked = []
fastq_files_list_checked = []

# only choosing typing files
for file in input_files_list_raw:
    if (file[-4:] == '.xls' or file[-5:] == '.xlsx') and 'TGSDEV' in file:
        typing_files_list_checked.append(file)

# only choosing FASTQ files
for file in input_files_list_raw:
    if file[-6:] == '.fastq' and 'TGSDEV' in file:
        fastq_files_list_checked.append(file)

number_of_files = len(typing_files_list_checked)

# giving the user an update 
print(f"Attempting to analyse {number_of_files} typing reports:")
print("")
for file in typing_files_list_checked:
    print(file)
print("")
print('###################################')
print("")

# matching FASTQ files to typing reports (based on the first part of the file name)
list_of_matched_files = []
list_of_typing_files_matched = []

counter = 0
for typing_file in typing_files_list_checked:
    for fastq_file in fastq_files_list_checked:
        if typing_file[:22] == fastq_file[:22]:
            list_of_matched_files.append([typing_file, fastq_file])
            list_of_typing_files_matched.append(typing_file)
            counter += 1

# what files haven't been matched?
non_matched_typing_reports = [file for file in typing_files_list_checked if file not in list_of_typing_files_matched]

# if every file has a match, counter should equal the number of typing files you have
# giving the user an update
if counter == number_of_files:
    print("Every typing report matched with its corresponding FASTQ file")
else: 
    print("Some typing reports have not been matched to FASTQ files:")
    print("")
    for file in non_matched_typing_reports:
        print(file)

# importing log and previous typing files
imported_fmp_types = 'previous_typing.csv'
log_of_past_uploads = 'log.csv'



for typing_file, fastq_file in list_of_matched_files:
    set_off_typing_report(imported_fmp_types, log_of_past_uploads, typing_file, fastq_file)



