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
import sys
from population_functions import *



def make_all_dataframes(typing_report_name, imported_fmp_name, past_log_name):
    
    # locating files
    log_relative_path = "import_files" + "/" + past_log_name
    past_typing_relative_path = "import_files" + "/" + imported_fmp_name
    typing_report_relative_path = "import_files" + "/" + typing_report_name

    #making pandas dataframe of typing file
    typing_df = pd.read_excel(typing_report_relative_path, index_col=None) 
    typing_df['SampleID'] = typing_df['SampleID'].astype(str)

    # making a dataframe of only the end of the sheet, to add back on at the end
    start_of_run_metrics_position = -26 # should be the same everytime
    run_metrics_df = typing_df[start_of_run_metrics_position:]
    run_metrics_df = run_metrics_df.fillna('')

    # this will delete any rows without info for Barcode or Num reads - ie blank rows and all the info at the bottom of the sheet
    typing_df = typing_df[:-28]
    typing_df.dropna(subset=['Barcode', 'NumReads'], inplace=True)

    #making pandas DF of log of past uploads. This has barcode, cDNA, gDNA suggested ANRI code and reasons (should change suggested to submitted anri code?)
    past_log = pd.read_csv(log_relative_path)
    past_log['SampleID'] = past_log['SampleID'].astype(str)


    # making pandas dataframe of imported_fmp_types file
    fmp_raw = pd.read_csv(past_typing_relative_path)
    fmp_raw['sampleID'] = fmp_raw['sampleID'].astype(str)
    cohort_1_list = fmp_raw['sampleID'].to_list()
    fmp_raw = fmp_raw.set_index('sampleID')


    return typing_df, run_metrics_df, fmp_raw, past_log


def novels_in_upload(typing_df): 
        
        # for making a list later on 
        # then use this to check if a sample is confirmed by another sample in the typing file
        drop_not_novel = typing_df.dropna(subset=['MismatchDesc_cDNA', 'MismatchDesc_gDNA'], how='all')
        filling_na = drop_not_novel.fillna('empty')
        novs = filling_na.groupby(['Allele_cDNA','MismatchDesc_cDNA', 'Allele_gDNA', 'MismatchDesc_gDNA' ]).SampleID.agg(['count']).reset_index()
        
        # making a list
        novs = novs[novs['count'] > 1 ]
        novs_list = novs.values.tolist()

        return novs_list




def edit_imported_typing_file(dataframe, fmp_raw, past_log_raw, list_of_novels, row_counts):

        past_log_no_nan = past_log_raw.fillna('empty')

        # sets off all the new columns in one go
        dataframe['Automated_Quality_Checks'] = dataframe.apply(quality_check, axis=1)
        dataframe['none_available'] = dataframe.apply(none_available, axis=1)
        dataframe = previous_type_match(dataframe, fmp_raw)
        dataframe['perfect_match_to_reference'] = dataframe.apply(perfect_match_reference, axis=1)
        dataframe = already_accepted(dataframe, past_log_raw)
        dataframe['novel_type'] = dataframe.apply(novel_type, axis=1)
        dataframe = alleles_per_sample_count(dataframe, row_counts)
        dataframe['enough_alleles_per_sample'] = dataframe.apply(too_many_alleles, axis=1)
        dataframe = novel_confirmation(dataframe, past_log_no_nan)

        dataframe['internal_novel_confirmation'] = dataframe.apply(internal_novel_confirmation, axis=1)
        dataframe[['ANRI_code','manual_check_importance', 'ANRI_comment']] = dataframe.apply(ANRI_code_and_comment ,axis=1)


        return dataframe


def save_file(edited_typing_df, run_metrics_df, typing_report_name):


    # concatenating the typing bit with the run metrics
    edited_typing_df = edited_typing_df.append(pd.Series(), ignore_index=True)
    edited_typing_df = edited_typing_df.append(pd.Series(), ignore_index=True)
    entire_analysis = pd.concat([edited_typing_df, run_metrics_df])


    ###### saving 1) the entire analysis and 2) a typing report with suggestions

    # columns we dont want in the typing report
    # we are just keeping the suggested ANRI, manual check importance and reasons why this is suggested
    columns_to_drop = ['Automated_Quality_Checks', 'none_available', 'match_previous_type', 'previous_type_zygosity', 'perfect_match_to_reference', 'already_accepted','novel_type','alleles_per_sample','enough_alleles_per_sample','novel_confirmation', 'internal_novel_confirmation']
    typing_report_edited = entire_analysis.drop(columns_to_drop, axis=1)

    # sorting out what we're naming the files
    no_end_name = typing_report_name[:-4]
    edited_typing_report_name = f"{no_end_name}_typing_report.csv"
    entire_analysis_name = f"{no_end_name}_entire_analysis.csv"

    # saving
    typing_report_edited.to_csv(edited_typing_report_name, index=False)
    entire_analysis.to_csv(entire_analysis_name, index=False)


