
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



#1
def quality_check(row):
        if (row['NumReads'] < 40) or (row['PredictedAccuracy'] < 0.99) or (row['QV'] < 0.93):
            count_fails = 'FAILED'
        else:
            count_fails = 'PASS'
        return count_fails

#2 
def none_available(row):
    
    # is this row a pile of shit ie both references are 'none available'
    
    Allele_gDNA = row['Allele_gDNA']
    Allele_cDNA = row['Allele_cDNA']

    if Allele_gDNA == 'None available' and Allele_cDNA == 'None available':
        none_available = 'both_references_none_available'
    else:
        none_available = 'at_least_1_reference_available'

    return none_available

#3
def previous_type_match(dataframe, fmp_raw):

    def previous_type_match_inner(row):
        # does typing fit with previous known type? And was the old type homozygous?
        fit_previous_typing = 'unassigned'
        previous_zygosity = 'unassigned'

        sample = (row['SampleID'])
        gene = row['Gene']
        sample = (row['SampleID'])
        column_1_to_search = "hla_" + gene[4:].lower() + "_tgs_1"
        column_2_to_search = "hla_" + gene[4:].lower() + "_tgs_2"
        gDNA_full = (row['Allele_gDNA']).strip('()')
        gDNA_3_field = gDNA_full[:13]
        gDNA_2_field = gDNA_full[:10]

        cDNA_full = (row['Allele_cDNA']).strip('()')
        cDNA_3_field = cDNA_full[:13]
        cDNA_2_field = cDNA_full[:10]


        try:
            col_1_sample_in_fmp = fmp_raw.loc[sample, column_1_to_search]
            # get key error if sample name isn't in that previous type file
        except KeyError:
            fit_previous_typing = 'no_previous_data'
            previous_zygosity = 'no_previous_data'


        try:
            col_2_sample_in_fmp = fmp_raw.loc[sample, column_2_to_search]
        except KeyError:
            fit_previous_typing = 'no_previous_data'
            previous_zygosity = 'no_previous_data'

        # aka if we actually have previous data for this sample
        if (fit_previous_typing != 'no_previous_data') and (previous_zygosity != 'no_previous_data'):


            # checking for matching at various levels in the gDNA
            if gDNA_full != 'None available':
                if gDNA_full in col_1_sample_in_fmp or gDNA_full in col_2_sample_in_fmp:
                    fit_previous_typing = "full_length_match"
                elif gDNA_3_field in col_1_sample_in_fmp or gDNA_3_field in col_2_sample_in_fmp:
                    fit_previous_typing = "3_field_match"
                elif gDNA_2_field in col_1_sample_in_fmp or gDNA_3_field in col_2_sample_in_fmp:
                    fit_previous_typing = "2_field_match"
                    # NB for the below bit, it's already checked that the allele isn't in either column, 
                    # so the bit below means there's no data for this allele specifically. There could be data for the other allele
                elif col_1_sample_in_fmp == "No Data" or col_2_sample_in_fmp == "No Data":
                    fit_previous_typing = "no_previous_data"
                else:
                    fit_previous_typing = "not_a_match"

            # checking for matching at various levels in the cDNA if 'None Available' for gDNA (ie may be an extension)
            if gDNA_full == 'None available':
                if cDNA_full in col_1_sample_in_fmp or cDNA_full in col_2_sample_in_fmp:
                    fit_previous_typing = "full_length_match"
                elif cDNA_3_field in col_1_sample_in_fmp or cDNA_3_field in col_2_sample_in_fmp:
                    fit_previous_typing = "3_field_match"
                elif cDNA_2_field in col_1_sample_in_fmp or cDNA_3_field in col_2_sample_in_fmp:
                    fit_previous_typing = "2_field_match"
                else:
                    fit_previous_typing = "not_a_match" 


            # checking for previous zygosity
            if fit_previous_typing == "no_previous_data":
                previous_zygosity = "no_previous_data"
            elif col_1_sample_in_fmp == col_2_sample_in_fmp:
                previous_zygosity = "homozygote"
            elif col_1_sample_in_fmp != col_2_sample_in_fmp:
                previous_zygosity = "heterozygote"
            else:
                previous_zygosity = "error"



                    # returns two values, and when it's called later on it creates two separate columns per row
        return  pd.Series([fit_previous_typing, previous_zygosity])

    dataframe[['match_previous_type', 'previous_type_zygosity']] = dataframe.apply(previous_type_match_inner ,axis=1)

    return dataframe

#4
def perfect_match_reference(row):
    
    # does it match its reference sequence perfectlY ~~
    Mismatches_gDNA = row['Mismatches_gDNA']
    Mismatches_cDNA = row['Mismatches_cDNA']
    Allele_gDNA = row['Allele_gDNA']
    Allele_cDNA = row['Allele_cDNA']
    Gaps_cDNA = row['Gaps_cDNA']
    Gaps_gDNA = row['Gaps_gDNA']

    if Gaps_cDNA == 0 and Gaps_gDNA == 0 and Mismatches_cDNA == 0 and Mismatches_gDNA == 0 and Allele_gDNA != 'None available' and Allele_cDNA != 'None available':
        perfect_match = 'perfect_match_to_reference'
    else:
        perfect_match = 'variation_from_reference'

    return perfect_match
    
#5
def already_accepted(dataframe, past_log):
    def already_accepted_inner(row):
        sample = row['SampleID'] 
        # not stripped because if the sample has a novel and a reference allele but was previously typed homozygous,
        # you may have accepted the other allele which would be a perfect match but not this one with () which is a novel
        gDNA_full = (row['Allele_gDNA'])
        cDNA_full = row['Allele_cDNA']

        # making a mini dataframe for this sample in the log file
        df_sample_log = past_log[(past_log.SampleID) == sample]



        # if this sample has a gDNA reference, taking a subset of the above dataframe
        # with the same reference as this sample has
        if gDNA_full != 'None Available':
            # taking only instances of this sample where the gDNA is the same as here
            gDNA_sample_log = df_sample_log[df_sample_log.Allele_gDNA == gDNA_full]
            # making a list of ANRI codes submitted for these matching examples
            log_ANRI = gDNA_sample_log['ANRI_code'].to_list()

        else:
            # taking only instances of this sample where the cDNA is the same as here
            cDNA_sample_log = df_sample_log[df_sample_log.Allele_cDNA == cDNA_full]
            # making a list of ANRI codes submitted for these matching examples
            log_ANRI = cDNA_sample_log['ANRI_code'].to_list()



        # checking whether the list of codes for matching samples (by reference, not mismatches but they're implied
        # by ()s contains any 1s or 2s ie they were previously submitted.
        if '1' in log_ANRI or '2' in log_ANRI:
            already_accepted = True
        else: 
            already_accepted = False

        return already_accepted

    dataframe['already_accepted'] = dataframe.apply(already_accepted_inner, axis=1)
    return dataframe


#6
def novel_type(row):
    # used by other functions but basically if it's not an exact match, it suggests the kind of novel it may be. NB a lot of samples may just be typing errors but that's accounted for later on. having this separate just makes the later steps easier. making a mini dataframe for each sample in the log file

    gDNA_full = (row['Allele_gDNA']).strip('()')
    gDNA_mismatches = row['MismatchDesc_gDNA']
    cDNA_mismatches = row['MismatchDesc_cDNA']
    cDNA_nan = pd.isna(cDNA_mismatches)

    # gDNA is none available, cDNA has no mismatches. Probably an extension
    if (row.MismatchDesc_cDNA != 'None available') and (gDNA_full == 'None available') and (cDNA_nan == True):
        novel_type = 'extension'

    # gDNA is none available, cDNA has mismatches. cDNA novel
    elif pd.isna(cDNA_mismatches) == False and pd.isna(gDNA_mismatches) == False:
        novel_type = 'cDNA_novel'

    # gDNA has mismatches. gDNA novel
    elif pd.isna(gDNA_mismatches) == False and pd.isna(cDNA_mismatches) == True:
        novel_type = 'gDNA_novel'

    else:
        novel_type = 'N/A'

    return novel_type


#7 
def alleles_per_sample_count(dataframe, row_counts):
    def alleles_per_sample_count_inner(row):

     # the row_counts bit is made separately
    # groups whole typing file by sample ID and counts how many occurances there are

        sample = (row['SampleID'])

        # grouping by the sample ID and counting. Counts by Gene (originally counted by barcode 
        # but it screwed things up if different samples/genes were on the same bc)
        row_sample_ID = row_counts['SampleID']


        counter = row_counts.loc[row_sample_ID == sample, 'count'].iloc[0]
        return counter
    dataframe['alleles_per_sample'] = dataframe.apply(alleles_per_sample_count_inner, axis=1)
    return dataframe


#8
def too_many_alleles(row):
    # too many alleles yo
    # ie either > 2 or if it's previously typed as homozygous then > 1
    
    previous_type_zygosity = row['previous_type_zygosity']
    alleles_per_sample = row['alleles_per_sample']

    if previous_type_zygosity == 'homozygote':
        expected_alleles = 1
    elif previous_type_zygosity == 'heterozygote':
        expected_alleles = 2
    elif previous_type_zygosity == 'no_previous_data':
        expected_alleles = 'no_previous_data'
    else:
        print("error in too many alleles script part 1")

    if expected_alleles == 'no_previous_data':
        too_many_alleles = 'no_previous_data'
    elif alleles_per_sample == expected_alleles:
        too_many_alleles = 'correct_number_of_alleles'
    elif alleles_per_sample > expected_alleles:
        too_many_alleles = 'too_many_alleles'
    elif alleles_per_sample < expected_alleles:
        too_many_alleles = 'too_few_alleles'
    else: 
        print("error in too many alleles script part 2")

    if alleles_per_sample > 2:
        too_many_alleles = 'too_many_alleles'

    return too_many_alleles


#9
def novel_confirmation(dataframe, past_log_no_nan):

    def novel_confirmation_inner(row):
        # if the sample has a chance of being novel, this checks whether we have seen this combo before 
        # ie same gDNA/cDNA reference (might only be cDNA if eg extension) and same mismatch descriptions
        # then checks whether this occured on either a different sample or a different barcode and tells you
        # if this sample can confirm it or not 

        def isNaN(value):
            return value != value

        finished = False

        # values for sample
        novel_type = row['novel_type']
        Allele_cDNA_sample = row['Allele_cDNA']
        Allele_gDNA_sample = row['Allele_gDNA']
        MismatchDesc_gDNA_sample = row['MismatchDesc_gDNA']
        MismatchDesc_cDNA_sample = row['MismatchDesc_cDNA']
        sample = row['SampleID']
        barcode = row['Barcode']

        # to stop nan comparisons (nan doesnt equal itself)
        if isNaN(MismatchDesc_gDNA_sample) == True:
            MismatchDesc_gDNA_sample = 'empty'

        if isNaN(MismatchDesc_cDNA_sample) == True:
            MismatchDesc_cDNA_sample = 'empty'

        # to track if anything isn't analysed
        novel_confirmation = 'not_assigned'



        # if it's got no chance of being novel
        if novel_type == 'N/A':
            novel_confirmation = 'N/A'
            finished = True
        else:
            # if potentially novel, dataframe of matching gDNA, cDNA, gDNA MM and cDNA MMs
            df_matching_in_log = past_log_no_nan[((past_log_no_nan.Allele_cDNA) == Allele_cDNA_sample)]

            df_matching_in_log = df_matching_in_log[((df_matching_in_log.Allele_gDNA) == Allele_gDNA_sample)]
            df_matching_in_log = df_matching_in_log[((df_matching_in_log.MismatchDesc_gDNA) == MismatchDesc_gDNA_sample)]
            df_matching_in_log = df_matching_in_log[((df_matching_in_log.MismatchDesc_cDNA) == MismatchDesc_cDNA_sample)]

        if finished == False:

            # making a set of seen sample IDs and barcodes. NB sets are only unique values:
            unique_sightings = set()
            # add this sample
            unique_sightings.add((sample, barcode))
            # add other samples from log file
            for index, logrow in df_matching_in_log.iterrows():
                unique_sightings.add((logrow['SampleID'], logrow['Barcode']))


            unique_sightings_of_sequence = int(len(unique_sightings))


            # if it's only been seen once, it's in this typing report, so needs a repeat
            if unique_sightings_of_sequence == 1:
                novel_confirmation = f"{novel_type}_needs_1_repeat"

            # if it's two, then it's been seen twice on different samples and or barcodes
            if unique_sightings_of_sequence == 2:
                novel_confirmation = f"{novel_type}_confirmed_by_this_sample"

            # if more than two then it must have already been in the log twice
            if unique_sightings_of_sequence > 2:
                novel_confirmation = f"{novel_type}_previously_confirmed"

        return novel_confirmation

    dataframe['novel_confirmation'] = dataframe.apply(novel_confirmation_inner, axis=1)
    return dataframe


#10
def internal_novel_confirmation(row):
        
    # does this sequence (if potentially novel) confirm itself within the typing report?
    internal_novel_confirmation = 'unassigned'

    # variables which might be useful in working this out

    row = row.fillna('empty')
    Allele_cDNA = row['Allele_cDNA']
    Allele_gDNA = row['Allele_gDNA']
    MismatchDesc_cDNA = row['MismatchDesc_cDNA']
    MismatchDesc_gDNA = row['MismatchDesc_gDNA']
    novel_type = row['novel_type']

    # if it's got no chance of being novel
    try:
        if novel_type == 'N/A':
            internal_novel_confirmation = 'N/A'

        # using a list made previously which takes into account: 'Allele_cDNA','MismatchDesc_cDNA', 'Allele_gDNA', 'MismatchDesc_gDNA'
        # it counts all rows in this typing report that are potentially novel and matches on these columns
        # it then deletes anything that has only one count (ie not confirmed)
        # if this row matches something on this list, it has been confirmed by this typing report
        # warning is it doesn't check for QC so will include common naff typing
        # however the actual classification at the end of this script will check for QCs first
        # so at least one of the samples confirming each other will have a good QC - unlikely if just shit typing (except for RRAs)

        else:
            if confirmed_novel in novs_list:
                if (Allele_cDNA == confirmed_novel[0]) and (MismatchDesc_cDNA == confirmed_novel[1]) and (Allele_gDNA == confirmed_novel[2]) and (MismatchDesc_gDNA == confirmed_novel[3]):
                    internal_novel_confirmation = 'confirmed_internally'
                else:
                    internal_novel_confirmation = 'not_confirmed_internally'
            elif confirmed_novel not in novs_list:
                internal_novel_confirmation = 'not_confirmed_internally'
    except NameError:
        internal_novel_confirmation = 'not_confirmed_internally'

    return internal_novel_confirmation





def ANRI_code_and_comment(row):
    
    # these are here to assess if anything has slipped through the gaps - will still be these values at the end if so
    
    ANRI_score = 'unscored'
    checkpoint2 = 'not used'
    checkpoint3 = 'not used'
    checkpoint4 = 'not used'
    manual_check_importance = 'not used'
    
    # all the variables that may be useful in suggesting how you should analyse this row
    
    num_reads = row['NumReads']
    Automated_Quality_Checks = row['Automated_Quality_Checks']
    none_available = row['none_available']
    match_previous_type = row['match_previous_type']
    perfect_match_to_reference = row['perfect_match_to_reference']
    previous_type_zygosity = row['previous_type_zygosity']
    already_accepted = row['already_accepted']
    alleles_per_sample = row['alleles_per_sample']
    enough_alleles_per_sample = row['enough_alleles_per_sample']
    novel_confirmation = row['novel_confirmation']
    internal_novel_confirmation = row['internal_novel_confirmation']
    

    # should you class this row as 0 ie ignore. 
    #Yes, 1) if this allele has already been accepted (1 or 2) for this sample,
    # 2) if both cDNA and gDNA on the typing report are 'None Available'
    # 3) if it's not a match to the previous data you have for this allele, 
    # re number 3. if it also has too many alleles it will suggest this may be a poor quality version of an acceptable allele. 
    # NB if there is no previous data for this sample it will carry on through the chain
    
    
    def is_zero():
        
        zero = False
        reason = ""
        if already_accepted == 'TRUE':
            zero = True
            reason += "Already accepted this allele. "
        if none_available == 'both_references_none_available':
            zero = True
            reason += "Both cDNA and gDNA on the typing report are 'None Available'. "
        if match_previous_type == 'not_a_match':
            zero = True
            reason += "Not a match to the previous data you have for this allele. "
            if enough_alleles_per_sample == 'too_many_alleles':
                zero = True
                reason += "There are too many alleles for this sample, this is probably a low quality version of an acceptable one. "
                #check_manually = True
                
        return zero, reason

    
    # setting off the above function. Every row will pass through this initial test.
    zero, reason = is_zero()
    
    # checkpoint 1
    if zero == True:
        ANRI_score = 0
        reason = reason
        manual_check_importance = 'low'
        
    ########################### Every 0 should be dealt with now
    ########################### now we assess everything which has a perfect match to a reference 
    ########################### (1, 2 and some 3s mixed in)
        
    # checkpoint 2 - it passed checkpoint 1, but does it have a perfect match to a reference sequence?
    if ANRI_score == 'unscored':
        # stops 0s from going through further
        if perfect_match_to_reference == 'perfect_match_to_reference':
            checkpoint2 = True
        else: 
            checkpoint2 = False
        
    
    # checkpoint 3 - it has a perfect match to a reference sequence, but does it pass QC?
    
    if checkpoint2 == True:
        if Automated_Quality_Checks == 'PASS':
            checkpoint3 = True
        else: 
            checkpoint3 = False
    
    # if it has a perfect match to a reference sequence and passes QC, is it a 1 or a 2?
    def one_or_two():
        
        manual_check_importance = 'unassigned'
        reason = "Perfect match to reference sequence, passed QC. "
        
        # do we have previous type data for this allele?
        if match_previous_type == 'no_previous_data':
            previous_type = False
        else:
            previous_type = True
        
        # what was the previous zygosity of this sample in our previous data?
        if previous_type == True:
            length_of_match = f"Level of match to previous typing: {match_previous_type}. "
            if previous_type_zygosity == 'heterozygote':
                ANRI_score = 1
                reason += f" Previously typed as a heterozygote. {length_of_match}"
                manual_check_importance = 'low'
            if previous_type_zygosity == 'homozygote':
                if num_reads >= 100:
                    if alleles_per_sample == 2:
                        ANRI_score = "1/2"
                        read_description = f'. However another allele present. If the other allele is also a match to the expected typing, this sample is likely heterozygous at 4 fields. However, acceptable number of reads for a homozygote ({num_reads})'
                        manual_check_importance = 'high'
                    elif alleles_per_sample == 1:
                        ANRI_score = 2
                        read_description = f"Acceptable number of reads for a homozygote ({num_reads})"
                        manual_check_importance = 'low'
                    elif alleles_per_sample > 2:
                        ANRI_score = 1/2
                        read_description = f" More than 2 alleles present. Manually check. Acceptable number of reads for a homozygote ({num_reads})"
                        manual_check_importance = 'high'
                elif num_reads < 100:
                    ANRI_score = "2/3"
                    read_description = f"Only {num_reads}, which is low for a homozygote, might be worth rerunning. "
                    manual_check_importance = 'high'
                reason += f" Previously typed as a homozygote. {length_of_match}"
                reason += read_description
        
        
        # We don't have previous data for this allele, can we accept it based on number of alleles
        # in typing report and number of reads?
        if previous_type == False:
            reason += "No previous data for this allele, "
            
            # if there are 2 alleles in the typing report, we already know this one matches and has > 40 reads,
            # so can accept it as a 1 
            
            if alleles_per_sample == 2:
                ANRI_score = 1
                reason += " but there are 2 alleles for this sample on the typing report and it has passed quality control."
                manual_check_importance = 'medium'
            # if there are > 2 alleles in the typing report and we don't have previous typing data,
            # we need to look at this sample manually. It is probably either 1 or 3 depending on the other alleles
            
            elif alleles_per_sample > 2:
                ANRI_score = "1/3"
                reason += f" and there are 3 alleles for this sample on the typing report. Check manually."
                manual_check_importance = 'high'
            # if there is 1 allele in the typing report and reads are > 100, it's probably fair to assume it's homozygous
            # however if there's less than 100 reads, it's probably a 3 but manually check
            
            elif alleles_per_sample == 1:
                if num_reads >= 100:
                    ANRI_score = 2
                    reason = f" but acceptable number of reads for a homozygote ({num_reads})"
                    manual_check_importance = 'medium'
                elif num_reads < 100:
                    ANRI_score = 3
                    reason = f" and only {num_reads}, which is low for a homozygote. Rerun this sample"
                    manual_check_importance = 'high'

        return ANRI_score, reason, manual_check_importance
                

            
    # setting off the one_or_two function only for relevant rows
    if checkpoint2 == True and checkpoint3 ==True:
        ANRI_score, reason, manual_check_importance = one_or_two()
        
        
    # if a sample matches reference perfectly but hasn't passed QC we should check this manually, haven't yet
    # checked if homozygous or heterozygous though so have put three - should probs be 1/2/3
    if checkpoint2 == True and checkpoint3 == False:
        reason = "Check sample manually. Perfect match to reference sequence but failed quality checks."
        manual_check_importance = 'high'
        if previous_type_zygosity == 'homozygote':
            ANRI_score = '2/3'
        elif previous_type_zygosity == 'heterozygote':
            ANRI_score = '1/3'
        else:
            ANRI_score = '1/2/3'
        
    
    ########################### we have now assessed every row with a perfect match to a reference sequence or ANRI =0. 
    ########################### now we assess everything which is left (novel, sequencing is shit or RRAs)
    
    # checkpoint 2 was 'is it a perfect match to its reference' ie no mismatches at cDNA or gDNA level
    # does it pass QC (good indicator as to whether it's shit sequencing or a novel)
    if checkpoint2 == False:
        if Automated_Quality_Checks == 'PASS':
            checkpoint4 = True
        else: 
            checkpoint4 = False
            
    # if checkpoint 4 is false, it's failed QC, so will be given a 3. Even if it's actually novel, you should
    # run again anyway to make sure it's not just crap sequencing.
    if checkpoint4 == False:
        ANRI_score = 3
        reason = "Sample is not an exact match to the reference and failed QC. It may be novel, but more likely the sequencing is poor. If other alleles for this sample are acceptable and this one is similar to one of those, ignore this one. Ideally repeat."
        manual_check_importance = 'high'
    #### should consider elaborating on this. If you checked previous typing for homozygote/heterozygote
    #### and there are other acceptable alleles and this sample is similar to an acceptable allele, 
    #### you could change this to 0 rather than asking for a manual check (or at least manual check is less important)
    
    
    # if checkpoint4 is true, it's passed QC whilst not being a perfect match. This is probably a novel sample
    # should consider adding to this whether it matches at a low resolution - would bolster how much you can trust it
    # also would be worth adding in a way of checking for RRAs - checking against a list for known variants in particular alleles
    # for now though, putting a 5 should initiate a manual check (even just to see where the SNP or whatever is)
    if checkpoint4 == True:
        ANRI_score = 5
        reason = "Does not match reference perfectly but has passed QC. "
        manual_check_importance = 'high'
        
        log_novel_confirmation = 'Not assessed'
        internal_novel_confirmation = 'Not assessed'
        # has this potential novel been confirmed in this uploaded typing report?
        # NB if this is the case, whilst this sample has passed QC, we can't be sure that the other sample(s) confirming it have
        # however the likelihood is that it is genuine if this one has passed QC

        if (internal_novel_confirmation != 'not_confirmed_internally') and (internal_novel_confirmation != 'N/A'):
            internal_confirmation = True
            reason_adder1 = "The same sequence has been seen at least one other time in this typing report, and likely confirms your novel. Check whether it's the same sample on a different barcode or a different sample altogether. NB whilst this sequence has passed QC, the other hasn't necessarily so check that too."
            reason += reason_adder1
        else:
            internal_novel_confirmation = False # don't need to add any extra info to reason
        
        # has this potential novel been confirmed in anything we have previously run through the script?
        # NB this is currently based on the gDNA/cDNA reference sequence and the mismatches seen in past uploads. 
        # it also checks that you've seen the same thing on 2 or more barcodes
        # it doesn't however currently base it on whether you previously gave it a 5, or it passed QC or anything so i should 
        # improve this.
    
        if novel_confirmation != 'N/A':
            log_novel_confirmation = True
            reason_adder2 = f"After looking through your log of past script runs, the status of this potential novel is: {novel_confirmation}. This takes into account not only how many times this sequence has been seen, but whether it has been seen on multiple barcodes. Make sure you look at the mismatches before submitting though."
            reason += reason_adder2
            
        if (log_novel_confirmation == True) and (internal_confirmation == True):
            reason = "Does not match reference perfectly but has passed QC. Can be confirmed both by at least one other allele in this typing report (although check if it's the same sample on the same barcode) and your past log of uploads (which accounts for barcodes)."
    
    return pd.Series([ANRI_score, manual_check_importance, reason])
        

    