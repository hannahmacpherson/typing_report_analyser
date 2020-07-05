from utils import *


def set_off_typing_report(imported_fmp_types, log_of_past_uploads, typing_file, fastq_file):
    
    # making dataframes 
    typing_df, run_metrics_df, fmp_raw, past_log_raw = make_all_dataframes(typing_file, imported_fmp_types, log_of_past_uploads)

    # making row counts by number of alleles
    row_counts = typing_df.groupby('SampleID').Gene.agg(['count']).reset_index()

    # making a list of novels that may be confirmed within the typing report
    list_of_novels = novels_in_upload(typing_df)

    # populating the typing dataframe with your new analysis columns
    edited_typing_df = edit_imported_typing_file(typing_df, fmp_raw, past_log_raw, list_of_novels, row_counts)

    print(f"{typing_file} analysed")

    # save files
    save_file(edited_typing_df, run_metrics_df, typing_file)
    print(f"{typing_file} saved")