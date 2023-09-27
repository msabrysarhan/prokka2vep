#!/usr/bin/env python3
"""
Author: Mohamed S. Sarhan
Email: mohamed.sarhan@eurac.edu; m.sabrysarhan@gmail.com
Date: September 27, 2023
Description: This Python script convert the GFF3 file of prokka to fit it to the format needed by VEP annotation.
"""

import argparse
import pandas as pd
import csv
import os

def process_gff_file(input_file, output_file):
    """
    This function takes the gff3 file that comes from prokka
    annotation and remove the gff headers and the fasta sequences 
    at the end of the file.

    Returns:
        a temp 9-column temp file wihtout headers nor fasta sequences 
    """
    in_fasta_section = False

    with open(input_file, 'r') as input_f, open(output_file, 'w') as output_f:
        for line in input_f:
            line = line.strip()

            if line.startswith('>'):
                in_fasta_section = True
                continue

            if not in_fasta_section and not line.startswith('#'):
                output_f.write(line + '\n')


def read_gff_as_dataframe(gff_file):
    """
    This function reads the temp gff file and returns pandas dataframe
    """
    print("Reading GFF file as pandas dataframe ...\n")
    gff_data = []

    with open(gff_file, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            fields = line.split('\t')

            # Ensure there are at least 9 fields (standard GFF format)
            if len(fields) < 9:
                continue

            if fields[0].startswith('#'):
                continue  # Skip comment lines

            seqname = fields[0]
            source = fields[1]
            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            score = fields[5]
            strand = fields[6]
            phase = fields[7]
            attributes = fields[8]

            # Parse attributes as a dictionary
            attribute_dict = {}
            for attr in attributes.split(';'):
                key, value = attr.split('=')
                attribute_dict[key] = value

            feature_data = [
                seqname,
                source,
                feature_type,
                start,
                end,
                score,
                strand,
                phase,
                attribute_dict
            ]

            gff_data.append(feature_data)

    columns = [
        'SeqName',
        'Source',
        'FeatureType',
        'Start',
        'End',
        'Score',
        'Strand',
        'Phase',
        'Attributes'
    ]

    gff_df = pd.DataFrame(gff_data, columns=columns)
    return gff_df


def create_transcript_df(df):
    """
    This function go through the pandas df and copy each gene line and modify it to 
    work as be used as a transcript line.

    Returns: transcripts only pandas dataframe
    """
    print("Creating transcript records ...\n")
    modified_rows = []
    # Iterate through the DataFrame
    for index, row in df.iterrows():
        if row['FeatureType'] == 'gene':
            # Create a copy of the row
            new_row = row.copy()
            # Change the "FeatureType" in the copied row
            new_row['FeatureType'] = 'transcript'
            # Append the modified row to the list
            modified_rows.append(new_row)

    transcript_df = pd.DataFrame(modified_rows)

    # Now, let's modify the attributes of this df
    for index, row in transcript_df.iterrows():
        transcript_df.at[index, 'Attributes']['ID'] = transcript_df.at[index, 'Attributes']['ID'].replace('_gene', '_transcript')
        transcript_df.at[index, 'Attributes']['Parent'] = transcript_df.at[index, 'Attributes']['ID'].replace('_transcript', '_gene')
        transcript_df.at[index, 'Attributes']['biotype'] = 'protein_coding'
    
    return transcript_df


def modify_df(df):
    """
    This function takes the original gff df and change the attributes of the parents and the IDs
    Returns: modified pandas dataframe
    """
    # Now, let's modify the attributes of this df
    for index, row in df.iterrows():
        if row['FeatureType'] == 'mRNA':
            df.at[index, 'Attributes']['Parent'] = df.at[index, 'Attributes']['ID'].replace('_mRNA', '_transcript')
            df.at[index, 'Attributes']['ID'] = df.at[index, 'Attributes']['ID'].replace('_mRNA', '_exon')
            df.at[index, 'FeatureType'] = 'exon'

        elif row['FeatureType'] == 'CDS':
            df.at[index, 'Attributes']['Parent'] = df.at[index, 'Attributes']['ID']+'_transcript'
            df.at[index, 'Attributes']['ID'] = df.at[index, 'Attributes']['ID'] + '_cds'
        
    return df


def merge_gffs(df1, df2):
    """
    This function merges the transcript df and the modified df

    Returns: merged unsorted dataframe
    """
    print("Merging GFF dataframes ... \n")
    merged_df = pd.concat([df1, df2], ignore_index=True)
    return merged_df


def reorder_gff(df):
    """
    This function reorder the rows of the merged gff df
    """
    print("Reordering GFF rows ... \n")
    df['SeqName'] = pd.to_numeric(df['SeqName'])
    df['Start'] = pd.to_numeric(df['Start'])
    df['End'] = pd.to_numeric(df['End'])


    # Group rows by the 'Start' and 'End' columns
    grouped = df.groupby(['Start', 'End'])
    
    reordered_rows = []
    
    for _, group in grouped:
        # Sort rows within each group based on FeatureType
        group = group.sort_values(by='FeatureType', key=lambda x: x.map({'gene': 0, 'transcript': 1, 'exon': 2, 'CDS': 3}))
        
        # Append the sorted rows to the reordered list
        reordered_rows.extend(group.to_dict(orient='records'))
    
    # Create a new DataFrame with the reordered rows
    reordered_df = pd.DataFrame(reordered_rows)
    
    # Reorder the DataFrame by 'SeqName', 'Start', and 'End'
    reordered_df = reordered_df.sort_values(['SeqName', 'Start', 'End'], ascending=[True, True, True])
    
    return reordered_df


def non_coding_rna(df):
    """
    This function removes the transcript lines related to the non-coding RNA 
    genes and adds biotype as t/rRNA.

    Returns: modified gff dataframe
    """
    print("Processing the non-coding RNA ... \n")
    drop_list = []
    for index, row in df.iterrows():

        if row['FeatureType'] == 'tRNA':
            df.at[index, 'Attributes']['biotype'] = 'tRNA'
            df.at[index, 'Attributes']['Parent'] = df.at[index, 'Attributes']['ID']+'_gene'
            df.at[index-1, 'Attributes']['Parent'] = df.at[index, 'Attributes']['ID']+'_gene'
            drop_list.append(index-2)

        elif row['FeatureType'] == 'rRNA':
            df.at[index, 'Attributes']['biotype'] = 'rRNA'
            df.at[index, 'Attributes']['Parent'] = df.at[index, 'Attributes']['ID']+'_gene'
            df.at[index-1, 'Attributes']['Parent'] = df.at[index, 'Attributes']['ID']+'_gene'
            drop_list.append(index-1)
                
    return df.drop(drop_list)


def write_gff_to_file(gff_dataframe, output_file_path):
    """
    Write a GFF pandas DataFrame to a GFF3 file with formatted attributes.

    Parameters:
        gff_dataframe (pd.DataFrame): The GFF DataFrame to be written.
        output_file_path (str): The path to the output GFF3 file.

    Returns:
        None
    """
    print(f"Writing GFF file to {output_file_path} \n")
    # Create a copy of the DataFrame to avoid modifying the original
    formatted_df = gff_dataframe.copy()

    # Format the 'Attributes' column with key=value pairs separated by semicolons
    formatted_df['Attributes'] = formatted_df['Attributes'].apply(
        lambda x: ';'.join([f"{key}={value}" for key, value in x.items()])
    )

    # Write the formatted DataFrame to a GFF3 file
    formatted_df.to_csv(output_file_path, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)
    print("GFF conversion is done. Thanks for using prokka2vep!")


def main():
    parser = argparse.ArgumentParser(description='Convert prokka GFF to VEP-friendly GFF')
    parser.add_argument('--gff', required=True, help='Input GFF file path')
    parser.add_argument('--out', required=True, help='Output file path')

    args = parser.parse_args()

    input_gff = args.gff
    output_gff = args.out

    process_gff_file(input_gff, output_gff+'.tmp')

    gff_df = read_gff_as_dataframe(output_gff+'.tmp')

    transcript_df = create_transcript_df(gff_df)
    
    modified_df = modify_df(gff_df)

    merged_df = merge_gffs(transcript_df, modified_df)

    reordered_df = reorder_gff(merged_df)

    non_coding_adjusted = non_coding_rna(reordered_df)

    write_gff_to_file(non_coding_adjusted, output_gff)

    os.remove(output_gff+'.tmp')

if __name__ == "__main__":
    main()
