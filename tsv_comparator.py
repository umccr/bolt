import argparse
import pandas as pd
import logging
import os

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

def load_tsv(file_path):
    """Load a TSV file into a pandas DataFrame."""
    df = pd.read_csv(file_path, delimiter='\t', low_memory=False, dtype=str)
    # Strip whitespace from all string columns to avoid issues with trailing spaces
    for col in df.select_dtypes(include=['object']).columns:
        df[col] = df[col].str.strip()
    return df

def ensure_same_dtype(df1, df2):
    """Ensure that both DataFrames have the same data types for all columns."""
    for column in df1.columns:
        if column in df2.columns:
            df2[column] = df2[column].astype(df1[column].dtype)

def compare_annotation_order(value1, value2, delimiter=','):
    """Compare annotation order for fields with a given delimiter."""
    if pd.isna(value1) or pd.isna(value2):
        return value1 == value2
    list1 = value1.split(delimiter)
    list2 = value2.split(delimiter)
    return sorted(list1) == sorted(list2)

def split_and_sort_vep_all_csq(value):
    """Split the VEP_ALL_CSQ value on ',' and then sort each individual annotation on '&' for comparison."""
    if pd.isna(value):
        return []
    return sorted(["&".join(sorted(annotation.split('&'))) for annotation in value.split(',')])

def compare_dataframes(df1, df2, log_file_path):
    """Compare two DataFrames and return differences record by record."""
    # Define key columns to uniquely identify records
    key_columns = ['CHROM', 'POS', 'REF', 'ALT']
    
    # Ensure both DataFrames have the same data types
    ensure_same_dtype(df1, df2)
    
    # Set the key columns as index to facilitate direct comparison
    df1.set_index(key_columns, inplace=True)
    df2.set_index(key_columns, inplace=True)
    
    # Find keys that are only in df1 or df2
    only_in_df1_keys = df1.index.difference(df2.index)
    only_in_df2_keys = df2.index.difference(df1.index)
    common_keys = df1.index.intersection(df2.index)

    # Get records only in df1
    only_in_df1 = df1.loc[only_in_df1_keys].reset_index()
    # Get records only in df2
    only_in_df2 = df2.loc[only_in_df2_keys].reset_index()

    # Prepare logging to a file
    with open(log_file_path, 'w') as log_file:
        # Compare records with common keys
        differences = []
        order_differences = []
        for key in common_keys:
            row1 = df1.loc[key]
            row2 = df2.loc[key]
            different = False
            order_difference = False
            for col in df1.columns:
                if pd.isna(row1[col]) and pd.isna(row2[col]):
                    continue  # Skip comparison if both values are NaN
                if row1[col] != row2[col]:
                    if col in ['CONSEQUENCE', 'VEP_ALL_CSQ']:
                        if col == 'VEP_ALL_CSQ':
                            sorted_row1 = split_and_sort_vep_all_csq(row1[col])
                            sorted_row2 = split_and_sort_vep_all_csq(row2[col])
                            if sorted_row1 != sorted_row2:
                                log_file.write(f"CHROM: {key[0]}, POS: {key[1]}, REF: {key[2]}, ALT: {key[3]}\n")
                                log_file.write(f"Order annotation difference in column '{col}': file1='{row1[col]}' vs file2='{row2[col]}'\n")
                                order_difference = True
                        else:
                            if not compare_annotation_order(row1[col], row2[col]):
                                log_file.write(f"CHROM: {key[0]}, POS: {key[1]}, REF: {key[2]}, ALT: {key[3]}\n")
                                log_file.write(f"Order annotation difference in column '{col}': file1='{row1[col]}' vs file2='{row2[col]}'\n")
                                order_difference = True
                    else:
                        log_file.write(f"CHROM: {key[0]}, POS: {key[1]}, REF: {key[2]}, ALT: {key[3]}\n")
                        log_file.write(f"Annotation difference in column '{col}': file1='{row1[col]}' vs file2='{row2[col]}'\n")
                        different = True
                        break
            if different:
                differences.append(key)  # Add key columns back to the list of differences
            elif order_difference:
                order_differences.append(key)  # Add key columns back to the list of order annotation differences

    # Create DataFrame from differences
    if differences:
        differences_df = pd.DataFrame(differences, columns=key_columns).reset_index(drop=True)
    else:
        differences_df = pd.DataFrame(columns=key_columns)

    # Create DataFrame from order annotation differences
    if order_differences:
        order_differences_df = pd.DataFrame(order_differences, columns=key_columns).reset_index(drop=True)
    else:
        order_differences_df = pd.DataFrame(columns=key_columns)
    
    # Return the key columns, records unique to each DataFrame, and records with differences
    return key_columns, only_in_df1, only_in_df2, differences_df, order_differences_df

def main():
    parser = argparse.ArgumentParser(description="Compare two TSV files record by record.")
    parser.add_argument("tsv_file1", help="Path to the first TSV file (original)")
    parser.add_argument("tsv_file2", help="Path to the second TSV file (new version)")
    args = parser.parse_args()

    tsv_file1 = args.tsv_file1
    tsv_file2 = args.tsv_file2

    # Derive log file name from input files
    log_file_name = f"{os.path.splitext(os.path.basename(tsv_file1))[0]}_{os.path.splitext(os.path.basename(tsv_file2))[0]}_differences.log"
    log_file_path = os.path.join(os.getcwd(), log_file_name)

    df1 = load_tsv(tsv_file1)
    df2 = load_tsv(tsv_file2)

    key_columns, only_in_df1, only_in_df2, differences, order_differences = compare_dataframes(df1, df2, log_file_path)

    logger.info(f"Number of records in {tsv_file1} (original): {len(df1)}")
    logger.info(f"Number of records in {tsv_file2} (new version): {len(df2)}")
    logger.info(f"Number of common records: {len(df1) - len(only_in_df1)}")
    logger.info(f"Number of records only in {tsv_file1}: {len(only_in_df1)}")
    logger.info(f"Number of records only in {tsv_file2}: {len(only_in_df2)}")
    logger.info(f"Number of common records with annotation differences: {len(differences)}")
    logger.info(f"Number of common records with order annotation differences: {len(order_differences)}")
    logger.info(f"Detailed differences are logged in: {log_file_path}")

    # Optionally, save the differences to a file
    if not differences.empty:
        differences.to_csv('differences.tsv', sep='\t', index=False)
    if not order_differences.empty:
        order_differences.to_csv('order_differences.tsv', sep='\t', index=False)

if __name__ == "__main__":
    main()
