from typing import List, Tuple, Dict
import csv
import json
import os
import polars
from fastparquet import write



from functions import load_config, map_column_name, lookup_screen_slug, is_valid_peptide







def load_data(csv_file_path: str) -> Tuple[List[str], List[List[str]]]:
    """
    Load data from a CSV file and return it as a list of dictionaries.
    
    Args:
        csv_file_path (str): The path to the CSV file to load.  
    """
    data = []
    with open(csv_file_path, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        data = [row for row in reader if any(cell.strip() for cell in row)]
    # return as a polars DataFrame for easier manipulation
    raw_data = polars.DataFrame(data, schema=header, orient='row')
    return header, raw_data


def standardise_data(raw_data: List, header:List, config: Dict, screen_type: str) -> List:
    """

    
    """
    column_names = [map_column_name(config, screen_type, col.lower()) for col in header]
    # first we'll select only the columns that we have mapped, and rename them to their standardised names
    standardised_data = raw_data.select([
        polars.col(orig).alias(mapped)
        for orig, mapped in zip(header, column_names)
        if mapped is not None
    ])
    # filter out any rows where the peptide is not valid
    standardised_data = standardised_data.with_columns(
        polars.col('peptide')
        .map_elements(is_valid_peptide, return_dtype=polars.Boolean)
        .alias('is_valid_peptide')
    )
    cleaned_data = standardised_data.filter(polars.col('is_valid_peptide')).drop('is_valid_peptide')
    # finally we'll cast the round columns to integers
    cleaned_data = cleaned_data.with_columns([
        polars.col(col).cast(polars.Float64).cast(polars.Int64)
        for col in cleaned_data.columns
        if col.startswith('round_')
    ])
    max_round = max([int(col.split('_')[-1]) for col in column_names if col and col.startswith('round_')])
    # sort the cleaned data by the counts in the highest round column
    cleaned_data = cleaned_data.sort(by=f'round_{max_round}', descending=True)

    # finally now we have the ordering correct, we can add a rank column to the cleaned data
    cleaned_data = cleaned_data.with_row_index(name='rank', offset=1)

    # now we're going to set a column which is the delta between the count of the highest round of the previous row and the current row, to give an idea of how much the counts are dropping off
    # the delta for the first row will be set to 0, as there is no previous row to compare to
    # we'll make the delta always be positive, as we are only interested in the magnitude of the drop off, not the direction
    cleaned_data = cleaned_data.with_columns([
        (polars.col(f'round_{max_round}') - polars.col(f'round_{max_round}').shift(1)).abs().fill_null(0).alias('delta')
    ])
    # if the delta is over half the count of the previous row, then we will set a flag to indicate that this is a significant drop off
    cleaned_data = cleaned_data.with_columns([
        (polars.col('delta') > (polars.col(f'round_{max_round}').shift(1) * 0.5)).fill_null(False).alias('significant_drop_off')
    ])
    return column_names, cleaned_data, max_round


def filter_data_by_round(data: polars.DataFrame, specified_round: int, cutoff:int = 1) -> polars.DataFrame:
    """
    Filter the data to include only rows above the cutoff for the specified round.
    
    Args:
        data (polars.DataFrame): The input data to filter.
        specified_round (int): The round number to filter by.
        cutoff (int): The minimum value for the specified round to include a row in the filtered data.
        
    Returns:
        polars.DataFrame: The filtered data containing only rows for the specified round.
    """
    round_column = f'round_{specified_round}'
    if round_column in data.columns:
        return data.filter(polars.col(round_column) > cutoff)
    else:
        raise ValueError(f"Round column '{round_column}' does not exist in the data.")


def load_and_process_screen(local_config: Dict, file_name: str, screen_slug: str, availability: str, config: Dict, screen_type: str) -> None:
    """
    Load and process data a specific screen from the given folder and file name. This function reads the data, standardizes it, filters it by the maximum round, and prints relevant information.
    
    Args:
        local_config (Dict): The local configuration dictionary.
        file_name (str): The name of the data file to process.
        screen_slug (str): The slug identifier for the screen.
        availability (str): The availability status of the screen.

    Returns:
        None
    """

    input_folder = f"{local_config['input_root']}/yeast_display/{availability}"
    output_folder = f"{local_config['output_root']}/{availability}/yeast_display"
    file_path = f"{input_folder}/{file_name}"
    column_names = []
    header, raw_data = load_data(file_path) 

    print (f"Header: {header}")
    print (f"Raw data shape: {raw_data.shape}")
    print (f"Raw data columns: {raw_data.columns}")
    print (f"Raw data: {raw_data.head()}")
    column_names, cleaned_data, max_round = standardise_data(raw_data, header, config, screen_type)
    print (f"Cleaned data shape: {cleaned_data.shape}")
    print (f"Cleaned data columns: {cleaned_data.columns}")
    print (f"Cleaned data: {cleaned_data.head()}")
    filtered_data = filter_data_by_round(cleaned_data, max_round)
    print (f"Filtered data shape: {filtered_data.shape}")
    print (f"Filtered data columns: {filtered_data.columns}")
    print (f"Filtered data: {filtered_data.head()}")


    stats = {
        'max_round': max_round,
        'raw_data_count': raw_data.shape[0],
        'cleaned_data_count': cleaned_data.shape[0],
        'filtered_data_count': filtered_data.shape[0],
        'columns': ','.join(filtered_data.columns),
        'column_count': len(filtered_data.columns)
    }

    cleaned_file_path = f"{output_folder}/{screen_slug}__cleaned__brotli.parquet"
    print (f"Writing cleaned data to: {cleaned_file_path}")
    cleaned_data.write_parquet(cleaned_file_path, compression='brotli')  


    filtered_file_path = f"{output_folder}/{screen_slug}__filtered__brotli.parquet"
    print (f"Writing filtered data to: {filtered_file_path}")
    filtered_data.write_parquet(filtered_file_path, compression='brotli')


    onek_file_path = f"{output_folder}/{screen_slug}__1k__brotli.parquet"
    print (f"Writing 1k data to: {onek_file_path}")
    onek_data = filtered_data[:1000]
    print (f"1k data shape: {onek_data.shape}")
    onek_data.write_parquet(onek_file_path, compression='brotli')


    overview_file_path = f"{output_folder}/{screen_slug}__1c__brotli.parquet"
    print (f"Writing overview data to: {overview_file_path}")
    overview_data = filtered_data[:100]
    print (f"Overview data shape: {overview_data.shape}")
    overview_data.write_parquet(overview_file_path, compression='brotli')

    print ('')

    return stats









def load_raw_data() -> None:
    """
    Load and process raw data for all screens in the specified folder and availability.
    
    Args:
        local_config (Dict): The local configuration dictionary.
        screen_type (str): The type of screen to process (e.g., 'yeast_display').
        availability (str): The availability status of the screens (e.g., 'private').
        
    Returns:
        None
    """
    local_config = load_config("local_config.json")
    config = load_config("config.json")
    screen_type = 'yeast_display'
    availability = 'private'

    screens_csv = f"{local_config['output_root']}/screens.csv"

    with open(screens_csv, 'r') as f:
        reader = csv.DictReader(f)
        screens = {row['screen_slug']: row for row in reader if len(row['screen_slug']) > 0 and row['screen_type'] == screen_type}

    folder_name = f"{local_config['input_root']}/{screen_type}/{availability}"

    files = sorted([f for f in os.listdir(folder_name) if f.endswith('.csv')])

    augmented_screens = {}
    slugs = []
    for filename in files:
        print (f"Reading file: {filename}")
        screen_slug = lookup_screen_slug(local_config, filename, screen_type)
        print (f"Screen slug: {screen_slug}")
        slugs.append(screen_slug)
        stats = load_and_process_screen(local_config, filename, screen_slug, availability, config, screen_type)
        print (f"Stats: {stats} \n\n")
        augmented_screens[screen_slug] = {**screens[screen_slug], **stats}
        print (f"Updated screens dict for '{screen_slug}']: {screens[screen_slug]}")

    screen_list = list(augmented_screens.values())
    output_file = f"{local_config['output_root']}/yeast_display__with_stats__brotli.parquet"
    screen_polars = polars.DataFrame(screen_list)
    print (screen_polars.head())
    print (f"Writing screens with stats to: {output_file}")
    screen_polars.write_parquet(output_file, compression='brotli')

if __name__ == "__main__":
    load_raw_data()
