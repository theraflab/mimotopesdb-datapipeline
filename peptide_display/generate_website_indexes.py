import duckdb
import os
import polars
from fastparquet import write

from functions import load_config, kmerize_sequence, load_screen_slugs


def generate_search_index(overview_data: polars.DataFrame, screen_slug: str, depth: int = 20) -> None:
    """
    Generate a search index for the yeast display data.
    The search index is a parquet file containing the kmers and their positions for each peptide in the dataset.
    
    Args:
        overview_data (polars.DataFrame): The overview data for the yeast display screen.
        screen_slug (str): The slug of the yeast display screen.
        depth (int): The number of peptides to process for generating the search index.
    
    """

    peptides = overview_data['peptide'].to_list()

    kmer_positions = {}
    
    for peptide in peptides[:depth]:
        kmers = kmerize_sequence(peptide, 3)
        for i, kmer in enumerate(kmers):
            if i not in kmer_positions:
                kmer_positions[i] = set()
            kmer_positions[i].add(kmer)

    search_index = []
    for position, kmers in kmer_positions.items():
        for kmer in kmers:
            search_index.append({
                "position": position,
                "kmer": kmer,
                "screen_slug": screen_slug
            })
    return search_index, kmer_positions


local_config = load_config("local_config.json")
config = load_config("config.json")


availability = 'private'

screen_slugs = list(load_screen_slugs(local_config, 'yeast_display', availability).values())

search_indices = []
dataframes = []

depth = 20


for screen_slug in screen_slugs:
    overview_file_path = f"{local_config['output_root']}/{availability}/yeast_display/{screen_slug}__1c__brotli.parquet"

    overview_data = polars.read_parquet(overview_file_path)

    search_index, kmer_positions = generate_search_index(overview_data, screen_slug, depth)

    overview_data = overview_data.with_columns(polars.lit(screen_slug).alias('screen_slug'))

    dataframes.append(overview_data)
    search_indices.extend(search_index)


search_index_df = polars.DataFrame(search_indices)
search_index_file_path = f"{local_config['output_root']}/search_index__brotli.parquet"
print (f"Writing combined search index  data to: {search_index_file_path}")
search_index_df.write_parquet(search_index_file_path, compression='brotli')


combined_overview_data = polars.concat(dataframes, how='diagonal')

# set the column order to match the standardised columns, adding any missing columns with null values
columns = config['yeast_display']['columns']['standardised_columns'].copy()
for col in columns:
    if col not in combined_overview_data.columns:
        combined_overview_data = combined_overview_data.with_columns(polars.lit(None).alias(col))

combined_overview_file_path = f"{local_config['output_root']}/yeast_display_combined__1c__brotli.parquet"
print (f"Writing combined overview data to: {combined_overview_file_path}")
combined_overview_data.write_parquet(combined_overview_file_path, compression='brotli')
