import polars
import csv
from fastparquet import write

from functions import load_config, kmerize_sequence, load_screen_slugs


def generate_explore_index(explore_data: polars.DataFrame, screen_slug: str) -> polars.DataFrame:
    """
    Generate an explore index for the yeast display data.
    The explore index is a parquet file containing the the amino acids at each position and the peptide sequencefor each peptide in the dataset.
    
    Args:
        explore_data (polars.DataFrame): The overview data for the yeast display screen.
    """

    peptides = explore_data['peptide'].to_list()

    explore_index = []
    
    for peptide in peptides:
        explore_row = {"peptide": peptide, "screen_slug": screen_slug}
        residues = kmerize_sequence(peptide, 1)
        for i, aa in enumerate(residues):
            explore_row[f"p{i + 1}"] = aa
        explore_index.append(explore_row)

    explore_index_df = polars.DataFrame(explore_index)
    return explore_index_df


local_config = load_config("local_config.json")
config = load_config("config.json")


availability = 'private'

screen_slugs = list(load_screen_slugs(local_config, 'yeast_display', availability).values())

explore_indices = []

for screen_slug in screen_slugs:
    filtered_file_path = f"{local_config['output_root']}/{availability}/yeast_display/{screen_slug}__filtered__brotli.parquet"

    filtered_data = polars.read_parquet(filtered_file_path)
    explore_index_df = generate_explore_index(filtered_data, screen_slug)

    explore_index_file_path = f"{local_config['output_root']}/{availability}/yeast_display/{screen_slug}__explore_index__brotli.parquet"
    print (f"Writing explore index data to: {explore_index_file_path} \n")
    explore_index_df.write_parquet(explore_index_file_path, compression='brotli')

    explore_indices.append(explore_index_df)

combined_explore_index = polars.concat(explore_indices, how='diagonal')
combined_explore_index_file_path = f"{local_config['output_root']}/combined__explore_index__brotli.parquet"
print (f"Writing combined explore index data to: {combined_explore_index_file_path} \n")
combined_explore_index.write_parquet(combined_explore_index_file_path, compression='brotli')