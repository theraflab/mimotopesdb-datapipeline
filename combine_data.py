from fastparquet import write
import polars

# TODO this is ugly, refactor to a central config loader that can be used by both pipelines and the combined pipeline runner. For now, just load the config here to avoid circular imports.
from peptide_display.functions import load_config


def combine_data():
    """
    Run the data processing pipelines for both peptide display and CPL data.
    """

    local_config = load_config("local_config.json")

    print ("Combining the processed data from both pipelines into a single output file...")

    parquet_files = [
        'yeast_display__with_stats__brotli.parquet'
    ]

    screen_array = []

    for parquet_file in parquet_files:
        parquet_path = f"{local_config['output_root']}/{parquet_file}"
        print (f"Processed data file: {parquet_path}")
        # read the parquet file into a polars DataFrame to ensure it exists and is readable
        screen_type_df = polars.read_parquet(parquet_path)
        screen_type_df = screen_type_df.with_columns([
            polars.col('peptide_length').cast(polars.Int64)
        ])
        screen_array.append(screen_type_df)

    combined_df = polars.concat(screen_array, rechunk=True)
    combined_file_path = f"{local_config['output_root']}/screens__brotli.parquet"
    print (f"Writing combined screens data to: {combined_file_path} \n")
    combined_df.write_parquet(combined_file_path, compression='brotli')


if __name__ == "__main__":
    combine_data()