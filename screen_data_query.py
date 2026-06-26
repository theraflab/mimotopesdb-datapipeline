import duckdb
import polars as pl

from datetime import datetime

from peptide_display.functions import load_config

start_time = datetime.now()

local_config = load_config("local_config.json")

screen_name = 'ebv1__hla_a_02_01__yeast_display__9mer__001'

combined_overview_file_path = f"{local_config['output_root']}/yeast_display_combined__1c__brotli.parquet"
SQL_QUERY = f"select rank, peptide, naive, round_1, round_2, round_3, round_4 from parquet_scan('{combined_overview_file_path}') WHERE screen_slug = '{screen_name}'"

results = duckdb.query(SQL_QUERY).pl()

print (f"Result shape: {results.shape}")
print (f"Result columns: {results.columns}")
print (results)

elapsed_time = datetime.now() - start_time
print(f"Elapsed time: {elapsed_time}")