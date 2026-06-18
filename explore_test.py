import duckdb
import polars as pl

from datetime import datetime

from peptide_display.functions import load_config

start_time = datetime.now()

local_config = load_config("local_config.json")

screen_name = 'cmv1__hla_a_01_01__yeast_display__9mer__001'

combined_overview_file_path = f"{local_config['output_root']}/combined__explore_index__brotli.parquet"
SQL_QUERY = f"select peptide from parquet_scan('{combined_overview_file_path}') WHERE screen_slug = '{screen_name}' and p4!='H' and p5!='D' "

results = duckdb.query(SQL_QUERY).pl()

print (f"Result shape: {results.shape}")
print (f"Result columns: {results.columns}")
print (results)

elapsed_time = datetime.now() - start_time
print(f"Elapsed time: {elapsed_time}")

