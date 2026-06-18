import duckdb
import polars as pl

from datetime import datetime

from peptide_display.functions import load_config

start_time = datetime.now()

local_config = load_config("local_config.json")

screen_name = 'cmv1__hla_a_01_01__yeast_display__9mer__001'

explore_index_file_path = f"{local_config['output_root']}/combined__explore_index__brotli.parquet"
SQL_QUERY = f"select peptide from parquet_scan('{explore_index_file_path}') WHERE screen_slug = '{screen_name}' and p4!='H' and p5!='D' "

explore_results = duckdb.query(SQL_QUERY).pl()

print (f"Result shape: {explore_results.shape}")
print (f"Result columns: {explore_results.columns}")
print (explore_results)


filtered_screen_file_path = f"{local_config['output_root']}/private/yeast_display/{screen_name}__filtered__brotli.parquet"
SQL_QUERY = f"select * from parquet_scan('{filtered_screen_file_path}')"

screen_results = duckdb.query(SQL_QUERY).pl()


print (f"Result shape: {screen_results.shape}")
print (f"Result columns: {screen_results.columns}")
print (screen_results)

# join the two results on the peptide column
joined_results = explore_results.join(screen_results, on='peptide', how='inner')

print (f"Joined result shape: {joined_results.shape}")
print (f"Joined result columns: {joined_results.columns}")
print (joined_results)


elapsed_time = datetime.now() - start_time
print(f"Elapsed time: {elapsed_time}")

