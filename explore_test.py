import duckdb
import polars as pl

from datetime import datetime

from peptide_display.functions import load_config

start_time = datetime.now()

local_config = load_config("local_config.json")

screen_name = 'brenetafusp__hla_a_02_01__yeast_display__9mer__001'

filters = [{"p5": "!H"}, {"p7": "I"}]

explore_index_file_path = f"{local_config['output_root']}/combined__explore_index__brotli.parquet"
SQL_QUERY = f"select peptide from parquet_scan('{explore_index_file_path}') WHERE screen_slug = '{screen_name}' and "
for filter in filters:
    for position, aa in filter.items():
        if '!' in aa:
            aa = aa.replace('!', '')
            SQL_QUERY += f"{position} != '{aa}'"
        else:
            SQL_QUERY += f"{position} = '{aa}'"
    if filter != filters[-1]:
        SQL_QUERY += " and "

explore_results = duckdb.query(SQL_QUERY).pl()

print (f"Result shape: {explore_results.shape}")
print (f"Result columns: {explore_results.columns}")
print (explore_results)


filtered_screen_file_path = f"{local_config['output_root']}/private/yeast_display/{screen_name}__filtered__brotli.parquet"
SQL_QUERY = f"select * from parquet_scan('{filtered_screen_file_path}')"

screen_results = duckdb.query(SQL_QUERY).pl()

screen_results = screen_results.limit(1000)

print (f"Result shape: {screen_results.shape}")
print (f"Result columns: {screen_results.columns}")
print (screen_results)

# join the two results on the peptide column
joined_results = explore_results.join(screen_results, on='peptide', how='inner')

print (f"Joined result shape: {joined_results.shape}")
print (f"Joined result columns: {joined_results.columns}")
print (joined_results)

# filter the joined results to yield only the top 100 peptides

filtered_joined_results = joined_results.limit(100)

print (f"Filtered joined result shape: {filtered_joined_results.shape}")
print (f"Filtered joined result columns: {filtered_joined_results.columns}")
print (filtered_joined_results)


elapsed_time = datetime.now() - start_time
print(f"Elapsed time: {elapsed_time}")

