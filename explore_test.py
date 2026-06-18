import duckdb
import polars as pl

from datetime import datetime

from peptide_display.functions import load_config

start_time = datetime.now()

local_config = load_config("local_config.json")

screen_name = 'brenetafusp__hla_a_02_01__yeast_display__9mer__001'

filters = [{"p5": "!H"}, {"p7": "I"}]
filters = [{"p5": "!H"}]


explore_index_file_path = f"{local_config['output_root']}/explore_query_index__brotli.parquet"
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

print (f"SQL_QUERY: {SQL_QUERY} \n")

explore_results = duckdb.query(SQL_QUERY).pl()

explore_peptides = explore_results['peptide'].to_list()

print (f"Result shape: {explore_results.shape}")
print (f"Result columns: {explore_results.columns}")
print (explore_results)


filtered_screen_file_path = f"{local_config['output_root']}//yeast_display_combined__10k__brotli.parquet"
SQL_QUERY = f"select * from parquet_scan('{filtered_screen_file_path}') WHERE screen_slug = '{screen_name}' and peptide in ({', '.join([f'\'{peptide}\'' for peptide in explore_peptides])})"

screen_results = duckdb.query(SQL_QUERY).pl()

screen_results = screen_results.limit(1000)

print (f"Result shape: {screen_results.shape}")
print (f"Result columns: {screen_results.columns}")
print (screen_results)

# filter the joined results to yield only the top 100 peptides

filtered_joined_results = screen_results.limit(100)

print (f"Filtered joined result shape: {filtered_joined_results.shape}")
print (f"Filtered joined result columns: {filtered_joined_results.columns}")
print (filtered_joined_results)


elapsed_time = datetime.now() - start_time
print(f"Elapsed time: {elapsed_time}")

