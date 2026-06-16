import duckdb
import polars as pl
import os

from peptide_display.functions import load_config, kmerize_sequence


def fast_fuzzy_match(target, peptide, file_path):
    search_kmers = kmerize_sequence(peptide, 3)
    
    print(f"Search kmers for peptide {peptide} / {target}: {search_kmers}")

    kmer_sql = [f"kmer = '{kmer}' AND position = {i}" for i, kmer in enumerate(search_kmers)]

    sql_query = f"SELECT count(kmer) as count, screen_slug FROM '{file_path}' WHERE {' OR '.join(kmer_sql)} GROUP BY screen_slug ORDER BY count DESC"

    con = duckdb.connect()
    results = con.execute(sql_query).pl()

    results = results.filter(pl.col("count") >= 2)
    print(results)
    print ("\n\n")


peptides = [
    {'target': 'MAGE-A3', 'peptide': 'ESDPIGHLY'}, 
    {'target': 'MAGE-A6', 'peptide': 'EMDPFHYDY'},
    {'target': 'PLD5', 'peptide': 'ETDPLTFNF'},
    {'target': 'hFat2', 'peptide': 'ETDPVNHMV'},
    {'target': 'MAGE-B18', 'peptide': 'EVDPIRHYY'},
    {'target': 'CD166', 'peptide': 'EMDPVTQLY'},
    {'target': 'STK/MRCKA', 'peptide': 'ETDPVENTY'},
    {'target': 'Titin', 'peptide': 'ESDPIVAQY'}, 
    {'target': 'ANKRD16', 'peptide': 'EGDPLILQY'}, 
    {'target': 'COG4', 'peptide': 'ELDPILTEV'},
    {'target': 'FXYD6', 'peptide': 'EMDPFHYDY'},
    {'target': 'DCTN4', 'peptide': 'EVEPLPEDY'},
    {'target': 'PRAME', 'peptide': 'SLLQHLIGL'}, 
    {'target': 'GP100', 'peptide': 'YLEPGPVTA'},
    {'target': 'CMV', 'peptide': 'VTEHDTLLY'},
    {'target': 'EBV', 'peptide': 'GLCTLVAML'},
    {'target': 'NY-ESO-1', 'peptide': 'SLLMWITQC'},
    {'target': 'PSG5', 'peptide': 'GRLPLLNPI'}, 
    {'target': 'YEIH', 'peptide': 'LRVMMLAPF'}
]


local_config = load_config("local_config.json")

file_path = f"{local_config['output_root']}/search_index__brotli.parquet"

file_size = os.path.getsize(file_path) / (1024 * 1024)  # Size in MB


sql_query = f"SELECT position, kmer, screen_slug FROM '{file_path}'"
con = duckdb.connect()
search_df = con.execute(sql_query).pl()


print (f"\n\n\nTotal unique kmers in mini index: {search_df.select(pl.col('kmer').n_unique()).item()}" )
print (f"Total unique screen_slugs in mini index: {search_df.select(pl.col('screen_slug').n_unique()).item()} \n" )
print(f"Size of the mini index parquet file: {file_size:.2f} MB \n")


print ("Starting fuzzy matching of peptides against a combined mini index...\n")


for peptide in peptides:

    fast_fuzzy_match(peptide['target'], peptide['peptide'], file_path)

