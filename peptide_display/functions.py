import json
import csv
import re


## CONFIG LOADING ##
def load_config(config_file: str ="config.json") -> dict:
    with open(config_file) as f:
        config = json.load(f)
    return config   


## CONSTANTS ## 
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'


## UTILITY FUNCTIONS ##
def is_valid_peptide(peptide):
    """
    Check if the given peptide sequence is valid by ensuring it only contains standard amino acid characters.
    Args:
        peptide (str): The peptide sequence to validate.
    Returns:
        bool: True if the peptide is valid, False otherwise.
    """
    if len(peptide) == 0:
        return False
    return not re.search(rf'[^{AMINO_ACIDS}]', peptide)


def map_column_name(config: dict, screen_type: str, column: str) -> str:
    """
    Map the given column name to its corresponding name in the configuration for the specified screen type.
    Args:
        config (dict): The configuration dictionary.
        screen_type (str): The type of screen to use for mapping.
        column (str): The column name to map.
    Returns:
        str: The mapped column name if found, otherwise None.
    """
    if column in config[screen_type]['columns']['standardised_columns']:
        return column
    elif column in config[screen_type]['columns']['mappings']:
        return config[screen_type]['columns']['mappings'][column]
    else:
        return None


def slugify(text):
    """
    Convert text to a slug suitable for use as a column name.
    
    Args:
        text (str): The input text to slugify.
        
    Returns:
        str: The slugified text.
    """
    return text.lower().replace(' ', '_')


def load_screen_slugs(config: dict, screen_type: str, availability: str = 'private') -> list:
    """
    Load the screen slugs from the screens.csv file based on the specified availability.
    
    Args:
        config (dict): The configuration dictionary containing the output root.
        screen_type (str): The type of screen to load slugs for.
        availability (str): The availability type ('private' or 'public').
        
    Returns:
        list: A list of screen slugs.
    """
    screens_csv = f"{config['output_root']}/screens.csv"

    slug_lookup_dict = {}
    with open(screens_csv, 'r') as f:
        reader = csv.DictReader(f)
        slug_lookup_dict = {row['original_rawdata_file']: row['screen_slug'] for row in reader if len(row['screen_slug']) > 0 and row['screen_type'] == screen_type}
    return slug_lookup_dict


## DATA LOOKUP FUNCTIONS ##
def lookup_screen_slug(config: dict, raw_file_name: str, screen_type: str) -> str:
    """
    Lookup the slug for a given raw file name using the provided slug lookup dictionary.
    
    Args:
        config (dict): The configuration dictionary containing the slug lookup.
        raw_file_name (str): The name of the raw file to look up.
        screen_type (str): The type of screen to use for lookup.
        
    Returns:
        str: The corresponding slug for the given raw file name, or None if not found.
    """
    slug_lookup_dict = load_screen_slugs(config, screen_type)

    return slug_lookup_dict.get(raw_file_name, None)


## Kmerization ##

def kmerize_sequence(sequence:str, k:int) -> list:
    """
    Generate a list of k-mers for a given sequence.
    
    Args:
        sequence (str): A single sequence (string).
        k (int): The length of the k-mers to generate.
        
    Returns:
        list: A list of k-mers.
    """
    kmers = [sequence[i:i + k] for i in range(len(sequence) - k + 1)]
    return kmers    

