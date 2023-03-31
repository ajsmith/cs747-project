"""Utilities for working with Uniprot data and APIs.

"""
from collections import namedtuple
from pathlib import Path
from urllib.request import urlopen
import json
import pickle

from Bio.SeqIO.FastaIO import FastaIterator
import pandas as pd


INPUT_DIR = Path(r"data/uniprot_sprot.fasta")
TAXONOMY_DB_PATH = Path(r'data/taxonomy_db.pickle')


HeaderFields = namedtuple(
    'HeaderFields',
    [
        'db',
        'unique_id',
        'entry_name',
        'protein_name',
        'organism_name',
        'organism_id',
    ],
)


def parse_fasta_header(raw_header: str) -> HeaderFields:
    """Parse a Uniprot FASTA entry header.

    """
    id_block, remaining = raw_header.split(None, 1)
    db, uid, entry = id_block.split('|')

    os_idx = remaining.index('OS=')
    ox_idx = remaining.index('OX=')
    protein = remaining[:os_idx].rstrip()
    organism = remaining[os_idx:ox_idx].split('=')[1].rstrip()
    organism_id = remaining[ox_idx:].split(None, 1)[0].split('=')[1]

    result = HeaderFields(db, uid, entry, protein, organism, organism_id)
    return result


def parsed_fasta_df(input_dir:Path) -> pd.DataFrame:
    """
    Parse the uniprot fasta data and put it into a Dataframe.
    Paramters:
        -input_dir(Path): The directory of the fasta data.
    Return:
        - A dataframe with the following keys [db, unique_id, entry_name, sequence]
    """

    data = {
        'db': [],
        'unique_id': [],
        'entry_name': [],
        'protein_name': [],
        'organism_name': [],
        'organism_id': [],
        'sequence': [],
    }

    with open(input_dir) as handle:
        for index,value in enumerate(FastaIterator(handle)):
            fields = parse_fasta_header(value.description)
            data['db'].append(fields.db)
            data['unique_id'].append(fields.unique_id)
            data['entry_name'].append(fields.entry_name)
            data['protein_name'].append(fields.protein_name)
            data['organism_name'].append(fields.organism_name)
            data['organism_id'].append(fields.organism_id)
            data['sequence'].append(value.seq)

    output = pd.DataFrame(data)
    return output


def lookup_uniprot_organism(organism_id: str) -> dict:
    """Query Uniprot for an organism by ID.

    """
    result = None
    base_url = 'https://rest.uniprot.org/taxonomy/'
    query_url = f'{base_url}{organism_id}.json'
    response = urlopen(query_url)
    result = json.loads(response.read())
    return result


def populate_taxonomy_db(seq_df, tax_db):
    """Populate the taxonomy db from organisms in the sequence dataframe."""
    report_size = 10
    processed = []

    for organism_id in seq_df['organism_id']:

        if organism_id not in tax_db:
            entry = lookup_uniprot_organism(organism_id)
            tax_db[organism_id] = entry
            processed.append(organism_id)

        if len(processed) >= report_size:
            print(f'Added {processed}')
            processed = []


def load_taxonomy_db(db_file_path=TAXONOMY_DB_PATH):
    """Load a taxonomy DB from a backing file."""
    with open(db_file_path, 'rb') as db_file:
        result = pickle.load(db_file)
        return result


def build_taxonomy_db(db_file_path=TAXONOMY_DB_PATH):
    """Create the taxonomy database."""
    seq_df = pd.read_csv('data/seq.csv')

    try:
        tax_db = load_taxonomy_db(db_file_path)
    except FileNotFoundError:
        print('Initializing new taxonomy DB')
        tax_db = {}
    else:
        print(f'Loaded taxonomy DB from {db_file_path}')

    populate_taxonomy_db(seq_df, tax_db)

    with open(db_file_path, 'wb') as db_file:
        pickle.dump(tax_db, db_file)
        print(f'Saved taxonomy DB to {db_file_path}')


def import_fasta_to_csv():
    """Test function"""
    df = parsed_fasta_df(INPUT_DIR)
    entry_count = len(df)
    print(f'Parsed {entry_count} entries.')

    df.to_csv('data/seq.csv', index=False)
    print('Saved data to data/seq.csv')
