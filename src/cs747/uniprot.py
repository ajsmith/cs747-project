"""Utilities for working with Uniprot data and APIs.

"""
from collections import namedtuple
from pathlib import Path
from typing import Union
from urllib.request import urlopen
import json
import pickle
import random

from Bio.SeqIO.FastaIO import FastaIterator
import pandas as pd

FilePath = Union[str, Path]

FASTA_FILE_PATH = Path(r"data/uniprot_sprot.fasta")
TAXONOMY_DB_PATH = Path(r'data/taxonomy_db.pickle')
SEQUENCE_CSV_PATH = Path(r'data/seq.csv')

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


def parse_fasta_df(input_dir: Path) -> pd.DataFrame:
    """
    Parse the uniprot fasta data and put it into a Dataframe.
    Paramters:
        - input_dir(Path): The directory of the fasta data.
    Return:
        - A dataframe with the following columns:
            * db
            * unique_id
            * entry_name
            * protein_name
            * organism_name
            * organism_id
            * sequence
    """

    data: dict[str, list] = {
        'db': [],
        'unique_id': [],
        'entry_name': [],
        'protein_name': [],
        'organism_name': [],
        'organism_id': [],
        'sequence': [],
    }

    with open(input_dir) as handle:
        for index, value in enumerate(FastaIterator(handle)):
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


def load_taxonomy_db(
        db_file_path: FilePath = TAXONOMY_DB_PATH
) -> dict:
    """Load the taxonomy DB from a backing file and return it."""
    with open(db_file_path, 'rb') as db_file:
        result = pickle.load(db_file)
        return result


class TaxonomyDatabaseBuilder:
    """Utility to create a Taxonomy DB from Uniprot data."""

    def __init__(
            self,
            db_file_path=TAXONOMY_DB_PATH,
            seq_csv_path=SEQUENCE_CSV_PATH,
            recreate=False,
    ):
        """Create a TaxonomyDatabaseBuilder."""

        self.db_file_path = db_file_path
        self.seq_csv_path = seq_csv_path
        self.init_db(recreate=recreate)

    def populate(self, save_interval=100):
        """Populate DB from organisms in the sequence dataframe."""

        seq_df = pd.read_csv(self.seq_csv_path)
        n_seq = len(seq_df)
        count = 0
        prev_count = 0

        print(f'Populating taxonomy DB from {n_seq} sequence entries.')

        for organism_id in seq_df['organism_id']:

            if organism_id not in self.db:
                entry = lookup_uniprot_organism(organism_id)
                self.db[organism_id] = entry
                count += 1

            if count > prev_count and count % save_interval == 0:
                self.save()
                entry_count = count - prev_count
                prev_count = count
                print(f'{entry_count} organisms added to taxonomy DB.')

        self.save()
        db_size = len(self.db)
        print(f'{count} organisms added to taxonomy DB.')
        print(f'Taxonomy DB contains {db_size} total organisms.')

    def load(self):
        """Load DB from backing file and return it."""
        result = load_taxonomy_db(self.db_file_path)
        return result

    def init_db(self, recreate=False):
        """Instantiate the DB."""
        db = None

        if recreate is False:
            try:
                db = self.load()
            except FileNotFoundError:
                print('No database found.')
            else:
                print(f'Loaded taxonomy DB from {self.db_file_path}')

        if db is None:
            print('Initializing new taxonomy DB')
            db = {}

        self.db = db

    def save(self):
        """Save DB to backing file."""
        with open(self.db_file_path, 'wb') as db_file:
            pickle.dump(self.db, db_file)


def build_taxonomy_db(
        db_file_path=TAXONOMY_DB_PATH,
        seq_csv_path=SEQUENCE_CSV_PATH
):
    """Create the taxonomy database."""
    db_builder = TaxonomyDatabaseBuilder(db_file_path, seq_csv_path)
    db_builder.populate(save_interval=100)
    print(f'Saved taxonomy DB to {db_file_path}')


def import_fasta_to_csv():
    """Test function"""
    df = parse_fasta_df(FASTA_FILE_PATH)
    entry_count = len(df)
    print(f'Parsed {entry_count} entries.')

    df.to_csv('data/seq.csv', index=False)
    print('Saved data to data/seq.csv')


def create_test_data(data_dir='test/cs747/data') -> None:
    """Recreate test data from a sample FASTA file."""
    data_dir = Path(data_dir)
    fasta_file_path = data_dir / 'uniprot_sprot.fasta'
    seq_csv_path = data_dir / 'seq.csv'
    tax_db_file_path = data_dir / 'taxonomy_db.pickle'

    seq_df = parse_fasta_df(fasta_file_path)
    seq_df.to_csv(seq_csv_path, index=False)

    builder = TaxonomyDatabaseBuilder(
        tax_db_file_path, seq_csv_path, recreate=True)
    builder.populate()
    print(f"Test data created in {data_dir}")


def build_percentage_label_stats(data_dir:Path = FASTA_FILE_PATH, label:str = "label") -> dict:

    df = parse_fasta_df(data_dir)

    # TODO - Fake labels (delete following block)
    fake_labels = [f"label{x}" for x in range(8)]
    df[label] = df["db"].apply(lambda x: random.choice(fake_labels))
    
    # Count all of the label and divide it by population. 
    percentage_df = df[label].value_counts().apply(lambda x: x/len(df))

    # Convert df into dict. 
    percentage_dict = percentage_df.to_dict()

    # print(percentage_dict)

    return percentage_dict


