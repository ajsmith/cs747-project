"""Utilities for working with Uniprot data and APIs.

"""
from collections import namedtuple
from functools import cache
from pathlib import Path
from pprint import pprint
from urllib.request import urlopen
import json
import pickle
import random

from Bio.SeqIO.FastaIO import FastaIterator
import pandas as pd


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


def load_taxonomy_db(db_file_path=TAXONOMY_DB_PATH) -> dict[str, dict]:
    """Load the taxonomy DB from a backing file and return it."""
    with open(db_file_path, 'rb') as db_file:
        result = pickle.load(db_file)
        return result


def load_sequence_df(csv_path=SEQUENCE_CSV_PATH) -> pd.DataFrame:
    """Load a sequence dataframe from a CSV file and return it."""
    result = pd.read_csv(csv_path)
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


class Labeler:
    """Create labels for sequence data."""

    def __init__(
            self,
            db_file_path=TAXONOMY_DB_PATH,
            seq_csv_path=SEQUENCE_CSV_PATH,
    ):
        """Create a Labeler."""
        self.tax_db = load_taxonomy_db(db_file_path)
        self.seq_df = load_sequence_df(seq_csv_path)

    def lineage_by_name(self, organism_id):
        """Return the lineage as a list of scientific names.

        The returned list is ordered from most specific to least
        specific.
        """
        entry = self.tax_db[organism_id]
        result = [t['scientificName'] for t in entry['lineage']]
        return result

    def has_lineage(self, organism_id, tax_name, tax_rank=None):
        """Return True if the organism has the given taxonomy.

        If `tax_rank` is not None, returns True only if the taxonomic
        name is found at the given rank index.
        """
        lineage = self.lineage_by_name(organism_id)
        if tax_rank is None:
            result = tax_name in lineage
        else:
            result = lineage[tax_rank] == tax_name
        return result

    def is_virus(self, organism_id):
        """Return True if the organism is a virus."""
        result = self.has_lineage(organism_id, "Viruses", -1)
        return result

    def is_cellular_organism(self, organism_id):
        """Return True if the organism is a cellular organism."""
        result = self.has_lineage(organism_id, "cellular organisms", -1)
        return result

    def is_eukaryote(self, organism_id):
        """Return True if the organism is a eukaryote."""
        result = False
        if self.is_cellular_organism(organism_id):
            result = self.has_lineage(organism_id, "Eukaryota", -2)
        return result

    def is_bacteria(self, organism_id):
        """Return True if the organism is a bacteria."""
        result = False
        if self.is_cellular_organism(organism_id):
            result = self.has_lineage(organism_id, "Bacteria", -2)
        return result

    def is_archaea(self, organism_id):
        """Return True if the organism is an archaea."""
        result = False
        if self.is_cellular_organism(organism_id):
            result = self.has_lineage(organism_id, "Archaea", -2)
        return result

    def is_viridiplantae(self, organism_id):
        """Return True if the organism is in Viridiplantae."""
        result = False
        if self.is_eukaryote(organism_id):
            result = self.has_lineage(organism_id, "Viridiplantae")
        return result

    def is_fungi(self, organism_id):
        """Return True if the organism is in Fungi."""
        result = False
        if self.is_eukaryote(organism_id):
            result = self.has_lineage(organism_id, "Fungi")
        return result

    def is_chordata(self, organism_id):
        """Return True if the organism is in Chordata."""
        result = False
        if self.is_eukaryote(organism_id):
            result = self.has_lineage(organism_id, "Chordata")
        return result

    def is_metazoa(self, organism_id):
        """Return True if the organism is in Metazoa."""
        result = False
        if self.is_eukaryote(organism_id) and not self.is_chordata(organism_id):
            result = self.has_lineage(organism_id, "Metazoa")
        return result

    @cache
    def label_organism(self, organism_id):
        """Return the label for an organism."""
        result = None
        if self.is_virus(organism_id):
            result = "Viruses"
        elif self.is_bacteria(organism_id):
            result = "Bacteria"
        elif self.is_archaea(organism_id):
            result = "Archaea"
        elif self.is_viridiplantae(organism_id):
            result = "Viridiplantae"
        elif self.is_fungi(organism_id):
            result = "Fungi"
        elif self.is_chordata(organism_id):
            result = "Chordata"
        elif self.is_metazoa(organism_id):
            result = "Metazoa"
        else:
            result = "Eukaryota"
        return result

    def label_sequences(self):
        """Label sequence data."""
        labels = self.seq_df['organism_id'].map(self.label_organism)
        self.seq_df['label'] = labels


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


def label_sequences():
    """Label and balance the sequence data."""
    labeler = Labeler(TAXONOMY_DB_PATH, SEQUENCE_CSV_PATH)
    print(f"Labeling sequence data from {SEQUENCE_CSV_PATH}")
    labeler.label_sequences()

    print("Label statistics:")
    label_stats = build_percentage_label_stats(labeler.seq_df)
    pprint(label_stats)

    print("Balancing data")
    balanced = generate_balanced_data(labeler.seq_df)
    balanced.to_csv('data/labeled_sequences.csv', index=False)
    print("Wrote labeled data to data/labeled_sequences.csv")

    print("Balanced data statistics:")
    balanced_stats = build_percentage_label_stats(balanced)
    pprint(balanced_stats)


def generate_fake_lbl(df, header):
    """Delete this function"""
    # TODO - Delete function
    fake_labels = [f"label{x}" for x in range(8)]
    df[header] = df["db"].apply(lambda x: random.choice(fake_labels))

    return df


def build_percentage_label_stats(
        data: pd.DataFrame, header: str = "label"
) -> dict:

    # Count all of the label and divide it by population.
    n = len(data)
    percentage_df = data[header].value_counts().apply(
        lambda x: (x / n, n * x / n)
    )

    # Convert df into dict.
    percentage_dict = percentage_df.to_dict()

    return percentage_dict


def generate_balanced_data(
        data: pd.DataFrame,
        header: str = "label",
        frac_population: float = 0.01639,
) -> pd.DataFrame:
    """Genereate a balanced dataset based on the fraction of the population."""
    sample_size = round(frac_population * len(data))
    output_df = data.groupby(header).sample(sample_size)
    return output_df
