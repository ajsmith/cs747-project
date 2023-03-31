"""Utilities for working with Uniprot data and APIs.

"""
from collections import namedtuple
from pathlib import Path

from Bio.SeqIO.FastaIO import FastaIterator
import pandas as pd


INPUT_DIR = Path(r"data/uniprot_sprot.fasta")


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


def import_fasta_to_csv():
    """Test function"""
    df = parsed_fasta_df(INPUT_DIR)
    entry_count = len(df)
    print(f'Parsed {entry_count} entries.')

    df.to_csv('data/seq.csv', index=False)
    print('Saved data to data/seq.csv')
