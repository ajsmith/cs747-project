from pathlib import Path

import pandas as pd
import pytest

from cs747.uniprot import Labeler, load_sequence_df, load_taxonomy_db


@pytest.fixture
def sequence_csv_path() -> Path:
    """Return a path to a sample sequence CSV file."""
    return Path('test/cs747/data/seq.csv')


@pytest.fixture
def taxonomy_db_path() -> Path:
    """Return a path to a sample taxonomy DB backing file."""
    return Path('test/cs747/data/taxonomy_db.pickle')


@pytest.fixture
def sequence_df(sequence_csv_path: Path) -> pd.DataFrame:
    """Return a sample sequence dataframe."""
    result = load_sequence_df(sequence_csv_path)
    return result


@pytest.fixture
def taxonomy_db(taxonomy_db_path) -> dict[str, dict]:
    """Return a sample taxonomy database."""
    result = load_taxonomy_db(taxonomy_db_path)
    return result

@pytest.fixture
def labeler(taxonomy_db_path, sequence_csv_path) -> Labeler:
    """Return a labeler instance."""
    result = Labeler(taxonomy_db_path, sequence_csv_path)
    return result
