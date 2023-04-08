from pathlib import Path

import pandas as pd
import pytest

from cs747.uniprot import load_taxonomy_db


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
    result = pd.read_csv(sequence_csv_path)
    return result


@pytest.fixture
def taxonomy_db(taxonomy_db_path: Path) -> dict[str, dict]:
    """Return a sample taxonomy database."""
    result = load_taxonomy_db(taxonomy_db_path)
    return result
