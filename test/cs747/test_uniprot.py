"""Test Uniprot data utilities.

"""
from tempfile import NamedTemporaryFile

from cs747.uniprot import (
    Labeler,
    TaxonomyDatabaseBuilder,
    load_taxonomy_db,
    lookup_uniprot_organism,
    parse_fasta_header,
)


def test_parse_fasta_header():
    """Test the FASTA header parser."""

    s1 = 'sp|Q6GZX4|001R_FRG3G Putative transcription factor 001R OS=Frog virus 3 (isolate Goorha) OX=654924 GN=FV3-001R PE=4 SV=1'
    fields = parse_fasta_header(s1)
    assert fields.db == 'sp'
    assert fields.unique_id == 'Q6GZX4'
    assert fields.entry_name == '001R_FRG3G'
    assert fields.protein_name == 'Putative transcription factor 001R'
    assert fields.organism_name == 'Frog virus 3 (isolate Goorha)'
    assert fields.organism_id == '654924'


def test_lookup_uniprot_organism():
    """Test lookup of an organism entry on Uniprot."""
    organism_id = '402880'
    entry = lookup_uniprot_organism(organism_id)
    assert str(entry['taxonId']) == organism_id
    assert entry['scientificName'] == "Methanococcus maripaludis (strain C5 / ATCC BAA-1333)"
    assert entry['lineage'][-2]['scientificName'] == "Archaea"


def test_create_taxonomy_db(sequence_csv_path, taxonomy_db):
    """Test taxonomy DB creation."""
    tmp_file = NamedTemporaryFile(prefix='cs747.test-')
    db_file_path = tmp_file.name
    tmp_file.close()

    builder = TaxonomyDatabaseBuilder(db_file_path, sequence_csv_path)
    builder.populate(save_interval=3)
    new_tax_db = load_taxonomy_db(db_file_path)
    assert new_tax_db == taxonomy_db


def test_init_labeler(taxonomy_db_path, sequence_csv_path):
    """Test creating a labeler."""
    labeler = Labeler(taxonomy_db_path, sequence_csv_path)
    assert hasattr(labeler, 'seq_df')
    assert hasattr(labeler, 'tax_db')


def test_label_viruses(labeler):
    """Test labeling viruses."""
    seq_df = labeler.seq_df
    virus_idx = seq_df['organism_id'].map(labeler.is_virus)
    viruses = seq_df[virus_idx]

    assert len(viruses) == 1
    assert (viruses['organism_name'] == 'Frog virus 3 (isolate Goorha)').all()

    labels = seq_df['organism_id'].map(labeler.label_organism)
    assert labels[virus_idx].to_list() == ['Viruses']


def test_label_sequences(labeler):
    """Test labeling of sequence data."""
    labeler.label_sequences()
    seq_df = labeler.seq_df
    label = seq_df['label']
    organism_name = seq_df['organism_name']

    classes = (
        "Viruses",
        "Bacteria",
        "Archaea",
        "Viridiplantae",
        "Fungi",
        "Chordata",
        "Metazoa",
        "Eukaryota",
    )

    results = {c:(organism_name[label == c].to_list()) for c in classes}
    assert list(sorted(results.keys())) == list(sorted(classes)), (
        "All classes must be represented."
    )

    assert results["Viruses"] == ['Frog virus 3 (isolate Goorha)']
    assert results["Bacteria"] == ['Bacillus subtilis (strain 168)']
    assert results["Archaea"] == ['Metallosphaera sedula (strain ATCC 51363 / DSM 5348 / JCM 9185 / NBRC 15509 / TH2)']
    assert results["Viridiplantae"] == ['Oryza sativa subsp. japonica']
    assert results["Fungi"] == ['Saccharomyces cerevisiae (strain ATCC 204508 / S288c)']
    assert results["Chordata"] == ['Homo sapiens']
    assert results["Metazoa"] == ['Drosophila melanogaster']
    assert results["Eukaryota"] == ['Dictyostelium discoideum']
