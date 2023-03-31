"""Test Uniprot data utilities.

"""
from cs747.uniprot import lookup_uniprot_organism, parse_fasta_header


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
