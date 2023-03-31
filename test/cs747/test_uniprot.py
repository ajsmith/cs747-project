from cs747.uniprot import parse_fasta_header

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
