#!/usr/bin/env bash

# download
curl -LO https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz
mkdir -pv data
gzip -dv uniprot_sprot.fasta.gz
mv -v uniprot_sprot.fasta data/
