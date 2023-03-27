#/usr/bin/env bash ./download_data.sh

# download 
wget https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz
mkdir data
gzip -d uniprot_sprot.fasta.gz  
mv uniprot_sprot.fasta data
rm uniprotq_sprot.fasta