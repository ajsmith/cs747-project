# cs747-project

## Team Roster

- Alexander Smith
- Kelvin Lu
- Archange Giscard Destiné


# Requirements

- Python 3


# Installation

To install, run:

```shell
$ install.sh
```

This creates a Python virtual env in the directory `venv`.

Activate the environment with the command:

```shell
$ source venv/bin/activate
```

# Usage

## Downloading sequence data

Run the `download_data.sh` script.

This downloads the Uniprot Swiss-Prot FASTA data and unpacks it to
`data/uniprot_sprot.fasta`.


## Generating the sequence data CSV

Run `cs747-parse-uniprot-fasta`

This parses Uniprot FASTA data into a Pandas dataframe and saves it as
a CSV. This creates the file, `data/seq.csv`


## Building the taxonomy database

Run `cs747-build-taxonomy-db`

This populates the taxonomy database from sequence data contained in
the sequence CSV by looking up organism entries from the Uniprot
Taxonomy REST API. It then saves it as a Python Pickle file, named
`data/taxonomy_db.pickle`.


## Creating the labeled dataset

Run `cs747-label-data`

This labels the sequence data and balances the data for our classes.
It saves the labeled data as a new CSV, named
`data/labeled_sequences.csv`.
