"""cs747 package configuration"""

from setuptools import setup, find_packages


ENTRY_POINTS = {
    'console_scripts': [
        'cs747-parse-uniprot-fasta=cs747.uniprot:import_fasta_to_csv',
        'cs747-build-taxonomy-db=cs747.uniprot:build_taxonomy_db',
        'cs747-create-test-data=cs747.uniprot:create_test_data',
        'cs747-label-data=cs747.uniprot:label_sequences',
    ]
}

AUTHORS = [
    'Alexander Smith',
    'Kelvin Lu',
    'Archange Destine',
]

AUTHOR_EMAILS = [
    'asmitl@gmu.edu',
    'klu21@gmu.edu',
    'adestine@gmu.edu',
]

setup(
    name='cs747',
    version='0',
    author=', '.join(AUTHORS),
    author_email=', '.join(AUTHOR_EMAILS),
    url='https://github.com/ajsmith/cs747-project',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    entry_points=ENTRY_POINTS,
    install_requires=[
        'PyYAML',
        'biopython',
        'pandas',
    ],
)
