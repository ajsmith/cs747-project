"""cs747 package configuration"""

from setuptools import setup, find_packages


ENTRY_POINTS = {
    'console_scripts': [
        'cs747-parse-uniprot-fasta=cs747.io:main'
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
