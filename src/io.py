
from Bio.SeqIO.FastaIO import FastaIterator
from pathlib import Path
import pandas 

INPUT_DIR = Path(r"data/uniprot_sprot.fasta")

def parse_fasta(input_dir:Path):

    with open(input_dir) as handle:
        for index,value in enumerate(FastaIterator(handle)):
            print(value)
            if index == 5: 
                break
            



def main(): 
    """ Test function """

    print("Foo bar")
    parse_fasta(INPUT_DIR)

if __name__ == "__main__":
    main()
