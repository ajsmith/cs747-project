
from Bio.SeqIO.FastaIO import FastaIterator
from pathlib import Path
import pandas as pd

INPUT_DIR = Path(r"data/uniprot_sprot.fasta")

def uniprot_fasta_to_df(input_dir:Path) -> pd.DataFrame:
    """
        Parse the uniprot fasta data and put it into a Dataframe. 
    """
    # output = pd.DataFrame()
    db = []
    unique_id = []
    entry_name = []
    sequence = []

    with open(input_dir) as handle:
        for index,value in enumerate(FastaIterator(handle)):
            ids = value.name.split("|")

            db.append(ids[0])
            unique_id.append(ids[1])
            entry_name.append(ids[2])
            sequence.append(value.seq)
            
            # if index == 10:
                # break

    output = pd.DataFrame({
        "db":db, 
        "unique_id": unique_id, 
        "entry_name": entry_name,
        "sequence": sequence
    })

    # print(output)
    return output



def main(): 
    """ Test function """

    uniprot_fasta_to_df(INPUT_DIR)

if __name__ == "__main__":
    main()
