
from Bio.SeqIO.FastaIO import FastaIterator
from pathlib import Path
import pandas as pd

INPUT_DIR = Path(r"data/uniprot_sprot.fasta")

def parsed_fasta_df(input_dir:Path) -> pd.DataFrame:
    """
        Parse the uniprot fasta data and put it into a Dataframe. 
        Paramters: 
            -input_dir(Path): The directory of the fasta data. 
        Return: 
            - A dataframe with the following keys [db, unique_id, entry_name, sequence]
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
    """Test function"""
    df = parsed_fasta_df(INPUT_DIR)
    entry_count = len(df)
    print(f'Parsed {entry_count} entries.')

    df.to_csv('seq.csv')
    print('Saved data to seq.csv')


if __name__ == "__main__":
    main()
