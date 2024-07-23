import glob
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO

def send_sequence(file):
    seq_ID = file.split("/")[-1][:-6]
    for seq_rec in SeqIO.parse(file, "fasta"):
        sequence = seq_rec.seq
    os.system(f"curl -X POST --data '{sequence}' https://api.esmatlas.com/foldSequence/v1/pdb/ --insecure > {seq_ID}.pdb")
    return seq_ID

def main():
    files = glob.glob("*.fasta")
    max_workers = 8
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_file = {executor.submit(send_sequence, file): file for file in files}
        
        for future in as_completed(future_to_file):
            file = future_to_file[future]
            try:
                seq_ID = future.result()
                print(f"Sequence {seq_ID} processed successfully.")
            except Exception as exc:
                print(f"Sequence {file} generated an exception: {exc}")

if __name__ == "__main__":
    main()
