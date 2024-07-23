#Simple script to predict protein structures from sequence with ESMFold API
#Adjust the number of workers as needed!

import glob
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO
import time

def send_sequence(file):
    seq_ID = file.split("/")[-1][:-6]
    for seq_rec in SeqIO.parse(file, "fasta"):
        sequence = seq_rec.seq
    output_file = f"./{seq_ID}.pdb"
    
    # Send the sequence to the server and save the output
    os.system(f"curl -X POST --data '{sequence}' https://api.esmatlas.com/foldSequence/v1/pdb/ --insecure > {output_file}")
    
    return seq_ID, output_file

def check_file_size(file):
    """ Check if the file size is 23 bytes and return True if it is. """
    #This is needed, the 23b files indicate the server gave you temporary ban
    return os.path.exists(file) and os.path.getsize(file) == 23

def main(done):
    files = glob.glob("./*.fasta")
    
    # Filter out files that have already been processed
    files = [file for file in files if file.split("/")[-1][:-6] not in done]
    
    max_workers = 3 # adjust as needed
    delay_between_requests = 2  # seconds
    max_retries = 5
    retry_delay = 180  # seconds to wait before retrying
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_file = {executor.submit(send_sequence, file): file for file in files}
        
        for future in as_completed(future_to_file):
            file = future_to_file[future]
            try:
                seq_ID, output_file = future.result()
                print(f"Sequence {seq_ID} processed successfully.")
                
                # Check the size of the output file
                if check_file_size(output_file):
                    print(f"Output file {output_file} is 23 bytes. Pausing for {retry_delay} seconds...")
                    time.sleep(retry_delay)  # Pause if file is 23 bytes
                else:
                    time.sleep(delay_between_requests)  # Regular delay
                
            except Exception as exc:
                print(f"Sequence {file} generated an exception: {exc}")

if __name__ == "__main__":

    # List of finished predictions
    done = []

    # Read the output folder to populate the `done` list
    for file in glob.glob("./*.pdb"):
        if os.path.getsize(file) > 50:
            with open(file, "r") as f:
                header_found = any(line.startswith("HEADER") for line in f)
                if header_found:
                    seq_ID = file.split("/")[-1][:-4]
                    done.append(seq_ID)

    main(done=done)
