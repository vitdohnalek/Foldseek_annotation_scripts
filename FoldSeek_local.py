import glob
import requests
import time
import os


# Get all CIF files in the current directory
cif_files = glob.glob("*.cif") # pdb fromat works too

# Process files in batches of 5
batch_size = 5
for i in range(0, len(cif_files), batch_size):
    batch = cif_files[i:i + batch_size]
    print(f"\nUploading batch {i // batch_size + 1}:")

    for file in batch:
        seq_ID = file.split("/")[-1][:-4]
        print(f"  Uploading: {file}")
        with open(file, 'rb') as f:
            response = requests.post(
                "https://search.foldseek.com/api/ticket",
                files={"q": f},
                data={
                    "mode": "3diaa",
                    "database[]": ["BFVD", "afdb50", "afdb-swissprot", "afdb-proteome"]
                }
            )

        if response.status_code == 200:
            result = response.json()
            ticket_id = result["id"]
            print(f"    Ticket ID: {ticket_id}")

            # Polling for status
            status_url = f"https://search.foldseek.com/api/ticket/{ticket_id}"
            while True:
                status_response = requests.get(status_url)
                status_data = status_response.json()
                status = status_data["status"]
                print(f"    Status: {status}")
                if status == "COMPLETE":
                    print("Job complete.")
                    os.system(f"curl -L https://search.foldseek.com/api/result/download/{ticket_id} -o {seq_ID}.gz")
                    break
                elif status == "ERROR":
                    print("Error occurred.")
                    break
                time.sleep(10)  # wait X seconds before checking again
        elif response.status_code == 429:
            time.sleep(120)
        else:
            print(f"Failed to upload {file}, status code: {response.status_code}")

    time.sleep(20)
