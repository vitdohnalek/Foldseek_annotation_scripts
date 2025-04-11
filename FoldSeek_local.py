import glob
import requests
import time
import os


# Extracts results and deleted unnecessary files
def get_results(compressed_file=""):
    gz_file = compressed_file
    seq_ID = gz_file[:-3]

    os.system(f"tar -xvf {gz_file}")

    # Remove unnecesarry files
    os.system(f"rm {gz_file}")
    os.system("rm alis_afdb-swissprot_report.m8")
    os.system("rm alis_afdb50_report.m8")
    os.system("rm alis_afdb-proteome_report.m8")
    os.system("rm alis_BFVD_report.m8")
    os.system("rm alis_BFVD.m8")

    def check_db_results(file=""):

        best_hit = "None"
        all_hits = []

        with open(file, "r") as f:  
            for l in f:
                if not "\t" in l:

                    return best_hit, all_hits

                else:
                    line = l.split("\t")
                    annotation = " ".join(line[1].split()[1:])
                    probability = float(line[10])

                    if not "uncharacterized" in annotation.lower() and not "hypothetical" in annotation.lower() and not "putative" in annotation.lower() and not "predicted" in annotation.lower():
                        if probability >= 0.5 and best_hit == "None":
                            best_hit = annotation
                            all_hits.append(annotation)
                        elif probability >= 0.5:
                            all_hits.append(annotation)

        return best_hit, all_hits

    best_swissprot_hit, swiss_hits = check_db_results(file="alis_afdb-swissprot.m8")
    best_afdb50_hit, afdb50_hits = check_db_results(file="alis_afdb50.m8")
    best_proteomes_hit, proteome_hits = check_db_results(file="alis_afdb-proteome.m8")

    all_hits = swiss_hits + afdb50_hits + proteome_hits

    if len(all_hits) > 0:
        # Most common element finder
        most_common_hit = max(set(all_hits), key=all_hits.count)
        most_common_hit_n = str(all_hits.count(most_common_hit))
    else:
        most_common_hit = "None"
        most_common_hit_n = "0"

    os.system("rm alis_afdb50.m8")
    os.system("rm alis_afdb-swissprot.m8")
    os.system(f"rm alis_afdb-proteome.m8")

    tsv = "Protein ID\tSwiss-Prot\tUniProt\tAlphaFold-Proteomes\tMost frequent hit\tMost frequent hit n\n"
    tsv += seq_ID + "\t" + best_swissprot_hit + "\t"
    tsv += best_afdb50_hit + "\t" + best_proteomes_hit + "\t"
    tsv += most_common_hit + "\t" + most_common_hit_n + "\n"

    with open(f"{seq_ID}_foldseek.tsv", "w") as f:
        f.write(tsv)

# Get all CIF files in the current directory
cif_files = glob.glob("./done/*.cif") # pdb fromat works too

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
                    get_results(f"{seq_ID}.gz")
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
