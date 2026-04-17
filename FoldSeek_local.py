import glob
import requests
import time
import os


# Extracts results and returns best hits + all hits per database
def get_results(compressed_file="", all_hits_dir=""):
    gz_file = compressed_file
    seq_ID = gz_file[:-3]

    os.system(f"tar -xvf {gz_file}")

    # Remove unnecesarry files
    os.system(f"rm {gz_file}")
    os.system("rm alis_afdb-swissprot_report.m8")
    os.system("rm alis_afdb50_report.m8")
    os.system("rm alis_afdb-proteome_report.m8")

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

    # Write per-structure all-hits TSV
    all_hits_tsv = "Database\tAnnotation\n"
    for hit in swiss_hits:
        all_hits_tsv += "Swiss-Prot\t" + hit + "\n"
    for hit in afdb50_hits:
        all_hits_tsv += "UniProt\t" + hit + "\n"
    for hit in proteome_hits:
        all_hits_tsv += "AlphaFold-Proteomes\t" + hit + "\n"

    with open(os.path.join(all_hits_dir, f"{seq_ID}_all_hits.tsv"), "w") as f:
        f.write(all_hits_tsv)

    # Return best-hits row for the combined table
    return [seq_ID, best_swissprot_hit, best_afdb50_hit, best_proteomes_hit,
            most_common_hit, most_common_hit_n]


# Get all CIF files in the current directory
cif_files = glob.glob("./done/*.cif") # pdb fromat works too

# Set up output directories
results_dir = "./results"
all_hits_dir = os.path.join(results_dir, "all_hits")
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
if not os.path.exists(all_hits_dir):
    os.makedirs(all_hits_dir)

# Combined best-hits TSV
best_hits_path = os.path.join(results_dir, "best_hits.tsv")
with open(best_hits_path, "w") as f:
    f.write("Protein ID\tSwiss-Prot\tUniProt\tAlphaFold-Proteomes\tMost frequent hit\tMost frequent hit n\n")

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
                    "database[]": ["afdb50", "afdb-swissprot", "afdb-proteome"]
                }
            )

        if response.status_code == 200:
            result = response.json()
            ticket_id = result["id"]
            print(f"    Ticket ID: {ticket_id}")

            # Polling for status
            status_url = f"https://search.foldseek.com/api/ticket/{ticket_id}"
            while True:
                try:
                    status_response = requests.get(status_url)
                    status_data = status_response.json()
                    status = status_data.get("status", "UNKNOWN")
                    print(f"    Status: {status}")

                    if status == "COMPLETE":
                        print("Job complete.")
                        os.system(f"curl -L https://search.foldseek.com/api/result/download/{ticket_id} -o {seq_ID}.gz")
                        row = get_results(f"{seq_ID}.gz", all_hits_dir)
                        # Append to combined best-hits TSV
                        with open(best_hits_path, "a") as f:
                            f.write("\t".join(row) + "\n")
                        break
                    elif status == "ERROR":
                        print("Error occurred.")
                        break
                    elif status == "UNKNOWN":
                        print("Unknown status, will retry.")
                except Exception as e:
                    print(f"    Error while checking status: {e}")
                    print("    Will retry after 10 seconds...")

                time.sleep(10)  # wait before checking again to not spam the server
        elif response.status_code == 429:
            time.sleep(120)
        else:
            print(f"Failed to upload {file}, status code: {response.status_code}")

    time.sleep(20)

print(f"\nResults saved to {results_dir}/")
print(f"  Best hits: {best_hits_path}")
print(f"  All hits per structure: {all_hits_dir}/")
