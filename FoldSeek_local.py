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
                    target = line[1].split()[0]
                    # Header looks like "sp|P12345|NAME ...": take the accession
                    uniprot_id = target.split("|")[1] if "|" in target else target
                    annotation = " ".join(line[1].split()[1:])
                    probability = float(line[10])

                    if not "uncharacterized" in annotation.lower() and not "hypothetical" in annotation.lower() and not "putative" in annotation.lower() and not "predicted" in annotation.lower():
                        if probability >= 0.5 and best_hit == "None":
                            best_hit = annotation
                            all_hits.append((uniprot_id, annotation))
                        elif probability >= 0.5:
                            all_hits.append((uniprot_id, annotation))

        return best_hit, all_hits

    best_swissprot_hit, swiss_hits = check_db_results(file="alis_afdb-swissprot.m8")
    best_afdb50_hit, afdb50_hits = check_db_results(file="alis_afdb50.m8")
    best_proteomes_hit, proteome_hits = check_db_results(file="alis_afdb-proteome.m8")

    all_hits = swiss_hits + afdb50_hits + proteome_hits

    if len(all_hits) > 0:
        # Most common element finder (by annotation, ignoring the ID)
        annotations = [a for _, a in all_hits]
        most_common_hit = max(set(annotations), key=annotations.count)
        most_common_hit_n = str(annotations.count(most_common_hit))
    else:
        most_common_hit = "None"
        most_common_hit_n = "0"

    os.system("rm alis_afdb50.m8")
    os.system("rm alis_afdb-swissprot.m8")
    os.system(f"rm alis_afdb-proteome.m8")

    # Write per-structure all-hits TSV
    all_hits_tsv = "Database\tUniProt ID\tAnnotation\n"
    for uid, hit in swiss_hits:
        all_hits_tsv += "Swiss-Prot\t" + uid + "\t" + hit + "\n"
    for uid, hit in afdb50_hits:
        all_hits_tsv += "UniProt\t" + uid + "\t" + hit + "\n"
    for uid, hit in proteome_hits:
        all_hits_tsv += "AlphaFold-Proteomes\t" + uid + "\t" + hit + "\n"

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

# Skip already predicted structures
already_done = {f.split("/")[-1].replace("_all_hits.tsv", "") for f in glob.glob(os.path.join(all_hits_dir, "*_all_hits.tsv"))}
cif_files = [f for f in cif_files if f.split("/")[-1][:-4] not in already_done]
if already_done:
    print(f"Skipping {len(already_done)} already predicted structure(s).")
print(f"{len(cif_files)} structure(s) to process.")

# Combined best-hits TSV
best_hits_path = os.path.join(results_dir, "best_hits.tsv")
if not os.path.exists(best_hits_path):
    with open(best_hits_path, "w") as f:
        f.write("Protein ID\tSwiss-Prot\tUniProt\tAlphaFold-Proteomes\tMost frequent hit\tMost frequent hit n\n")

# Process files in batches
batch_size = 5
for i in range(0, len(cif_files), batch_size):
    batch = cif_files[i:i + batch_size]
    print(f"\nBatch {i // batch_size + 1}: uploading {len(batch)} structure(s)...")

    # Step 1: Submit all files in the batch
    pending = {}  # ticket_id -> seq_ID
    for file in batch:
        seq_ID = file.split("/")[-1][:-4]
        print(f"  Uploading: {seq_ID}")
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
            ticket_id = response.json()["id"]
            pending[ticket_id] = seq_ID
            print(f"    Ticket: {ticket_id}")
        elif response.status_code == 429:
            print("  Rate limited, waiting 120s...")
            time.sleep(120)
            # Retry this file
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
                ticket_id = response.json()["id"]
                pending[ticket_id] = seq_ID
                print(f"    Ticket: {ticket_id}")
            else:
                print(f"  Failed to upload {seq_ID}, status code: {response.status_code}")
        else:
            print(f"  Failed to upload {seq_ID}, status code: {response.status_code}")

    # Step 2: Poll all pending tickets until all are done
    while pending:
        time.sleep(10)
        done_tickets = []
        for ticket_id, seq_ID in pending.items():
            try:
                status_response = requests.get(f"https://search.foldseek.com/api/ticket/{ticket_id}")
                status = status_response.json().get("status", "UNKNOWN")

                if status == "COMPLETE":
                    print(f"  {seq_ID}: complete, downloading results...")
                    os.system(f"curl -sL https://search.foldseek.com/api/result/download/{ticket_id} -o {seq_ID}.gz")
                    row = get_results(f"{seq_ID}.gz", all_hits_dir)
                    with open(best_hits_path, "a") as f:
                        f.write("\t".join(row) + "\n")
                    done_tickets.append(ticket_id)
                elif status == "ERROR":
                    print(f"  {seq_ID}: server error.")
                    done_tickets.append(ticket_id)
            except Exception as e:
                print(f"  {seq_ID}: poll error ({e}), will retry...")

        for ticket_id in done_tickets:
            del pending[ticket_id]

        if pending:
            print(f"  Waiting on {len(pending)} job(s)...")

    time.sleep(5)  # brief pause between batches

print(f"\nResults saved to {results_dir}/")
print(f"  Best hits: {best_hits_path}")
print(f"  All hits per structure: {all_hits_dir}/")
