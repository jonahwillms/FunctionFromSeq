from pathlib import Path
import requests
import pandas as pd
import subprocess
import os

def download_pdb(pdb_id, out_dir = "structures"):
    pdb_id = pdb_id.upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    Path(out_dir).mkdir(exist_ok = True)

    outfile = Path(out_dir) / f"{pdb_id}.pdb"

    r = requests.get(url)

    if r.status_code != 200:
        raise ValueError(f"Could not download PDB {pdb_id}")

    outfile.write_text(r.text)

    return str(outfile)

def get_alphafold_url(uniprot_id):

    url = (
        f"https://alphafold.ebi.ac.uk/api/prediction/"
        f"{uniprot_id}"
    )

    r = requests.get(url)

    if r.status_code != 200:
        return None

    data = r.json()

    if len(data) == 0:
        return None

    return data[0]['pdbUrl']

def download_alphafold(uniprot_id, out_dir = "structures"):
    
    url = get_alphafold_url(uniprot_id)
    print(url)
    Path(out_dir).mkdir(exist_ok = True)

    outfile = Path(out_dir) / f"{uniprot_id}.pdb"

    r = requests.get(url)

    if r.status_code != 200:
        raise ValueError(
            f"Could not download Alphafold Model "
            f"for {uniprot_id}"
        )

    outfile.write_text(r.text)

    return str(outfile)

def get_structures():

    df = pd.read_excel(input("input path of dataframe: "))
    for row in df.itertuples():

        if row.struct_type == 1:
            download_pdb(row.pull_id)
        elif row.struct_type == 0:
            download_alphafold(row.pull_id[3:-3])
        t = row.pull_id

        print(f"Succesffuly downloaded {t}")


#get_structures()

def clean_structure_directory(win_dir):
    """
    Cleans the structure directory so only .pdb files remain.
    Deletes:
      - *_clean.pdb
      - *_out/ directories
      - all non-PDB files
    """

    # Convert Windows path to WSL path
    wsl_dir = win_dir.replace("C:/", "/mnt/c/").replace("\\", "/")

    for item in os.listdir(win_dir):
        win_path = os.path.join(win_dir, item)

        # Keep ONLY .pdb files (not *_clean.pdb)
        if item.endswith(".pdb") and not item.endswith("_clean.pdb"):
            continue

        # Delete *_clean.pdb
        if item.endswith("_clean.pdb"):
            os.remove(win_path)
            continue

        # Delete fpocket output directories
        if item.endswith("_out") and os.path.isdir(win_path):
            wsl_path = win_path.replace("C:/", "/mnt/c/").replace("\\", "/")
            subprocess.run(["wsl", "rm", "-rf", wsl_path], check=True)
            continue

        # Delete all other files (APBS, DX, PQR, JSON, logs, etc.)
        if os.path.isfile(win_path):
            os.remove(win_path)
        else:
            # Delete any other directories
            wsl_path = win_path.replace("C:/", "/mnt/c/").replace("\\", "/")
            subprocess.run(["wsl", "rm", "-rf", wsl_path], check=True)