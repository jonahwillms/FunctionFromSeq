from pathlib import Path
import requests
import pandas as pd

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

def download_alphafold(uniprot_id, out_dir = "structures"):
    url = (
        f"https://alphafold.ebi.ac.uk/files/"
        f"AF-{uniprot_id}-F1-model_v4.pdb"
    )
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

    df = input("input path of dataframe: ")
    for row in df.itertuples():

        if row.struct_type == 1:
            download_pdb(row.pull_id)
        elif row.struct_type == 0:
            download_alphafold(row.pull_id)
        t = row.ID

        print(f"Succesffuly downloaded {t}")


get_structures()