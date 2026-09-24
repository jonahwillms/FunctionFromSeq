'''organization of the active site, try to characterize the binding pocket of each structure, look into if theres possible structures that would
screw up the data and shouldnt be analyzed, use fpocket, also look into other programs that can be used to get different types of info

'''
import subprocess
import os
import pandas as pd

def run_fpocket(pdb_file):
    print(f"Running fpocket on {pdb_file}")

    subprocess.run(["fpocket", "-f", pdb_file], check=True)

    out_dir = pdb_file.replace(".pdb","_out")
    info_file = os.path.join(out_dir, pdb_file.replace(".pdb", "_info.txt"))

    pockets = []
    with open(info_file) as f:
        current = {}
        for line in f:
            if line.startswith("Pocket"):
                if current:
                    pockets.append(current)
                current = {"id": line.strip()}
            elif ":" in line:
                key, val = line.split(":", 1)
                current[key.strip()] = val.strip()
        if current:
            pockets.append(current)

        return pockets

def fpocket_features(pockets):
    if len(pockets) == 0:
        return {
            "num_pockets" : 0,
            "max_druggability" : 0,
            "largest_pocket_volume" : 0,
            "pocket_with_highest_hydrophobicity" : None,
            "pocket_with_highest_polarity" : None
        }

    for p in pockets:
        for key in ["Volume", "Hydrophobicity Score", "Polarity Score", "Druggability Score"]:
            if key in p:
                try:
                    p[key] = float(p[key])
                except:
                    p[key] = 0.0

    vols = [float(p['Volume1']) for p in pockets]
    drugs = [float(p["Druggability Score"]) for p in pockets]

    primary = max(pockets, key = lambda p: p.get("Druggability Score", 0.0))

    hydrophobic_pocket = max(pockets, key = lambda p: p.get("Hydrophobicity Score", 0.0))

    polar_pocket = max(pockets, key = lambda p: p.get("Polarity Score", 0.0))

    return {
        "num_pockets": len(pockets),
        "max_druggability": max(drugs),
        "largest_pocket_volume": max(vols),

        # New attributes
        "primary_pocket_volume": primary.get("Volume", 0.0),
        "pocket_with_highest_hydrophobicity": hydrophobic_pocket.get("id", None),
        "pocket_with_highest_polarity": polar_pocket.get("id", None)
    }


def apply_layer4(metadata):
    working_df = metadata.copy()
    features = []

    for row in working_df.itertuples():
        id = ""
        if row.pull_id[:3] == "AF-":
            id = row.pull_id[3:-3]
            id = f"structures/{id}.pdb"
        else:
            id = f"structures/{row.pull_id}.pdb"

        features.append(fpocket_features(run_fpocket(id)))

    feature_df = pd.DataFrame(features)

    concat = pd.concat([working_df.reset_index(drop = True), feature_df], axis = 1)

    with pd.ExcelWriter("layer4_analysis.xlsx") as writer:
        concat.to_excel(writer, sheet_name = '1', index = None)

    return concat