'''organization of the active site, try to characterize the binding pocket of each structure, look into if theres possible structures that would
screw up the data and shouldnt be analyzed, use fpocket, also look into other programs that can be used to get different types of info

'''
import subprocess
import os
import pandas as pd
import pullstructs

# ---------------------------------------------------------
# Run fpocket inside WSL with conda activation
# ---------------------------------------------------------
def run_fpocket(pdb_file):
    print(f"Running fpocket on {pdb_file}")

    # Directory and basename for the PDB
    pdb_dir = os.path.dirname(pdb_file)
    pdb_name = os.path.basename(pdb_file)

    # Remove old fpocket output in that directory (if any)
    out_dir = os.path.join(pdb_dir, pdb_name.replace(".pdb", "_out"))
    if os.path.exists(out_dir):
        subprocess.run(["wsl", "bash", "-c", f"cd {pdb_dir} && rm -rf {pdb_name.replace('.pdb', '_out')}"], check=True)

    # Activate conda, cd into the PDB directory, and run fpocket on the basename
    cmd = f"source ~/miniconda3/bin/activate && cd {pdb_dir} && fpocket -f {pdb_name}"
    subprocess.run(["wsl", "bash", "-c", cmd], check=True)

    # Now fpocket will have written <basename>_out in the same directory as the PDB
    info_file = os.path.join(out_dir, pdb_name.replace(".pdb", "_info.txt"))

    if not os.path.exists(info_file):
        print("No fpocket output found — returning empty pocket list.")
        return []

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




# ---------------------------------------------------------
# Extract fpocket features
# ---------------------------------------------------------
def fpocket_features(pockets):
    if len(pockets) == 0:
        return {
            "num_pockets": 0,
            "max_druggability": 0,
            "largest_pocket_volume": 0,
            "primary_pocket_volume": 0,
            "pocket_with_highest_hydrophobicity": None,
            "pocket_with_highest_polarity": None
        }

    # Convert numeric fields safely
    for p in pockets:
        for key in ["Volume", "Hydrophobicity Score", "Polarity Score", "Druggability Score"]:
            if key in p:
                try:
                    p[key] = float(p[key])
                except:
                    p[key] = 0.0

    vols = [p.get("Volume", 0.0) for p in pockets]
    drugs = [p.get("Druggability Score", 0.0) for p in pockets]

    # Primary pocket = highest druggability
    primary = max(pockets, key=lambda p: p.get("Druggability Score", 0.0))

    # Pocket with highest hydrophobicity
    hydrophobic_pocket = max(pockets, key=lambda p: p.get("Hydrophobicity Score", 0.0))

    # Pocket with highest polarity
    polar_pocket = max(pockets, key=lambda p: p.get("Polarity Score", 0.0))

    return {
        "num_pockets": len(pockets),
        "max_druggability": max(drugs),
        "largest_pocket_volume": max(vols),

        # New attributes
        "primary_pocket_volume": primary.get("Volume", 0.0),
        "pocket_with_highest_hydrophobicity": hydrophobic_pocket.get("id", None),
        "pocket_with_highest_polarity": polar_pocket.get("id", None)
    }


# ---------------------------------------------------------
# Apply Layer 4 to dataframe
# ---------------------------------------------------------
def apply_layer4(metadata):

    working_df = metadata.copy()
    features = []

    for row in working_df.itertuples():
        # Build WSL path
        if row.pull_id.startswith("AF-"):
            real_id = row.pull_id[3:-3]
            pdb_path = f"C:/Users/willm/OneDrive/Desktop/FunctionFromSeq/494/structures/{real_id}.pdb"
        else:
            pdb_path = f"C:/Users/willm/OneDrive/Desktop/FunctionFromSeq/494/structures/{row.pull_id}.pdb"

        wsl_path = pdb_path.replace("C:/", "/mnt/c/").replace("\\", "/")
        
        pockets = run_fpocket(wsl_path)
        features.append(fpocket_features(pockets))

    feature_df = pd.DataFrame(features)
    concat = pd.concat([working_df.reset_index(drop=True), feature_df], axis=1)

    with pd.ExcelWriter("layer4_analysis.xlsx") as writer:
        concat.to_excel(writer, sheet_name='1', index=None)

    pullstructs.clean_structure_directory("C:/Users/willm/OneDrive/Desktop/FunctionFromSeq/494/structures")

    return concat


# ---------------------------------------------------------
# Run Layer 4
# ---------------------------------------------------------
df = pd.read_excel("initial.xlsx")
apply_layer4(df)
