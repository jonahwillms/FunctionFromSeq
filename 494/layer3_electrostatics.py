from gridData import Grid
import freesasa
import subprocess
import numpy as np
import pandas as pd
import pullstructs

def strip_remarks(pdb_file, cleaned_file):
    skip_prefixes = (
        "HEADER", "TITLE", "COMPND", "SOURCE", "KEYWDS", "EXPDTA",
        "REVDAT", "DBREF", "SEQADV", "SEQRES", "HETNAM", "FORMUL",
        "HELIX", "SHEET", "CISPEP", "SITE", "CRYST1",
        "ORIGX1", "ORIGX2", "ORIGX3",
        "SCALE1", "SCALE2", "SCALE3",
        "REMARK", "AUTHOR", "JRNL", "HET", "CONECT", "MASTER", "TER"
    )

    with open(pdb_file) as f, open(cleaned_file, "w") as out:
        for line in f:
            if any(line.startswith(prefix) for prefix in skip_prefixes):
                continue
            out.write(line)


def analyze_structure(pdb_file):

    print(f"Layer 3 analyzing {pdb_file}")
    feature_vector = {}
    clean_pdb = pdb_file.replace(".pdb", "_clean.pdb")
    strip_remarks(pdb_file, clean_pdb)
    pqr_file = pdb_file.replace(".pdb", ".pqr")
    
    cmd = ["pdb2pqr30", "--ff=AMBER", clean_pdb, pqr_file]
    subprocess.run(cmd, check = True)

    apbs_input = pdb_file.replace(".pdb", ".in")
    dx_output = pdb_file.replace(".pdb", "")

    with open(apbs_input, "w") as f:
        f.write(f"""
read
    mol pqr {pqr_file}
end

elec
    mg-auto
    dime 161 161 161
    cglen 150 150 150
    fglen 120 120 120
    cgcent mol 1
    fgcent mol 1
    mol 1
    lpbe
    bcfl sdh
    pdie 2.0
    sdie 78.0
    chgm spl2
    srfm smol
    sdens 10.0
    srad 1.4
    temp 298.15
    write pot dx {dx_output}
end

quit
""")
    subprocess.run(["apbs", apbs_input], check = True)
    grid = Grid(dx_output + ".dx")
    potential_map = grid.grid

    feature_vector['mean_potential'] = float(np.mean(potential_map))
    feature_vector['std_potential'] = float(np.std(potential_map))

    net_charge = 0.0
    with open(pqr_file) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                parts = line.split()
                charge = float(parts[-2])
                net_charge += charge

    feature_vector['net_charge'] = float(net_charge)

    print("")

    return feature_vector

def apply_layer3(metadata):
    working_df = metadata.copy()
    features = []

    for row in working_df.itertuples():
        id = ""
        if row.pull_id[:3] == "AF-":
            id = row.pull_id[3:-3]
            id = f"structures/{id}.pdb"
        else:
            id = f"structures/{row.pull_id}.pdb"

        features.append(analyze_structure(id))

    feature_df = pd.DataFrame(features)

    concat = pd.concat([working_df.reset_index(drop = True), feature_df], axis = 1)

    with pd.ExcelWriter("layer3_analysis.xlsx") as writer:
        concat.to_excel(writer, sheet_name = '1', index = None)

    pullstructs.clean_structure_directory("C:/Users/willm/OneDrive/Desktop/FunctionFromSeq/494/structures")

    return concat

df = pd.read_excel("initial.xlsx")
apply_layer3(df)