from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import MDAnalysis as mda
import numpy as np



def analyze_structure(pdb_file):
    u = mda.Universe(pdb_file)

    atoms = u.select_atoms("protein")

    feature_vector = {}

    #radius of gyration
    feature_vector['radius_of_gyration'] = (atoms.radius_of_gyration())

    coords = atoms.positions

    #centering coordinates
    coords = coords - coords.mean(axis = 0)

    cov = np.cov(coords.T)

    eigenvalues, eigenvectors = np.linalg.eigh(cov)

    eigenvalues = np.sort(eigenvalues)[::-1]

    feature_vector['principal_axis_1'] = (eigenvalues[0])
    feature_vector['principal_axis_2'] = (eigenvalues[1])
    feature_vector['principal_axis_3'] = (eigenvalues[2])

    feature_vector['axis_ratio_1'] = eigenvalues[0] / eigenvalues[1]
    feature_vector['axis_ratio_2'] = eigenvalues[1] / eigenvalues[2]

    #moment of inertia

    inertia = atoms.moment_of_inertia()

    inertia_eigs = np.linalg.eigvalsh(inertia)
    inertia_eigs = np.sort(inertia_eigs)[::-1]

    feature_vector['moment_1'] = (inertia_eigs[0])
    feature_vector['moment_2'] = (inertia_eigs[1])
    feature_vector['moment_3'] = (inertia_eigs[2])
    feature_vector['moment_ratio_1'] = (inertia_eigs[0] / inertia_eigs[1])
    feature_vector['moment_ratio_2'] = (inertia_eigs[1]) / inertia_eigs[2]

    return feature_vector

def apply_layer2(metadata):
    working_df = metadata.copy()
    features = []

    for row in working_df.itertuples():
        id = ""
        if row.pull_id[:3] == "AF-":
            id = row.pull_id[3:-3]
            id = f"{id}.pdb"
        else:
            id = f"{row.pull_id}.pdb"

        features.append(analyze_structure(id))

    feature_df = pd.DataFrame(features)

    concat = pd.concat([working_df.reset_index(drop = True), feature_df], axis = 1)

    with pd.ExcelWriter("layer2_analysis.xlsx") as writer:
        concat.to_excel(writer, sheet_name = '1', index = None)

    return concat
    