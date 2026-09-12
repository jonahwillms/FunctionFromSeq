from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import dssp
import freesasa

def analyze_structure(pdb_file):

    u = mda.Universe(pdb_file)

    print(f"Now analyzing {pdb_file}")

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

    #computing SASA
    structure = freesasa.Structure(pdb_file)
    result = freesasa.calc(structure)

    total_sasa = result.totalArea()
    feature_vector['total_sasa'] = float(total_sasa)


    #DSSP analysis
    protein = u.select_atoms("protein")

    #alt locations in pdb files screw up analysis, need to remove
    protein = protein.select_atoms("(not altloc B)")

    valid_atoms = []
    for res in protein.residues:
        bb = res.atoms.select_atoms("name N CA C O")
        if len(bb) == 4: #backbone is complete for residue
            valid_atoms.extend(bb)

    if len(valid_atoms) == 0:
        feature_vector['helix_fraction'] = np.nan
        feature_vector['sheet_fraction'] = np.nan
        feature_vector['coil_fraction'] = np.nan
        return feature_vector

    backbone = mda.core.groups.AtomGroup(valid_atoms)
    d = dssp.DSSP(backbone)
    d.run()

    ss = d.results['dssp'].flatten()
    ss_str = "".join(ss)
    #print(ss_str)
    
    feature_vector['helix_fraction'] = ss_str.count("H") / len(ss_str)
    feature_vector['sheet_fraction'] = ss_str.count("E") / len(ss_str)
    feature_vector['coil_fraction'] = ss_str.count("-") / len(ss_str)

    return feature_vector

def apply_layer2(metadata):
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

    with pd.ExcelWriter("layer2_analysis.xlsx") as writer:
        concat.to_excel(writer, sheet_name = '1', index = None)

    return concat



df = pd.read_excel("initial.xlsx")
apply_layer2(df)
    