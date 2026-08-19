import requests
import pandas as pd
import pullstructs
from io import StringIO
import subprocess

def add_structures(metadata):

    for idx, row in metadata.iterrows():

        if row.experimental_structure:
            structure_file = (
                pullstructs.download_pdb(row.pdb_id)
            )
        else:
            structure_file = (
                pullstructs.download_alphafold(row.uniprot_id)
            )

        metadata.loc[idx, "structure_file"] = structure_file
        print(f"Downloaded {row.uniprot_id}")

    return metadata


def write_fasta(metadata, outfile = "sequences.fasta"):

    with open(outfile, "w") as f:

        for row in metadata.itertuples():

            f.write(
                f">{row.uniprot_id}\n"
            )
            f.write(
                f"{row.sequence}\n"
            )
    return outfile

def run_mmseqs(fasta, identity):
    #assumes mmseq is set up in path, need to figure this out before running the pipeline
    cmd = [
        "mmseqs",
        "easy-cluster",
        fasta,
        "clusterDB",
        "tmp",
        "--min-seq-id",
        str(identity)
    ]
    subprocess.run(cmd, check = True)

    


def remove_redundancy(metadata, identity, tmp_dir = "tmp"):
    fasta = write_fasta(metadata)

    run_mmseqs(fasta, identity)

    clusters = pd.read_csv(
            'clusterDB_cluster.tsv',
            sep = "\t",
            header = None,
            names = [
                "cluster_rep",
                'member'
            ]
        )
    metadata = metadata.merge(
        clusters,
        left_on = 'uniprot_id',
        right_on = 'member'
    )

    metadata['is_representative'] = (metadata['uniprot_id'] == metadata['cluster_rep'])

    return metadata

def get_pdb_ids(uniprot_id):

    url = (
        f"https://rest.uniprot.org/"
        f"uniprotkb/{uniprot_id}.json"
    )

    data = requests.get(url).json()

    pdb_ids = []

    for ref in data["uniProtKBCrossReferences"]:

        if ref["database"] == "PDB":
            pdb_ids.append(ref["id"])

    if len(pdb_ids) == 0:
        return None

    return pdb_ids[0]
    
    
def uniprot_search(terms,identity = 0.70, operator = "AND"):
    query = f" {operator} ".join(terms)
    url = (
    "https://rest.uniprot.org/uniprotkb/search"
    f"?query={query}"
    "&format=tsv"
    "&fields=accession,protein_name,organism_name,length,sequence,family"
    "&size=500"
    )

    response = requests.get(url)

    results = pd.read_csv(StringIO(response.text), sep ="\t")
    print(results.columns)

    metadata_rows = []

    for row in results.itertuples():
        uniprot_id = row.Entry
        pdb_id = get_pdb_ids(uniprot_id)
        if pdb_id is None:
            experimental_structure = False
        else:
            experimental_structure = True
        metadata_rows.append({
            "uniprot_id" : uniprot_id,
            "pdb_id" : pdb_id,
            "experimental_structure": experimental_structure,
            "protein_name" : row._2,
            "organism" : row.Organism,
            "length" : row.Length,
            "sequence" : row.Sequence,
            "query_source" : query,
            "structure_file" : None,
            "Family" : row.family
        })

    metadata = pd.DataFrame(metadata_rows)
    metadata = remove_redundancy(metadata, identity = identity)

    #removing redundancy from metadata
    clean_metadata = metadata[metadata['is_representative']].copy()

    #start pulling structures. If no crystal structure exists, pull the AF structure
    clean_metadata = add_structures(clean_metadata)

    return clean_metadata



   

    



    