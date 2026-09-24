import requests
import pandas as pd
import pullstructs
from io import StringIO
import subprocess
from Bio.PDB import PDBParser, PDBIO, Select
import shutil

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

def get_chain_length(pdb_file):
    parser = PDBParser(QUIET = True)

    structure = parser.get_structure("protein", pdb_file)

    model = structure[0]

    lengths = {}

    for chain in model:
        lengths[chain.id] = sum(residue.id[0] == " " for residue in chain)

    return lengths

def clean_pdb(pdb_file, expected_length, output_file):

    lengths = get_chain_length(pdb_file)

    if len(lengths) == 1:
        shutil.copy(pdb_file, output_file)

        return output_file
    else:
        rep = None

        for chain, length in lengths.items():

            if (abs(length - expected_length) / expected_length) < 0.15:
                rep = chain
                break


    if rep == None:
        return None

    class ChainSelect(Select):

        def __init__(self, chain_id):
            self.chain_id = chain_id

        def accept_chain(self, chain):
            return (chain.id == self.chain_id)

    parser = PDBParser(QUIET = True)
    structure = parser.get_structure("protein", pdb_file)
    io = PDBIO()
    io.set_structure(structure)
    print(
    f"{pdb_file}: "
    f"selected chain {rep}"
    )
    io.save(output_file, ChainSelect(rep))
    return {
    "selected_chain": rep,
    "chain_lengths": lengths
        }

def clean_structure_folder(metadata):

    metadata['clean_structure_file'] = None

    for idx, row in metadata.iterrows():

        infile = row.structure_file

        outfile = (
            f"cleaned_structures/"
            f"clean_{row.structure_file}"
        )

        clean_pdb(infile, row.length, outfile)

        metadata.loc[idx, "clean_structure_file"] = outfile

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

    #start working on this next, creating the redundancy net
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

    
def choose_best_rep(df):
    #1 prefer experimentally determined pdb structures
    exp = df[df.experimental_structure]
    if len(exp) > 0:
        df = exp

    #2 longest sequence
    df = df.sort_values("length", ascending = False)

    #3 reviewed entries
    df['is_reviewed'] = ~df.uniprot_id.str.startswith("A0A")
    df = df.sort_values("is_reviewed", ascending = False)

    #if none of these things are true, just pick top uniprot id
    df = df.sort_values("uniprot_id")

    return df.iloc[0].uniprot_id
    

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

    new_reps = (metadata.groupby("cluster_rep").apply(choose_best_rep).reset_index(name = "new_rep"))

    metadata = metadata.merge(new_reps, on = "cluster_rep")
    
    metadata['is_representative'] = (metadata.uniprot_id == metadata.new_rep)

    clean_metadata = metadata[metadata.is_representative].copy()

    return clean_metadata

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
   
    #removing redundancy from metadata
    clean_metadata = remove_redundancy(metadata, identity = identity)

    #start pulling structures. If no crystal structure exists, pull the AF structure
    clean_metadata = add_structures(clean_metadata)

    #filter pdb structures so that all inputted chains are monomers (this will match better with alphafold stuff)
    clean_metadata = clean_structure_folder(clean_metadata)

    return clean_metadata



    
        
   

    



    