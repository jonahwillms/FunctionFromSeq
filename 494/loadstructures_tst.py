import pandas as pd
import layer2_structgeometry
import pullstructs
from Bio.PDB import PDBIO, PDBParser, Select
import shutil
#function to download and clean all of the pdb structures, can pull all the structures with pullstructs, then I need to clean them, then comment out the function,
#then run layer 2


def structure_setup():

    #downloading structures
    

    #cleaning structures
    def clean_pdb(pdb_file, expected_length):

        if pdb_file[:3] == "AF-":
            pdb_file = pdb_file[3:-3]

        pdb_file = f'{pdb_file}.pdb'

        def get_chain_length(pdb_file):
            parser = PDBParser(QUIET = True)

            structure = parser.get_structure("protein", pdb_file)

            model = structure[0]

            lengths = {}

            for chain in model:
                lengths[chain.id] = sum(residue.id[0] == " " for residue in chain)

            return lengths
        lengths = get_chain_length(pdb_file)

        if len(lengths) == 1:
            return pdb_file
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
        io.save(pdb_file, ChainSelect(rep))
        return {
        "selected_chain": rep,
        "chain_lengths": lengths
            }

    df = pd.read_excel("initial.xlsx")

    for row in df.itertuples():
        clean_pdb(row.pull_id, len(row.sequence))


structure_setup()

'''df = pd.read_excel("initial.xlsx")

layer2_structgeometry.apply_layer2(df)'''