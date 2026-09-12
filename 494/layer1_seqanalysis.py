from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd

#stuff to individually analyze an inputted sequence
def residue_fraction(sequence, residues):
    count = sum(
        sequence.count(r)
        for r in residues
    )
    return count / len(sequence)

def max_basic_run(sequence):
    longest = 0
    current = 0
    basic = {"K", "R", "H"}
    for aa in sequence:
        if aa in basic:
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return longest

def longest_disorder(scores):
    longest = 0
    current = 0

    for score in scores:
        if score > 0.5:
            current += 1
            longest = max(longest, current)
        else:
            current = 0

    return longest

def analyze_sequence(sequence):
    sequence = str(sequence)
    sequence = sequence.strip()
    sequence = sequence.replace(" ", "")
    sequence = sequence.replace("\n", "")
    sequence = sequence.replace("\r", "")
    sequence = sequence.upper()

    

    analysis = ProteinAnalysis(sequence)

    feature_vector = {
        "length" : len(sequence),
        "MW" : analysis.molecular_weight(),
        "PI" : analysis.isoelectric_point(),
        "aromaticity" : analysis.aromaticity(),
        "instability_index" : analysis.instability_index(),
        "gravy" : analysis.gravy(),
        "charge_at_7pH" : analysis.charge_at_pH(7.0),
    }

    aa_frequencies = analysis.amino_acids_percent
    for aa in aa_frequencies.keys():
        feature_vector[f"aa_{aa}"] = aa_frequencies[aa]

    positive = "KRH"
    negative = "DE"
    polar = "STNQ"
    hydrophobic = "AVILMFWY"

    feature_vector['positive_aa_freq'] = residue_fraction(sequence, positive)
    feature_vector['negative_aa_freq'] = residue_fraction(sequence, negative)
    feature_vector['polar_aa_freq'] = residue_fraction(sequence, polar)
    feature_vector['hydrophobic_aa_freq'] = residue_fraction(sequence, hydrophobic)

    feature_vector['LTA_motif_freq'] = (sequence.count("KK") + sequence.count("KR") + sequence.count("RK") + sequence.count("RR")) / feature_vector['length']
    feature_vector['max_basic_run'] = max_basic_run(sequence)

    feature_vector['gly_freq'] = sequence.count("G") / feature_vector['length']
    feature_vector['pro_freq'] = sequence.count("P") / feature_vector['length']
    feature_vector['cys_freq'] = sequence.count("C") / feature_vector['length']
    feature_vector['aromatic_freq'] = (sequence.count("F") + sequence.count("Y") + sequence.count("W")) / feature_vector['length']

    '''#clone IUpred 3 repository and run it locally, get the disorder score per residue, put it in scores variable
    scores = ""
    feature_vector['fraction_disordered'] = sum(score > 0.5 for score in scores) / feature_vector['length']
    feature_vector['mean_disorder'] = sum(scores) / feature_vector['length']
    feature_vector['longest_disorder'] = longest_disorder(scores)'''

    return feature_vector

#working with entire inputted dataframe of sequences, assuming sequence column name is "sequence"

def apply_layer1(metadata):
    working_df = metadata.copy()

    features = []

    for row in working_df.itertuples():
        features.append(analyze_sequence(row.sequence))

    feature_df = pd.DataFrame(features)

    concat = pd.concat([working_df.reset_index(drop = True), feature_df], axis = 1)

    with pd.ExcelWriter("layer1_analysis.xlsx") as writer:
        concat.to_excel(writer, sheet_name = "1", index = None)


    return concat


df = pd.read_excel("initial.xlsx")

apply_layer1(df)


