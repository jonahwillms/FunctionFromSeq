import pandas as pd
import plotly.express as px

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def analyze_df():

    df = pd.read_excel("layer2_analysis.xlsx")
    df1 = pd.read_excel("layer1_analysis.xlsx")

    df = df.merge(df1, on = ["Source", "Function", "Class", "sequence"], how = 'left')


    plot_info = df[
        ["Source", "Function", "Class"]
    ].copy()


    '''feature_cols = [
    "radius_of_gyration",
    "principal_axis_1",
    "principal_axis_2",
    "principal_axis_3",
    "axis_ratio_1",
    "axis_ratio_2",
    "moment_1",
    "moment_2",
    "moment_3",
    "moment_ratio_1",
    "moment_ratio_2"
        ]'''
    feature_cols = [

    # ===== Layer 1: Sequence Features =====

    "length",
    "MW",
    "PI",
    "aromaticity",
    "instability_index",
    "gravy",
    "charge_at_7pH",

    "positive_aa_freq",
    "negative_aa_freq",
    "polar_aa_freq",
    "hydrophobic_aa_freq",

    "LTA_motif_freq",
    "max_basic_run",

    "aromatic_freq",

    "aa_A",
    "aa_C",
    "aa_D",
    "aa_E",
    "aa_F",
    "aa_G",
    "aa_H",
    "aa_I",
    "aa_K",
    "aa_L",
    "aa_M",
    "aa_N",
    "aa_P",
    "aa_Q",
    "aa_R",
    "aa_S",
    "aa_T",
    "aa_V",
    "aa_W",
    "aa_Y",

    # ===== Layer 2: Structure Geometry =====

    "radius_of_gyration",

    "principal_axis_1",
    "principal_axis_2",
    "principal_axis_3",

    "axis_ratio_1",
    "axis_ratio_2",

    "moment_1",
    "moment_2",
    "moment_3",

    "moment_ratio_1",
    "moment_ratio_2",

    # ===== Layer 2: Secondary Structure =====

    
        ]   


    X = df[feature_cols]

    scaler = StandardScaler()

    X_scaled = scaler.fit_transform(X)

    pca = PCA(n_components=3)

    X_pca = pca.fit_transform(X_scaled)

    pca_df = pd.DataFrame(
        X_pca,
        columns=["PC1", "PC2", "PC3"]
    )

    pca_df = pd.concat(
        [pca_df, plot_info.reset_index(drop=True)],
        axis=1
    )

    # custom colouring
    pca_df["Colour"] = "Other"

    pca_df.loc[
        pca_df["Class"] == "CAP/SCP",
        "Colour"
    ] = "CAP/SCP"

    pca_df.loc[
        pca_df["Class"] == "Polymer-binding",
        "Colour"
    ] = "Polymer-binding"

    pca_df.loc[
        pca_df["Class"] == "Hydrolase",
        "Colour"
    ] = "Hydrolase"

    pca_df.loc[
        pca_df["Source"] == "O31398",
        "Colour"
    ] = "YkwD"

    pca_df.loc[
        pca_df["Source"] == "H7C6X6",
        "Colour"
    ] = "SalB"

    fig = px.scatter_3d(
        pca_df,
        x="PC1",
        y="PC2",
        z="PC3",
        color="Colour",
        hover_data=[
            "Source",
            "Function",
            "Class"
        ],
        color_discrete_map={
            "YkwD": "green",
            "SalB" : "orange",
            "CAP/SCP": "red",
            "Polymer-binding": "blue",
            "Hydrolase": "purple",
            "Other": "gray"
        }
    )

    fig.update_traces(
        marker=dict(size=6)
    )

    fig.show(renderer = "browser")

    print(
        "Explained variance:"
    )

    print(
        pca.explained_variance_ratio_
    )
    loadings = pd.DataFrame(
    pca.components_.T,
    columns=["PC1", "PC2", "PC3"],
    index=feature_cols
)

    print("\nPCA Loadings:")
    print(loadings)
    

    geom_df = df.copy()

    geom_df["Colour"] = "Other"

    geom_df.loc[
        geom_df["Class"] == "CAP/SCP",
        "Colour"
    ] = "CAP/SCP"

    geom_df.loc[
        geom_df["Class"] == "Polymer-binding",
        "Colour"
    ] = "Polymer-binding"

    geom_df.loc[
        geom_df["Class"] == "Hydrolase",
        "Colour"
    ] = "Hydrolase"

    geom_df.loc[
        geom_df["Source"] == "O31398",
        "Colour"
    ] = "YkwD"

    geom_df.loc[
            pca_df["Source"] == "H7C6X6",
            "Colour"
        ] = "SalB"

    fig = px.scatter_3d(
        geom_df,
        x="axis_ratio_1",
        y="axis_ratio_2",
        z="moment_ratio_1",
        color="Colour",
        hover_data=[
            "Source",
            "Function",
            "Class"
        ],
        color_discrete_map={
                    "YkwD": "green",
                    "SalB" : "orange",
                    "CAP/SCP": "red",
                    "Polymer-binding": "blue",
                    "Hydrolase": "purple",
                    "Other": "gray"
                }
    )

    fig.update_traces(
        marker=dict(size=6)
    )

    fig.show(renderer="browser")

    return pca_df, pca
    
analyze_df()