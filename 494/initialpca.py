import pandas as pd
import plotly.express as px

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def analyze_df():

    df = pd.read_excel("layer2_analysis.xlsx")

    plot_info = df[
        ["Source", "Function", "Class"]
    ].copy()


    feature_cols = [
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