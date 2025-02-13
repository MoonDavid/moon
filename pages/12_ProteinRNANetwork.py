import streamlit as st
import pandas as pd
from pyvis.network import Network
import streamlit.components.v1 as components
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def df_from_uniprots(df, uniprots):
    """Create a DataFrame from a list of UniProt IDs."""
    return df[df['UniProtKB-AC'].isin(uniprots)].reset_index(drop=True)
def analyze_network(data):
    """
    Analyze network properties and return statistics.

    :param data: DataFrame with source and target columns
    :return: Dictionary of network statistics
    """
    # Create NetworkX graph
    G = nx.from_pandas_edgelist(data, "source", "target")

    stats = {
        "Nodes": G.number_of_nodes(),
        "Edges": G.number_of_edges(),
        "Density": nx.density(G),
        "Average Degree": sum(dict(G.degree()).values()) / G.number_of_nodes(),
        "Connected Components": nx.number_connected_components(G),
    }

    # Get unique RBPs and their RNA counts
    rbps = data["source"].unique()
    rna_counts = data.groupby("source").size()

    return stats, rbps, rna_counts


def plot_score_distribution(data):
    """
    Plot the distribution of scores.

    :param data: DataFrame with a score column
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.histplot(data=data, x="score", bins=30)
    plt.title("Distribution of Interaction Scores")
    plt.xlabel("Score")
    plt.ylabel("Count")
    return fig


def plot_rna_counts(rna_counts):
    """
    Plot the number of RNAs associated with each protein.

    :param rna_counts: Series with protein-RNA counts
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    rna_counts.plot(kind="bar")
    plt.title("Number of RNAs Associated with Each Protein")
    plt.xlabel("Protein")
    plt.ylabel("Number of Associated RNAs")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    return fig

def initialize():
    st.set_page_config(
        page_title="Protein-RNA Network",
        page_icon="🧬",
        layout="wide"
    )
    st.title("Protein-RNA Network Explorer")

@st.cache_data
def load_data(file_path):
    """Load CSV data into a pandas DataFrame."""
    try:
        data = pd.read_csv(file_path)
        return data
    except Exception as e:
        st.error(f"Error loading data: {e}")
        return None


def create_network(data, physics=True):
    """
    Build a PyVis Network object from the provided DataFrame.

    :param data: DataFrame with at least 'source' and 'target' columns.
    :param physics: Boolean to enable/disable physics simulation in the network.
    :return: A PyVis Network object.
    """
    net = Network(height="750px", width="50%", notebook=False)

    # Add nodes and edges from DataFrame
    for _, row in data.iterrows():
        if 'source' in row and 'target' in row:
            net.add_node(row['source'], label=row['source'])
            net.add_node(row['target'], label=row['target'])
            net.add_edge(row['source'], row['target'])

    # Show PyVis configuration buttons (filter only for 'physics')
    net.show_buttons(filter_=['node','edge','physics'])

    # Optionally disable physics if user unchecks the box
    if not physics:
        net.set_options("""
            var options = {
                "physics": {
                    "enabled": false
                }
            }
        """)

    return net


def render_pyvis_network(nt: Network, filename="network.html"):
    """
    Render the PyVis Network in a Streamlit app and provide a download button.

    :param nt: PyVis Network object
    :param filename: Name of the HTML file to save and render.
    """
    try:
        # Save the network to an HTML file
        nt.save_graph(filename)

        # Read and display the HTML
        with open(filename, 'r', encoding='utf-8') as file:
            html_content = file.read()
        components.html(html_content, height=1600, scrolling=True)

        # Provide a download button for the HTML
        with open(filename, 'rb') as file:
            st.download_button(
                label="Download Network Visualization",
                data=file,
                file_name=filename,
                mime="text/html"
            )
    except Exception as viz_error:
        st.error(f"Visualization Error: {viz_error}")
def plot_rbp_counts_st(data):
    # Get counts
    rbp_counts = data.groupby("rbp").size().reset_index(name="count")

    # Create bar chart using Streamlit
    st.bar_chart(data=rbp_counts, x="rbp", y="count")
    return rbp_counts
def plot_rbp_counts_altair(data):
    # Get counts
    rbp_counts = data.groupby("rbp").size().reset_index(name="count")

    # Create Altair chart
    import altair as alt

    chart = (
        alt.Chart(rbp_counts)
        .mark_bar()
        .encode(
            x=alt.X("rbp", sort="-y", title="RBP"),  # Sort by count descending
            y=alt.Y("count", title="Number of Interactions"),
            tooltip=["rbp", "count"],
        )
        .properties(title="Number of Interactions per RBP", width=600, height=400)
    )

    # Display in Streamlit
    st.altair_chart(chart, use_container_width=True)
    return rbp_counts



def main():
    initialize()
    if "gene_lists" not in st.session_state:
        st.warning(
            "No pre-loaded MPs list found in session state. Please load first main page Moonlighting Proteins Analysis to add some."
        )
        st.stop()
    st.markdown("""
        This section explore MPs-RNA interactions using the [POSTAR3/ClipDB](http://111.198.139.65/RBP.htm) database.

        ### Overview
        POSTAR3 represents the most comprehensive collection of RBP-RNA interactions, featuring:
        - **1,499** CLIP-seq datasets
        - **348** RBPs (**219** human)
        - **7** species coverage 
        """)
    col1, col2 = st.columns(2)

    # Left column - Technologies
    with col1:
        st.markdown("### Technologies")
        st.markdown("""
        - HITS-CLIP
        - PAR-CLIP
        - iCLIP
        - eCLIP
        - iCLAP
        - urea-iCLIP
        - 4sU-iCLIP
        - BrdU-CLIP
        - Fr-iCLIP
        - PIP-seq
        """)

    # Right column - Data Processing
    with col2:
        st.markdown("### Data Processing")
        st.markdown("""
        - **HITS-CLIP**: CLIPper (human) / CTK (other species)
        - **PAR-CLIP**: MiClip
        - **iCLIP-related**: PureCLIP
        """)

    # Load data from CSV
    file_path = r"C:\Users\david\PycharmProjects\dataframes\clipseq.csv"
    data = load_data(file_path)
    data.rename(columns={"source": "experiment"}, inplace=True)
    data["source"] = data["rbp"]
    data["target"] = (
        data["chromosome"]
        + ":"
        + data["start"].astype(str)
        + "-"
        + data["end"].astype(str)
    )

    # Sidebar for selection and controls
    with st.sidebar:
        st.title("Controls")

        # Gene list selection
        available_gene_lists = {
            key: value
            for key, value in st.session_state["gene_lists"].items()
            if isinstance(value, list)
        }

        if not available_gene_lists:
            st.warning("No gene list available in session state.")
            return

        selected_gene_list = st.selectbox(
            "Select the gene list:",
            options=[""] + list(available_gene_lists.keys()),  # Empty option first
        )

        if not selected_gene_list:  # If no selection
            st.warning("Please select a gene list to proceed")
            return

    # Main content

    # Process selected gene list
    gene_list = available_gene_lists[selected_gene_list]
    entrez = df_from_uniprots(st.session_state["df"], gene_list)["Entrez"].to_list()
    intersection = set(entrez).intersection(set(data["source"]))

    # Filter data
    if "SMN1" in entrez or "SMN2" in entrez:
        filtered_data = data[data["source"].isin(list(intersection) + ["SMN"])]
    else:
        filtered_data = data[data["source"].isin(list(intersection))]

    # Data Overview
    st.header("Data Overview")
    st.write(
        f"ClipDB intersection with {selected_gene_list} has {filtered_data.shape[0]} interactions "
        f"and {filtered_data['source'].nunique()} unique RBPs."
    )

    with st.expander("Show Raw Data", expanded=False):
        st.dataframe(filtered_data)

    # Score Analysis
    st.header("Score Analysis")

    # Calculate score distributions
    less_than_400 = (filtered_data["score"] < 400).sum()
    exactly_400 = (filtered_data["score"] == 400).sum()
    more_than_400 = (filtered_data["score"] > 400).sum()

    df_bins = pd.DataFrame(
        {
            "Score_range": [
                "0 to 400 (excluded)",
                "Exactly 400",
                f"400 to max score ({filtered_data['score'].max():.2f})",
            ],
            "Count": [less_than_400, exactly_400, more_than_400],
        }
    )

    col1, col2 = st.columns(2)
    with col1:
        st.write("**Score Statistics**")
        st.dataframe(filtered_data["score"].describe(), use_container_width=True)
    with col2:
        st.write("**Score Distribution**")
        st.dataframe(df_bins, hide_index=True, use_container_width=True)

    # RBP Analysis
    st.header("RBP Analysis")
    plot_rbp_counts_st(filtered_data)

    # Network Generation
    st.header("Network Generation")
    with st.form("network_parameters"):
        physics_enabled = st.checkbox("Enable Physics Simulation", value=True)
        min_score = st.number_input("Minimum Score Threshold", value=400.0)
        generate_button = st.form_submit_button(label="Generate Network")

    if generate_button:
        network_data = filtered_data[filtered_data["score"] >= min_score]
        st.write(
            f"Network filtered to score >= {min_score} contains {len(network_data)} interactions "
            f"and {network_data['source'].nunique()} unique RBPs."
        )

        net = create_network(network_data, physics=physics_enabled)
        render_pyvis_network(net)


if __name__ == "__main__":
    main()
