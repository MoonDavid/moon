import streamlit as st
import requests
import pandas as pd
import urllib.parse
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
import networkx as nx
import tempfile
import base64
import os
import numpy as np
# Configuration parameters
DEFAULT_FORMAT = 'json'
DEFAULT_SPECIES = 9606  # Homo sapiens
DEFAULT_LIMIT = 10
BASE_URL = "https://string-db.org/api"

def construct_url(endpoint, format=DEFAULT_FORMAT):
    """
    Constructs the full API URL.

    :param endpoint: Specific API endpoint
    :param format: Response format (json, xml, etc.)
    :return: Full URL as a string
    """
    return f"{BASE_URL}/{format}/{endpoint}"


def handle_response(response):
    """
    Handles the API response, checking for errors and parsing JSON.

    :param response: Response object from requests
    :return: Parsed response content as pandas DataFrame
    :raises Exception: If the response contains an error or unexpected structure
    """
    if response.status_code != 200:
        raise Exception(f"API request failed with status code {response.status_code}: {response.text}")
    try:
        data = response.json()
        if isinstance(data, list):
            return pd.DataFrame(data)
        elif isinstance(data, dict):
            return pd.DataFrame([data])
        else:
            raise Exception("Unexpected JSON structure")
    except ValueError:
        raise Exception("Failed to parse JSON response")

@st.cache_data
def get_interactions(proteins, required_score=400, format=DEFAULT_FORMAT, species=DEFAULT_SPECIES, limit=DEFAULT_LIMIT):
    """
    Retrieves protein-protein interactions.

    :param proteins: List of protein identifiers (e.g., UniProt IDs)
    :param required_score: Minimum required score (default: 400)
    :param format: Response format
    :param species: NCBI taxonomy identifier
    :param limit: Maximum number of results to return
    :return: DataFrame of interactions
    """
    endpoint = "network"
    proteins_str = "%0d".join(proteins)  # STRING API uses newline as separator encoded as %0d
    params = {
        "identifiers": proteins_str,
        "species": species,
        "required_score": required_score,
        "limit": limit
    }
    url = construct_url(endpoint, format)
    response = requests.post(url, data=params)
    return handle_response(response)

@st.cache_data
def get_string_ids(identifiers, method='get_string_ids', format=DEFAULT_FORMAT, species=DEFAULT_SPECIES, limit=DEFAULT_LIMIT):
    """
    Maps external identifiers to STRING identifiers.

    :param identifiers: List of external identifiers (e.g., gene symbols)
    :param method: Method to map identifiers (default: get_string_ids)
    :param format: Response format
    :param species: NCBI taxonomy identifier
    :param limit: Maximum number of results to return
    :return: DataFrame mapping identifiers
    """
    endpoint = method
    identifiers_str = "%0d".join(identifiers)
    params = {
        "identifiers": identifiers_str,
        "species": species,
        "limit": limit
    }
    url = construct_url(endpoint, format)
    response = requests.post(url, data=params)
    return handle_response(response)

@st.cache_data
def get_protein_info(protein_id, format=DEFAULT_FORMAT, species=DEFAULT_SPECIES):
    """
    Retrieves detailed information about a specific protein.

    :param protein_id: STRING identifier of the protein
    :param format: Response format
    :param species: NCBI taxonomy identifier
    :return: DataFrame with protein information
    """
    endpoint = "info"
    params = {
        "identifiers": protein_id,
        "species": species
    }
    url = construct_url(endpoint, format)
    response = requests.post(url, data=params)
    return handle_response(response)

@st.cache_data
def get_enrichment(proteins, format=DEFAULT_FORMAT, species=DEFAULT_SPECIES):
    """
    Performs enrichment analysis on a set of proteins.

    :param proteins: List of STRING protein identifiers
    :param format: Response format
    :param species: NCBI taxonomy identifier
    :return: DataFrame with enrichment analysis results
    """
    endpoint = "enrichment"
    proteins_str = "%0d".join(proteins)
    params = {
        "identifiers": proteins_str,
        "species": species
    }
    url = construct_url(endpoint, format)
    response = requests.post(url, data=params)
    return handle_response(response)

def main():
    st.set_page_config(page_title="STRING API Wrapper", layout="wide")
    st.title("🔗 STRING Database API Explorer")

    st.markdown("""
    This application allows you to interact with the [STRING database](https://string-db.org/) using the STRING API.
    You can map external protein identifiers to STRING IDs, retrieve protein-protein interactions, get detailed protein information, and perform enrichment analysis.
    """)

    # Initialize the API wrapper
    #api = StringAPIWrapper()

    # Sidebar for user inputs
    st.sidebar.header("Input Parameters")

    # Select API method
    api_method = st.sidebar.selectbox(
        "Choose API Method",
        ("Get Protein-Protein Interactions", "Get Protein Information", "Enrichment Analysis")
    )

    # Common parameters
    species = st.sidebar.text_input("Species NCBI ID", value="9606", help="Default is 9606 for Homo sapiens.")
    limit = st.sidebar.number_input("Result Limit", min_value=1, max_value=1000, value=10, step=1)

    # Update API wrapper with common parameters
    species = int(species)
    limit = limit

    if api_method == "Get Protein-Protein Interactions":
        st.header("🔗 Get Protein-Protein Interactions")
        available_gene_lists = {key: value for key, value in st.session_state.items() if isinstance(value, list)}
        if available_gene_lists:
            sorted_gene_list_keys = list(available_gene_lists.keys())
            selected_gene_list = st.selectbox(
                "Select the gene list to convert:",
                options=sorted_gene_list_keys,
                format_func=lambda x: x if x else "No list available"
            )
        proteins_input = ','.join(available_gene_lists[selected_gene_list])
        required_score = st.sidebar.slider("Required Score", min_value=0, max_value=1000, value=400, step=50,
                                           help="Minimum required score for interactions.")

        if st.button("Get Interactions"):
            proteins = [p.strip() for p in proteins_input.replace(',', '\n').split('\n') if p.strip()]
            if proteins:
                try:
                    with st.spinner("Fetching interactions..."):
                        interactions_df = get_interactions(proteins, required_score=required_score)
                        # Store the data in session state
                        st.session_state.interactions_df = interactions_df
                        st.session_state.G = nx.from_pandas_edgelist(
                            interactions_df,
                            source=interactions_df.columns[0],
                            target=interactions_df.columns[1],
                            edge_attr=interactions_df.columns[2],
                            create_using=nx.Graph()
                        )
                except Exception as e:
                    st.error(f"Error: {e}")
                    return

        # Only show visualization options and results if we have data
        if 'interactions_df' in st.session_state and 'G' in st.session_state:
            interactions_df = st.session_state.interactions_df
            G = st.session_state.G
            source_col, target_col, weight_col = interactions_df.columns[:3]

            st.success("Interactions retrieved!")
            st.dataframe(interactions_df)
            st.write(f"**Number of Interactions:** {len(interactions_df)}")

            # Network Analysis
            density = nx.density(G)
            st.write(f"**Network Density:** {density:.4f}")

            centrality_df = pd.DataFrame({
                "Degree": nx.degree_centrality(G),
                "Betweenness": nx.betweenness_centrality(G),
                "Closeness": nx.closeness_centrality(G),
                "Eigenvector": nx.eigenvector_centrality(G, max_iter=1000)
            }).round(4)

            st.dataframe(centrality_df)
            top_degree = centrality_df['Degree'].sort_values(ascending=False).head(5).index.tolist()
            st.write(f"**Top 5 Nodes by Degree Centrality:** {', '.join(top_degree)}")

            # Visualization Options (moved to main area for better UX)
            st.subheader("Network Visualization Options")
            physics_enabled = st.checkbox("Enable Physics", value=False)
            node_size = st.slider("Node Size", min_value=10, max_value=50, value=20)
            col1, col2 = st.columns(2)
            with col1:
                edge_width = st.slider("Edge Width", min_value=1, max_value=10, value=2)
            with col2:
                color_scheme = st.selectbox(
                    "Color Scheme",
                    options=["Default", "Degree-based", "Community-based"]
                )

            # Network Visualization
            st.subheader("Interactive Network Visualization")

            # Create Pyvis network
            nt = Network(
                height="600px",
                width="100%",
                bgcolor="#ffffff",
                font_color="black"
            )

            # Configure physics
            nt.toggle_physics(physics_enabled)

            # Add nodes and edges with custom styling
            if color_scheme == "Degree-based":
                degree_dict = dict(G.degree())
                max_degree = max(degree_dict.values())
                for node in G.nodes():
                    node_color = f"#{int((degree_dict[node] / max_degree) * 255):02x}0000"
                    nt.add_node(node, label=node, size=node_size, color=node_color)
            elif color_scheme == "Community-based":
                communities = nx.community.greedy_modularity_communities(G)
                color_map = plt.cm.get_cmap('Set3')(np.linspace(0, 1, len(communities)))
                node_colors = {}
                for i, comm in enumerate(communities):
                    for node in comm:
                        rgb = color_map[i][:3]
                        node_colors[node] = f"#{int(rgb[0] * 255):02x}{int(rgb[1] * 255):02x}{int(rgb[2] * 255):02x}"
                for node in G.nodes():
                    nt.add_node(node, label=node, size=node_size, color=node_colors.get(node, "#97C2FC"))
            else:
                for node in G.nodes():
                    nt.add_node(node, label=node, size=node_size)

            # Add edges with weights
            for edge in G.edges(data=True):
                weight = edge[2].get(weight_col, 1)
                nt.add_edge(edge[0], edge[1], value=weight, width=edge_width)

            # Set other visualization options
            nt.set_options("""
                var options = {
                    "nodes": {
                        "font": {
                            "size": 12
                        }
                    },
                    "edges": {
                        "color": {
                            "inherit": true
                        },
                        "smooth": false
                    },
                    "physics": {
                        "barnesHut": {
                            "gravitationalConstant": -5000,
                            "springLength": 95
                        },
                        "minVelocity": 0.75
                    }
                }
            """)

            # Save and display the network
            try:
                path = "network.html"
                nt.save_graph(path)
                with open(path, 'r', encoding='utf-8') as file:
                    html_content = file.read()
                st.components.v1.html(html_content, height=600)
            except Exception as viz_error:
                st.error(f"Visualization Error: {viz_error}")

    elif api_method == "Get Protein Information":
        st.header("ℹ️ Get Detailed Protein Information")
        protein_id = st.text_input(
            "Enter STRING Protein ID",
            "9606.ENSP00000269305",
            help="Enter a single STRING protein ID."
        )
        if st.button("Get Information"):
            if protein_id:
                try:
                    with st.spinner("Fetching protein information..."):
                        info_df = get_protein_info(protein_id)
                    st.success("Protein information retrieved!")
                    st.dataframe(info_df)
                except Exception as e:
                    st.error(f"Error: {e}")
            else:
                st.warning("Please enter a protein ID.")

    elif api_method == "Enrichment Analysis":
        st.header("📈 Perform Enrichment Analysis")
        proteins_input = st.text_area(
            "Enter STRING Protein IDs",
            "9606.ENSP00000269305,9606.ENSP00000269305",
            help="Enter STRING protein IDs separated by commas or newlines."
        )
        if st.button("Run Enrichment Analysis"):
            proteins = [p.strip() for p in proteins_input.replace(',', '\n').split('\n') if p.strip()]
            if proteins:
                try:
                    with st.spinner("Performing enrichment analysis..."):
                        enrichment_df = get_enrichment(proteins)
                    st.success("Enrichment analysis completed!")
                    st.dataframe(enrichment_df)
                except Exception as e:
                    st.error(f"Error: {e}")
            else:
                st.warning("Please enter at least one protein ID.")



if __name__ == "__main__":
    main()
