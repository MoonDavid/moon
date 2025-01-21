import streamlit as st
import requests
import pandas as pd
import urllib.parse
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
import tempfile
import base64
import os
import numpy as np
import streamlit.components.v1 as components

# =========================================================
#                   GLOBAL CONFIG
# =========================================================
DEFAULT_FORMAT = 'json'
DEFAULT_SPECIES = 9606  # Homo sapiens
DEFAULT_LIMIT = 10
BASE_URL = "https://string-db.org/api"

# Configure Streamlit page
st.set_page_config(page_title="STRING API Wrapper", layout="wide")


# =========================================================
#                   API HELPER FUNCTIONS
# =========================================================
def df_from_uniprots(df,uniprots):
    """Create a DataFrame from a list of UniProt IDs."""
    return df[df['UniProtKB-AC'].isin(uniprots)].reset_index(drop=True)
def string_to_uniprots(string_ids):
    """Convert a list of STRING IDs to a list of UniProt IDs."""
    data=st.session_state['df']
    #create dict from columns UniProtKB-AC and STRING if STRING not NaN
    d=dict(zip(data['STRING'],data['UniProtKB-AC']))
    return d.get(string_ids)
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
    :return: Parsed response content as a pandas DataFrame
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
    proteins_str = "%0d".join(proteins)  # STRING API uses newline as separator, encoded as %0d
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


# =========================================================
#        NETWORK ANALYSIS & VISUALIZATION FUNCTIONS
# =========================================================
def calculate_centralities(G: nx.Graph) -> pd.DataFrame:
    """
    Calculate various centralities of a NetworkX graph.

    :param G: NetworkX graph
    :return: Pandas DataFrame containing centralities
    """
    centrality_df = pd.DataFrame({
        "Degree": nx.degree_centrality(G),
        "Betweenness": nx.betweenness_centrality(G),
        "Closeness": nx.closeness_centrality(G),
        "Eigenvector": nx.eigenvector_centrality(G, max_iter=1000)
    }).round(4)
    return centrality_df


def create_pyvis_network2(
    G: nx.Graph,
    physics_enabled: bool,
    node_size: int,
    edge_width: int,
    color_scheme: str,
    weight_col
) -> Network:
    """
    Build a PyVis Network visualization from a NetworkX graph.

    :param G: NetworkX Graph to visualize
    :param physics_enabled: Whether to enable physics simulation
    :param node_size: Size of the nodes
    :param edge_width: Width of the edges
    :param color_scheme: Visualization color scheme
    :return: PyVis Network object
    """


    # Configure physics
    nt = Network(
        height="600px",
        width="100%",
        bgcolor="#ffffff",
        font_color="black"
    )
    nt.repulsion()

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
        "font": {"size": 12}
      },
      "edges": {
        "smooth": false
      },
      "physics": {
        "barnesHut": {
          "gravitationalConstant": -5000,
          "springLength": 95
        },
        "minVelocity": 0.75,
        "stabilization": {
          "iterations": 200
        }
      }
    }
    """)

    return nt
def create_pyvis_network_claude(
    G: nx.Graph,
    physics_enabled: bool,
    node_size: int,
    edge_width: int,
    color_scheme: str,
    weight_col: str
) -> Network:
    """
    Build a PyVis Network visualization from a NetworkX graph with improved stabilization.

    :param G: NetworkX Graph to visualize
    :param physics_enabled: Whether to enable physics simulation
    :param node_size: Size of the nodes
    :param edge_width: Width of the edges
    :param color_scheme: Visualization color scheme
    :param weight_col: Which column in edge_attr to treat as weight
    :return: PyVis Network object
    """
    nt = Network(
        height="600px",
        width="100%",
        bgcolor="#ffffff",
        font_color="black",
        notebook=False
    )

    # Disable physics by default - we'll only use it for initial layout
    nt.toggle_physics(False)

    # More conservative repulsion settings to prevent 'exploding'
    nt.repulsion(
        node_distance=100,      # Reduced from 120
        central_gravity=0.5,    # Increased from 0.2 for more centralization
        spring_length=100,      # Reduced from 150
        spring_strength=0.15,   # Increased from 0.05 for more stability
        damping=0.95           # Increased from 0.09 for more dampening
    )

    # Color scheme implementation remains the same
    if color_scheme == "Degree-based":
        degree_dict = dict(G.degree())
        max_degree = max(degree_dict.values()) if degree_dict else 1
        for node in G.nodes():
            deg_val = int((degree_dict[node] / max_degree) * 255)
            node_color = f"#{deg_val:02x}0000"
            nt.add_node(node, label=node, size=node_size, color=node_color)
    elif color_scheme == "Community-based":
        from networkx.algorithms import community
        communities = community.greedy_modularity_communities(G)
        color_map = plt.cm.get_cmap("Set3")(np.linspace(0, 1, len(communities)))
        node_colors = {}
        for i, comm in enumerate(communities):
            rgb = color_map[i][:3]
            color_hex = "#{:02x}{:02x}{:02x}".format(
                int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255)
            )
            for node in comm:
                node_colors[node] = color_hex
        for node in G.nodes():
            color = node_colors.get(node, "#97C2FC")
            nt.add_node(node, label=node, size=node_size, color=color)
    else:
        for node in G.nodes():
            nt.add_node(node, label=node, size=node_size)

    # Add edges
    for edge in G.edges(data=True):
        weight = edge[2].get(weight_col, 1)
        nt.add_edge(edge[0], edge[1], value=weight, width=edge_width)

    # Updated physics and stabilization options
    nt.set_options("""
    var options = {
      "nodes": {
        "font": {"size": 12},
        "fixed": {
          "x": false,
          "y": false
        }
      },
      "edges": {
        "smooth": {
          "type": "continuous",
          "forceDirection": "none"
        }
      },
      "physics": {
        "enabled": false,
        "solver": "forceAtlas2Based",
        "forceAtlas2Based": {
          "gravitationalConstant": -50,
          "centralGravity": 0.01,
          "springLength": 100,
          "springConstant": 0.08,
          "damping": 0.95,
          "avoidOverlap": 1
        },
        "stabilization": {
          "enabled": true,
          "iterations": 1000,
          "updateInterval": 25,
          "onlyDynamicEdges": false,
          "fit": true
        },
        "minVelocity": 0.75,
        "maxVelocity": 30
      },
      "interaction": {
        "navigationButtons": true,
        "keyboard": true,
        "dragNodes": true,
        "dragView": true,
        "zoomView": true
      }
    }
    """)

    # Only enable physics if explicitly requested
    if physics_enabled:
        nt.toggle_physics(True)

    return nt
def create_pyvis_network(
    G: nx.Graph,
    physics_enabled: bool,
    node_size: int,
    edge_width: int,
    color_scheme: str,
    weight_col: str
) -> Network:
    """
    Build a PyVis Network visualization from a NetworkX graph.

    :param G: NetworkX Graph to visualize
    :param physics_enabled: Whether to enable physics simulation
    :param node_size: Size of the nodes
    :param edge_width: Width of the edges
    :param color_scheme: Visualization color scheme
    :param weight_col: Which column in edge_attr to treat as weight
    :return: PyVis Network object
    """
    nt = Network(
        height="600px",
        width="100%",
        bgcolor="#ffffff",
        font_color="black",
        notebook=False
    )

    # If your network is large, you might want to tweak the repulsion
    # settings (node_distance, spring_length, etc.) so it doesn’t
    # ‘explode’ or keep moving:
    nt.repulsion(
        node_distance=120,      # more distance = more spread out
        central_gravity=0.2,    # 0.0 to 1.0
        spring_length=150,      # length of the edges
        spring_strength=0.05,   # higher = more rigid edges
        damping=0.09
    )

    # Toggle physics only if user checks the box,
    # otherwise forcibly disable it:
    if physics_enabled:
        nt.toggle_physics(True)
    else:
        nt.toggle_physics(False)

    # (Optional) Show a physics GUI in the final HTML,
    # so you can manually tweak the layout in the browser:
    # nt.show_buttons(filter_=['physics'])

    # Color scheme: example of degree-based or community-based
    if color_scheme == "Degree-based":
        degree_dict = dict(G.degree())
        max_degree = max(degree_dict.values()) if degree_dict else 1
        for node in G.nodes():
            # Simple gradient from 0 to 255 in red channel
            deg_val = int((degree_dict[node] / max_degree) * 255)
            node_color = f"#{deg_val:02x}0000"
            nt.add_node(node, label=node, size=node_size, color=node_color)
    elif color_scheme == "Community-based":
        from networkx.algorithms import community
        communities = community.greedy_modularity_communities(G)
        color_map = plt.cm.get_cmap("Set3")(np.linspace(0, 1, len(communities)))
        node_colors = {}
        for i, comm in enumerate(communities):
            rgb = color_map[i][:3]
            color_hex = "#{:02x}{:02x}{:02x}".format(
                int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255)
            )
            for node in comm:
                node_colors[node] = color_hex
        for node in G.nodes():
            color = node_colors.get(node, "#97C2FC")
            nt.add_node(node, label=node, size=node_size, color=color)
    else:
        # Default color scheme
        for node in G.nodes():
            nt.add_node(node, label=node, size=node_size)

    # Add edges
    for edge in G.edges(data=True):
        weight = edge[2].get(weight_col, 1)
        nt.add_edge(edge[0], edge[1], width=edge_width)

    # Provide advanced config for physics/stabilization
    # Increase stabilization iterations for larger graphs:
    nt.set_options(""" 
    var options = {
      "nodes": {
        "font": {"size": 12}
      },
      "edges": {
        "smooth": false
      },
      "physics": {
        "enabled": true,
        "solver": "barnesHut",
        "barnesHut": {
          "gravitationalConstant": -8000,
          "centralGravity": 0.2,
          "springLength": 150,
          "springConstant": 0.04,
          "damping": 0.09,
          "avoidOverlap": 1
        },
        "stabilization": {
          "enabled": true,
          "iterations": 800
        }
      }
    }
    """)

    return nt
def render_pyvis_network(nt: Network):
    """
    Render the PyVis Network in a Streamlit app and provide a download button.

    :param nt: PyVis Network object
    """
    try:
        path = "network.html"
        nt.save_graph(path)
        with open(path, 'r', encoding='utf-8') as file:
            html_content = file.read()
        components.html(html_content, height=600)

        # Download button
        with open(path, 'rb') as file:
            st.download_button(
                label="Download Network Visualization",
                data=file,
                file_name="network.html",
                mime="text/html"
            )
    except Exception as viz_error:
        st.error(f"Visualization Error: {viz_error}")


# =========================================================
#                 UI LOGIC FUNCTIONS
# =========================================================
def show_interactions_ui(species: int, limit: int):
    """
    UI and logic for retrieving and visualizing protein-protein interactions.
    """
    st.header("🔗 Get Protein-Protein Interactions")

    # Attempt to fetch available gene lists from session state
    available_gene_lists = {
        key: value for key, value in st.session_state.items() if isinstance(value, list)
    }

    if not available_gene_lists:
        st.warning("No pre-loaded gene lists found in session state. Please add some.")
        return

    # Let user pick which gene list to use
    sorted_gene_list_keys = list(available_gene_lists.keys())
    selected_gene_list = st.selectbox(
        "Select the gene list to convert:",
        options=sorted_gene_list_keys,
        format_func=lambda x: x if x else "No list available"
    )

    proteins_input = ','.join(available_gene_lists[selected_gene_list])

    # Sidebar parameters
    required_score = st.sidebar.slider(
        "Required Score",
        min_value=0,
        max_value=1000,
        value=400,
        step=50,
        help="Minimum required score for interactions."
    )

    if st.button("Get Interactions"):
        proteins = [
            p.strip()
            for p in proteins_input.replace(',', '\n').split('\n')
            if p.strip()
        ]
        if proteins:
            try:
                with st.spinner("Fetching interactions..."):
                    interactions_df = get_interactions(
                        proteins,
                        required_score=required_score,
                        species=species,
                        limit=limit
                    )

                    interactions_df['UniProtSource'] = interactions_df.iloc[:, 0].map(string_to_uniprots)
                    interactions_df['UniProtTarget'] = interactions_df.iloc[:, 1].map(string_to_uniprots)
                    st.dataframe(interactions_df)
                    # Store data in session
                    st.session_state["interactions_df"] = interactions_df
                    st.session_state["G"] = nx.from_pandas_edgelist(
                        interactions_df,
                        source=interactions_df.columns[2],
                        target=interactions_df.columns[3],
                        edge_attr=interactions_df.columns[2],
                        create_using=nx.Graph()
                    )
            except Exception as e:
                st.error(f"Error: {e}")
                return

    # If data is available, show it
    if "interactions_df" in st.session_state and "G" in st.session_state:
        interactions_df = st.session_state.interactions_df
        G = st.session_state.G

        source_col, target_col, weight_col = interactions_df.columns[:3]

        st.success("Interactions retrieved!")
        st.dataframe(interactions_df)
        st.write(f"**Number of Interactions:** {len(interactions_df)}")

        # Network Analysis
        density = nx.density(G)
        st.write(f"**Network Density:** {density:.4f}")

        centrality_df = calculate_centralities(G)
        st.dataframe(centrality_df)

        top_degree = centrality_df["Degree"].sort_values(ascending=False).head(5).index.tolist()
        st.write(f"**Top 5 Nodes by Degree Centrality:** {', '.join(top_degree)}")


        st.subheader("Networkx Visualization")
        plt.figure(figsize=(10, 6))

        # Draw the graph
        nx.draw(G,
                with_labels=True,
                node_color='lightblue',
                edge_color='gray',
                node_size=500,
                font_size=16,
                font_weight='bold')

        # Display the graph in Streamlit
        st.pyplot(plt)
        # Create and render the network
        st.subheader("Pyvis Interactive Network Visualization")

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
        nt = create_pyvis_network(
            G,
            physics_enabled=physics_enabled,
            node_size=node_size,
            edge_width=edge_width,
            color_scheme=color_scheme,
            weight_col=weight_col
        )

        render_pyvis_network(nt)


def show_protein_info_ui(species: int):
    """
    UI and logic for retrieving detailed protein information.
    """
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
                    info_df = get_protein_info(protein_id, species=species)
                st.success("Protein information retrieved!")
                st.dataframe(info_df)
            except Exception as e:
                st.error(f"Error: {e}")
        else:
            st.warning("Please enter a protein ID.")


def show_enrichment_ui(species: int):
    """
    UI and logic for performing enrichment analysis on a set of proteins.
    """
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
                    enrichment_df = get_enrichment(proteins, species=species)
                st.success("Enrichment analysis completed!")
                st.dataframe(enrichment_df)
            except Exception as e:
                st.error(f"Error: {e}")
        else:
            st.warning("Please enter at least one protein ID.")


# =========================================================
#                        MAIN APP
# =========================================================
def main():
    st.title("🔗 STRING Database API Explorer")
    st.markdown("""
    This application allows you to interact with the [STRING database](https://string-db.org/) 
    using the STRING API. You can map external protein identifiers to STRING IDs, retrieve 
    protein-protein interactions, get detailed protein information, and perform enrichment analysis.
    """)

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

    species = int(species)

    # Route to appropriate UI function
    if api_method == "Get Protein-Protein Interactions":
        show_interactions_ui(species, limit)
    elif api_method == "Get Protein Information":
        show_protein_info_ui(species)
    else:  # "Enrichment Analysis"
        show_enrichment_ui(species)


if __name__ == "__main__":
    main()
