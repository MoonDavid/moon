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


class StringAPIWrapper:
    """
    A Python wrapper for the STRING API.
    """

    BASE_URL = "https://string-db.org/api"

    def __init__(self, format='json', species=9606, limit=10, **kwargs):
        """
        Initializes the StringAPIWrapper.

        :param format: Response format (json, xml, etc.)
        :param species: NCBI taxonomy identifier (default: 9606 for Homo sapiens)
        :param limit: Maximum number of results to return
        :param kwargs: Additional parameters
        """
        self.format = format
        self.species = species
        self.limit = limit
        self.params = kwargs

    def _construct_url(self, endpoint):
        """
        Constructs the full API URL.

        :param endpoint: Specific API endpoint
        :return: Full URL as a string
        """
        return f"{self.BASE_URL}/{self.format}/{endpoint}"

    def _handle_response(self, response):
        """
        Handles the API response, checking for errors.

        :param response: Response object from requests
        :return: Parsed response content as pandas DataFrame
        :raises Exception: If the response contains an error
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

    def get_interactions(self, proteins, required_score=400):
        """
        Retrieves protein-protein interactions.

        :param proteins: List of protein identifiers (e.g., UniProt IDs)
        :param required_score: Minimum required score (default: 400)
        :return: DataFrame of interactions
        """
        endpoint = "network"
        proteins_str = "%0d".join(proteins)  # STRING API uses newline as separator encoded as %0d
        params = {
            "identifiers": proteins_str,
            "species": self.species,
            "required_score": required_score,
            "limit": self.limit
        }
        url = self._construct_url(endpoint)
        response = requests.post(url, data=params)
        return self._handle_response(response)

    def get_string_ids(self, identifiers, method='get_string_ids'):
        """
        Maps external identifiers to STRING identifiers.

        :param identifiers: List of external identifiers (e.g., gene symbols)
        :param method: Method to map identifiers (default: get_string_ids)
        :return: DataFrame mapping identifiers
        """
        endpoint = method
        identifiers_str = "%0d".join(identifiers)
        params = {
            "identifiers": identifiers_str,
            "species": self.species,
            "limit": self.limit
        }
        url = self._construct_url(endpoint)
        response = requests.post(url, data=params)
        return self._handle_response(response)

    def get_protein_info(self, protein_id):
        """
        Retrieves detailed information about a specific protein.

        :param protein_id: STRING identifier of the protein
        :return: DataFrame with protein information
        """
        endpoint = "info"
        params = {
            "identifiers": protein_id,
            "species": self.species
        }
        url = self._construct_url(endpoint)
        response = requests.post(url, data=params)
        return self._handle_response(response)

    def get_enrichment(self, proteins):
        """
        Performs enrichment analysis on a set of proteins.

        :param proteins: List of STRING protein identifiers
        :return: DataFrame with enrichment analysis results
        """
        endpoint = "enrichment"
        proteins_str = "%0d".join(proteins)
        params = {
            "identifiers": proteins_str,
            "species": self.species
        }
        url = self._construct_url(endpoint)
        response = requests.post(url, data=params)
        return self._handle_response(response)

    # Add more methods as needed based on STRING API endpoints

def main():
    st.set_page_config(page_title="STRING API Wrapper", layout="wide")
    st.title("🔗 STRING Database API Explorer")

    st.markdown("""
    This application allows you to interact with the [STRING database](https://string-db.org/) using the STRING API.
    You can map external protein identifiers to STRING IDs, retrieve protein-protein interactions, get detailed protein information, and perform enrichment analysis.
    """)

    # Initialize the API wrapper
    api = StringAPIWrapper()

    # Sidebar for user inputs
    st.sidebar.header("Input Parameters")

    # Select API method
    api_method = st.sidebar.selectbox(
        "Choose API Method",
        ("Map Identifiers to STRING IDs", "Get Protein-Protein Interactions", "Get Protein Information", "Enrichment Analysis")
    )

    # Common parameters
    species = st.sidebar.text_input("Species NCBI ID", value="9606", help="Default is 9606 for Homo sapiens.")
    limit = st.sidebar.number_input("Result Limit", min_value=1, max_value=1000, value=10, step=1)

    # Update API wrapper with common parameters
    api.species = int(species)
    api.limit = limit

    if api_method == "Map Identifiers to STRING IDs":
        st.header("🔍 Map External Identifiers to STRING IDs")
        identifiers_input = st.text_area(
            "Enter Protein Identifiers",
            "TP53, BRCA1, EGFR",
            help="Enter protein identifiers separated by commas or newlines."
        )
        if st.button("Map Identifiers"):
            identifiers = [id.strip() for id in identifiers_input.replace(',', '\n').split('\n') if id.strip()]
            if identifiers:
                try:
                    with st.spinner("Mapping identifiers..."):
                        mapping_df = api.get_string_ids(identifiers)
                    st.success("Mapping completed!")
                    st.dataframe(mapping_df)
                except Exception as e:
                    st.error(f"Error: {e}")
            else:
                st.warning("Please enter at least one identifier.")

    elif api_method == "Get Protein-Protein Interactions":
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
        required_score = st.sidebar.slider("Required Score", min_value=0, max_value=1000, value=400, step=50, help="Minimum required score for interactions.")
        if st.button("Get Interactions"):
            proteins = [p.strip() for p in proteins_input.replace(',', '\n').split('\n') if p.strip()]
            if proteins:
                try:
                    with st.spinner("Fetching interactions..."):
                        interactions_df = api.get_interactions(proteins, required_score=required_score)
                    st.success("Interactions retrieved!")
                    st.dataframe(interactions_df)
                    source_col, target_col, weight_col = interactions_df.columns[:3]
                    st.write(f"**Number of Interactions:** {len(interactions_df)}")
                    G = nx.from_pandas_edgelist(
                        interactions_df,
                        source=source_col,
                        target=target_col,
                        edge_attr=weight_col,
                        create_using=nx.Graph()
                    )

                    st.subheader("Protein Interaction Network")


                    # Density
                    density = nx.density(G)
                    st.write(f"**Network Density:** {density:.4f}")

                    # Centrality Measures
                    st.write("**Centrality Measures:**")
                    centrality_df = pd.DataFrame({
                        "Degree": nx.degree_centrality(G),
                        "Betweenness": nx.betweenness_centrality(G),
                        "Closeness": nx.closeness_centrality(G),
                        "Eigenvector": nx.eigenvector_centrality(G, max_iter=1000)
                    }).round(4)

                    st.dataframe(centrality_df)

                    # Optionally, highlight top central nodes
                    top_degree = centrality_df['Degree'].sort_values(ascending=False).head(5).index.tolist()
                    st.write(f"**Top 5 Nodes by Degree Centrality:** {', '.join(top_degree)}")
                except Exception as e:
                    st.error(f"Error: {e}")

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
                        info_df = api.get_protein_info(protein_id)
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
                        enrichment_df = api.get_enrichment(proteins)
                    st.success("Enrichment analysis completed!")
                    st.dataframe(enrichment_df)
                except Exception as e:
                    st.error(f"Error: {e}")
            else:
                st.warning("Please enter at least one protein ID.")



if __name__ == "__main__":
    main()
