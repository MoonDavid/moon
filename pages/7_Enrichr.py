import streamlit as st
import gseapy as gp
import pandas as pd
from io import StringIO


def df_from_uniprots(df,uniprots):
    """Create a DataFrame from a list of UniProt IDs."""
    return df[df['UniProtKB-AC'].isin(uniprots)].reset_index(drop=True)
# Function to perform enrichment analysis
def perform_enrichment(gene_list, background_list, libraries):
    """Performs enrichment analysis using gseapy."""

    try:
        enr = gp.enrichr(gene_list=gene_list,
                         gene_sets=libraries,
                         background=background_list,
                         no_plot=True,
                         )
        return enr.results
    except Exception as e:
        st.error(f"Error during enrichment analysis: {e}")
        return None


# Function to read gene lists from text area, returns set
def read_gene_list(text):
    if not text:
        return set()
    return {line.strip() for line in text.splitlines()}


# Function to read background from text area, returns set
def read_background(text):
    if not text:
        return None
    return {line.strip() for line in text.splitlines()}


# Main Streamlit app
def main():
    st.title("Enrichr Gene Set Enrichment Analysis")

    # Initialize session state for gene lists
    if 'gene_lists' not in st.session_state:
        st.session_state['gene_lists'] = {}

    if 'background_lists' not in st.session_state:
        st.session_state['background_lists'] = {}

    # Gene list management
    st.header("Gene List Input")
    with st.form(key="gene_list_form"):


        # Display gene lists in selectbox
        st.header("Select Gene List")
        available_gene_lists = {key: value for key, value in st.session_state.items() if isinstance(value, list)}
        if available_gene_lists:
            # Sort the gene lists to show the oldest first based on insertion order
            sorted_gene_list_keys = list(available_gene_lists.keys())
            selected_gene_list = st.selectbox(
                "Select the gene list to use:",
                options=sorted_gene_list_keys,
                format_func=lambda x: x if x else "No list available"
            )

        else:
            st.warning("No gene list available in session state.")
            selected_gene_list = None

        # Background option (optional)
        st.header("Background Genes (Optional)")
        background_option = st.radio("Background Gene Input", ("None", "Custom List"), key="bg_radio")

        background_list = None
        if background_option == "Custom List":
            # Display background lists in selectbox
            available_bg_lists = {key: value for key, value in st.session_state.items() if isinstance(value, list)}
            if available_bg_lists:
                sorted_bg_list_keys = list(available_bg_lists.keys())
                background_list = st.selectbox(
                    "Select the background list to use:",
                    options=sorted_bg_list_keys,
                    format_func=lambda x: x if x else "No list available",
                    key='select_bg_list'
                )
            else:
                st.warning("No background list available in session state.")
                background_list = None
            submit_button = st.form_submit_button(label="Run Analysis Enrichr")
        # Enrichment libraries
        st.header("Enrichr Libraries")
        available_libraries = gp.get_library_name()
        selected_libraries = st.multiselect("Select Enrichment Libraries", options=available_libraries,
                                            default=["GO_Biological_Process_2023"], key='select_lib')
        submit_button = st.form_submit_button(label="Run Analysis")
    # Run analysis button
    if submit_button:
        uniprots = available_gene_lists.get(selected_gene_list, [])
        df = df_from_uniprots(st.session_state['df'], uniprots)
        selected_gene_list = df.loc[df['Entrez'].notna(), 'Entrez'].values.tolist()
        if not selected_gene_list:
            st.error("Please add and/or select at least one gene list")
            return
        if not selected_libraries:
            st.error("Please select at least one library")
            return

        with st.spinner("Running enrichment analysis..."):
            results_df = perform_enrichment(list(selected_gene_list), background_list, selected_libraries)

        if results_df is not None and not results_df.empty:
            st.success("Enrichment analysis completed!")
            st.subheader("Enrichment Results")
            st.dataframe(results_df)

            # download button for results
            csv = results_df.to_csv(index=False).encode('utf-8')
            st.download_button(label="Download CSV",
                               data=csv,
                               file_name='enrichr_results.csv',
                               mime='text/csv')
        else:
            st.error("No results obtained or an error occurred during analysis.")

        if results_df is not None and not results_df.empty:
            st.session_state['enrichment_results_enrichr'] = results_df  # Save results in session state

            st.success("Enrichment analysis completed!")
        else:
            st.error("No results obtained or an error occurred during analysis.")
    # Display existing enrichment results if available
    if 'enrichment_results_enrichr' in st.session_state:
        st.subheader("Enrichment Results")
        st.dataframe(st.session_state['enrichment_results_enrichr'])

        # Download button for results
        csv = st.session_state['enrichment_results'].to_csv(index=False).encode('utf-8')
        st.download_button(label="Download CSV",
                           data=csv,
                           file_name='enrichr_results.csv',
                           mime='text/csv')
    else:
        st.info("No enrichment results to display. Run a new analysis to see results.")


if __name__ == "__main__":
    main()