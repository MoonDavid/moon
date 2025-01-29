import streamlit as st
import gseapy as gp
import pandas as pd
from io import StringIO

import pandas as pd
import plotly.graph_objects as go
import plotly.colors
import numpy as np
import io
def create_enrichment_plot(df):
    """
    Generates an interactive enrichment plot with a custom legend
    based on Gene_set from a pandas DataFrame.

    Args:
      df: A pandas DataFrame containing enrichment data with columns:
        'Gene_set': Gene set category
        'Term': The specific term enriched
        'P-value': The p-value for the enrichment
        'Genes': A semicolon-separated string of enriched genes

    Returns:
      A plotly Figure object.
    """
    df["-log10(P-value)"] = -df["P-value"].apply(lambda p: p if p == 0 else np.log10(p))

    # Map each Gene_set to a unique color
    gene_sets = df["Gene_set"].unique()
    gene_set_colors = {
        gene_set: color
        for gene_set, color in zip(
            gene_sets, plotly.colors.qualitative.Plotly
        )
    }

    # Create the bar chart
    fig = go.Figure(
        go.Bar(
            y=df["Term"],
            x=df["-log10(P-value)"],
            orientation="h",
            marker=dict(color=df["Gene_set"].map(gene_set_colors)),
            hovertemplate=(
                "<b>Term:</b> %{y}<br>"
                "<b>Gene Set:</b> %{customdata[0]}<br>"
                "<b>-log10(P-value):</b> %{x:.2f}<br>"
                "<b>Genes:</b> %{customdata[1]}<extra></extra>"
            ),
            customdata=df[["Gene_set", "Genes"]].values,
            showlegend=False,  # Don't show default legend
        )
    )
    # Add custom legend entries as dummy scatter traces
    for gene_set, color in gene_set_colors.items():
      fig.add_trace(go.Scatter(
          x=[None],
          y=[None],
          mode='markers',
          marker=dict(color=color, size=10),
          name=gene_set,
          showlegend=True # Show the legend
      ))


    fig.update_layout(
        title="Enrichment Results",
        yaxis_title="Enriched Term",
        xaxis_title="-log10(P-value)",
        hovermode="closest",
        legend=dict(
            title="Gene Sets",
             orientation="h", # Horizontal legend
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
           )
    )

    return fig

def df_from_uniprots(df,uniprots):
    """Create a DataFrame from a list of UniProt IDs."""
    return df[df['UniProtKB-AC'].isin(uniprots)].reset_index(drop=True)
# Function to perform enrichment analysis
@st.cache_data
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
            st.session_state['enrichment_results_enrichr'] = results_df  # Save results in session state

            st.success("Enrichment analysis completed!")
        else:
            st.error("No results obtained or an error occurred during analysis.")
    # Display existing enrichment results if available
    if 'enrichment_results_enrichr' in st.session_state:
        st.subheader("Enrichment Results")
        st.dataframe(st.session_state['enrichment_results_enrichr'])

        # Download button for results
        csv = st.session_state['enrichment_results_enrichr'].to_csv(index=False).encode('utf-8')
        st.download_button(label="Download CSV",
                           data=csv,
                           file_name='enrichr_results.csv',
                           mime='text/csv')
        fig=create_enrichment_plot(st.session_state['enrichment_results_enrichr'])
        st.plotly_chart(fig)
    else:
        st.info("No enrichment results to display. Run a new analysis to see results.")


if __name__ == "__main__":
    main()