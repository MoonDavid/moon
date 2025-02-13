import streamlit as st
import gseapy as gp
import pandas as pd
from io import StringIO
import plotly.graph_objects as go
import plotly.colors
import numpy as np
import io
from gseapy import barplot, dotplot
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def create_enrichment_plot(df, top_n=10):
    """
    Generates an interactive enrichment plot with a custom legend
    based on Gene_set from a pandas DataFrame, ensuring top_n terms
    per Gene_set and grouping them accordingly.

    Args:
      df: A pandas DataFrame containing enrichment data
      top_n: Number of top terms to display per Gene_set
    """
    # Sort and filter top N terms per Gene_set
    df_sorted = df.sort_values("P-value")
    df_top = df_sorted.groupby("Gene_set").head(top_n)

    # Transform p-values to -log10 scale
    df_top["-log10(P-value)"] = -np.log10(df_top["P-value"])

    # Order by Gene_set first, then by -log10(P-value)
    df_top = df_top.sort_values(["Gene_set", "-log10(P-value)"], ascending=[True, False])

    gene_sets = df_top["Gene_set"].unique()
    gene_set_colors = {
        gene_set: color for gene_set, color in zip(gene_sets, plotly.colors.qualitative.Plotly)
    }

    fig = go.Figure()

    for gene_set in gene_sets:
        df_subset = df_top[df_top["Gene_set"] == gene_set]
        fig.add_trace(go.Bar(
            y=df_subset["Term"],
            x=df_subset["-log10(P-value)"],
            orientation="h",
            marker=dict(color=gene_set_colors[gene_set]),
            hovertemplate=(
                "<b>Term:</b> %{y}<br>"
                "<b>Gene Set:</b> %{customdata[0]}<br>"
                "<b>-log10(P-value):</b> %{x:.2f}<br>"
                "<b>Genes:</b> %{customdata[1]}<extra></extra>"
            ),
            customdata=df_subset[["Gene_set", "Genes"]].values,
            name=gene_set,
        ))

    fig.update_layout(
        title=f"Top {top_n} Enrichment Results per Gene Set",
        yaxis_title="Enriched Term",
        xaxis_title="-log10(P-value)",
        hovermode="closest",
        barmode="group",
        legend=dict(
            title="Gene Sets",
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )
    return fig
def create_gseapy_dotplot(results_df, top_n=10):
    """Create a GSEA dotplot using gseapy"""

    ax = dotplot(results_df,
                 column="Adjusted P-value",
                 x='Gene_set',
                 size=top_n*len(results_df['Gene_set'].unique()),
                 top_term=top_n,
                 title="Enrichment Dotplot",
                 xticklabels_rot=45,
                 show_ring=False,
                 marker='o',
                 figsize=(10, 20))
    return ax.figure


def create_gseapy_barplot(results_df, top_n=10):
    """Create a GSEA barplot using gseapy"""
    cmap = plt.cm.viridis # You can change this to tab20, Set1, Paired, etc.
    colors = cmap(np.linspace(0, 1, len(results_df['Gene_set'].unique()) + 1))
    colors= [mcolors.rgb2hex(c) for c in colors]

    ax = barplot(results_df,
                 column="Adjusted P-value",
                 group='Gene_set',
                 size=top_n*len(results_df['Gene_set'].unique()),
                 top_term=top_n,
                 color=colors,
                 title="Enrichment Barplot",
                 figsize=(10, 20))
    return ax.figure


def df_from_uniprots(df, uniprots):
    """Create a DataFrame from a list of UniProt IDs."""
    return df[df['UniProtKB-AC'].isin(uniprots)].reset_index(drop=True)


@st.cache_data
def perform_enrichment(gene_list, libraries, background_list=None):
    """Performs enrichment analysis using gseapy."""
    try:
        enr = gp.enrichr(gene_list=gene_list,
                         gene_sets=libraries,
                         background=background_list,
                         no_plot=True)
        return enr.results
    except Exception as e:
        st.error(f"Error during enrichment analysis: {e}")
        return None


def read_gene_list(text):
    if not text:
        return set()
    return {line.strip() for line in text.splitlines()}


def read_background(text):
    if not text:
        return None
    return {line.strip() for line in text.splitlines()}


def main():

    st.title("Enrichr Gene Set Enrichment Analysis")

    st.write("This tool uses [Enrichr](https://maayanlab.cloud/Enrichr/)  to perform gene set enrichment analysis. "
             "[Enrichr](https://maayanlab.cloud/Enrichr/) is a comprehensive gene set enrichment analysis tool and search engine that provides access to curated gene set libraries, advanced visualization features, and resources for biological discovery.")
    if 'gene_lists' not in st.session_state:
        st.warning(
            "No pre-loaded MPs list found in session state. Please load first main page Moonlighting Proteins Analysis to add some."
        )
        st.stop()
    available_gene_lists = {key: value for key, value in st.session_state['gene_lists'].items() if isinstance(value, list)}

    if 'gene_lists' not in st.session_state:
        st.session_state['gene_lists'] = {}

    if 'background_lists' not in st.session_state:
        st.session_state['background_lists'] = {}

    st.header("Gene List Input")
    with st.form(key="gene_list_form"):
        available_gene_lists = {key: value for key, value in st.session_state['gene_lists'].items() if isinstance(value, list)}
        if available_gene_lists:
            sorted_gene_list_keys = list(available_gene_lists.keys())
            selected_gene_list = st.selectbox(
                "Select the gene list to use:",
                options=sorted_gene_list_keys,
                key='select_gene_list_enrichr'
            )
        else:
            st.warning("No gene list available in session state.")
            selected_gene_list = None


        if "enrichr_libraries" not in st.session_state:
            available_libraries = gp.get_library_name()
            st.session_state["enrichr_libraries"] = available_libraries
        selected_libraries = st.multiselect("Select Enrichment Libraries",
                                            options=st.session_state["enrichr_libraries"],
                                            default=None,
                                            key='select_lib')
        submit_button = st.form_submit_button(label="Run Analysis")

    if submit_button:
        #st.write(st.session_state.keys())
        st.info(f"Enrichment analysis for '{selected_gene_list}'")
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
            results_df = perform_enrichment(list(selected_gene_list), selected_libraries)

        if results_df is not None and not results_df.empty:
            st.session_state['enrichment_results_enrichr'] = results_df
            st.success("Enrichment analysis completed!")
        else:
            st.error("No results obtained or an error occurred during analysis.")

    if 'enrichment_results_enrichr' in st.session_state:
        st.subheader("Enrichment Results")
        st.dataframe(st.session_state['enrichment_results_enrichr'])

        csv = st.session_state['enrichment_results_enrichr'].to_csv(index=False).encode('utf-8')
        st.download_button(label="Download CSV",
                           data=csv,
                           file_name='enrichr_results.csv',
                           mime='text/csv')

        # Add slider for selecting number of top genes
        top_n = st.slider("Select number of top terms to display",
                          min_value=5,
                          max_value=50,
                          value=10)

        # Create tabs for different plot types
        tab1, tab2, tab3 = st.tabs(["Interactive Plot", "Dotplot", "Barplot"])

        with tab1:
            fig = create_enrichment_plot(st.session_state['enrichment_results_enrichr'], top_n)
            st.plotly_chart(fig)

        with tab2:
            dot_fig = create_gseapy_dotplot(st.session_state['enrichment_results_enrichr'], top_n)
            st.pyplot(dot_fig)

        with tab3:
            bar_fig = create_gseapy_barplot(st.session_state['enrichment_results_enrichr'], top_n)
            st.pyplot(bar_fig)
    else:
        st.info("No enrichment results to display. Run a new analysis to see results.")


if __name__ == "__main__":
    main()