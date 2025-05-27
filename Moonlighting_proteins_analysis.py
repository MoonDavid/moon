import ast
import itertools
from collections import Counter
import altair as alt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import streamlit as st
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from matplotlib_venn import venn3
from pygments.lexer import default
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from scipy.stats import hypergeom, binomtest
from statsmodels.stats.multitest import multipletests
import altair as alt
import json
import plotly.express as px
import plotly.graph_objects as go
from pyvis.network import Network
import streamlit.components.v1 as components
import streamlit as st
import scipy.stats as stats
from scipy.stats import binomtest, chisquare, fisher_exact
from statsmodels.stats.proportion import proportions_ztest
import os






# =======================
# Helper Functions
# =======================

@st.cache_data
def load_data(file_path: str) -> pd.DataFrame:
    """Load the dataset from a CSV file."""
    df = pd.read_csv(file_path)
    df = df.replace("", np.nan)
    return df


def enrichment_analysis(query_uniprots, dict, correction='fdr_bh',N=None):
    """
    Perform over-representation analysis on a set of UniProt IDs ('query_uniprots')
    against a background defined by 'domain_dict' (domain -> list of UniProt IDs).
    Returns a DataFrame with enrichment results for each domain, sorted by corrected p-value.

    Parameters
    ----------
    query_uniprots : iterable of str
        List or set of UniProt IDs for which we want to test enrichment.
    domain_dict : dict
        Dictionary mapping domain (str) -> list of UniProt IDs (str) having that domain.
    correction : str
        Type of p-value correction (e.g., 'fdr_bh', 'bonferroni', 'holm', 'sidak', etc.).
        Passed to statsmodels.stats.multitest.multipletests.

    Returns
    -------
    pd.DataFrame
        Columns: ['domain', 'k', 'K', 'n', 'N', 'pval', 'pval_corrected'] sorted by 'pval_corrected'.
    """
    query_set = set(query_uniprots)

    # Build the set of all background UniProt IDs
    all_uniprots = set()
    for d, uniprot_list in dict.items():
        all_uniprots.update(uniprot_list)
    if N is None:
        N = len(all_uniprots)
    else:
        N = N

    n = len(query_set)  # size of the query

    terms = []
    pvals = []
    k_vals = []
    K_vals = []

    for term, uniprot_list in dict.items():
        K = len(uniprot_list)  # number of background proteins with this domain
        k = len(query_set.intersection(uniprot_list))  # overlap with query


        # If there's no overlap, p-value is 1 (not enriched),
        # but we still record it for completeness.
        if K == 0 or k == 0:
            pval = 1.0
        else:
            # hypergeom.sf(k-1, N, K, n) => P(X >= k)
            pval = hypergeom.sf(k - 1, N, K, n)

        terms.append(term)
        pvals.append(pval)
        k_vals.append(k)
        K_vals.append(K)

    # Multiple test correction
    reject, pvals_corrected, _, _ = multipletests(pvals, alpha=0.05, method=correction)

    # Create a DataFrame
    results_df = pd.DataFrame({
        'term': terms,
        'k': k_vals,
        'K': K_vals,
        'n': n,
        'N': N,
        'pval': pvals,
        'pval_corrected': pvals_corrected,
        'reject': reject
    })

    # Sort by corrected p-value ascending
    results_df.sort_values('pval_corrected', inplace=True)
    results_df.reset_index(drop=True, inplace=True)

    return results_df

def get_lists_from_session():
    """Retrieve keys from session_state that contain list objects."""
    lists = []
    for key, value in st.session_state['gene_lists'].items():
        if isinstance(value, list):
            lists.append(key)
    return lists
def df_from_uniprots(df,uniprots):
    """Create a DataFrame from a list of UniProt IDs."""
    return df[df['UniProtKB-AC'].isin(uniprots)].reset_index(drop=True)
def df_to_uniprots(df):
    """Convert a DataFrame to a list of UniProt IDs."""
    return df['UniProtKB-AC'].tolist()
def select_organisms(data: pd.DataFrame) -> list:
    """Allow user to select organisms of interest."""
    organisms = data['Organism'].unique()
    select_all = st.checkbox("Select all organisms", key='select_all_organisms')

    if select_all:
        selected = st.multiselect("Select organism", organisms, default=organisms, key='multiselect_all')
    else:
        selected = st.multiselect("Select organism", organisms, default=["Homo sapiens (Human)"],
                                  key='multiselect_specific')

    st.write(f"Selected organisms: {selected}")
    return selected


def plot_venn_diagram(data: pd.DataFrame):
    """Plot a Venn diagram showing the distribution of proteins across databases."""
    moondb_set = set(data[data['MoonDB'] == True]['UniProtKB-AC'])
    moonprot_set = set(data[data['MoonProtDB'] == True]['UniProtKB-AC'])
    multitask_set = set(data[data['MultiTaskProtDB'] == True]['UniProtKB-AC'])

    plt.figure(figsize=(3, 3))
    venn = venn3([moondb_set, moonprot_set, multitask_set],
                 set_labels=('MoonDB', 'MoonProt', 'MultiTaskProtDB'))
    for text in venn.set_labels:
        text.set_fontsize(8)
    plt.title("Venn diagram showing the distribution of human MPs across the 3 main databases",fontsize=8)
    st.pyplot(plt)


import pandas as pd
import streamlit as st


def filter_proteins(data: pd.DataFrame) -> pd.DataFrame:
    """
    Filter proteins based on selected databases and operation using a form layout.
    Stores results in session state with descriptive keys indicating the operation.
    Also initializes three default filters in session state.
    """

    # Initialize default filters if not already present in session state
    if 'default_filters_initialized' not in st.session_state:
        # Create the sets dictionary
        sets_dict = {
            "MoonDB": set(data[data['MoonDB'] == True]['UniProtKB-AC']),
            "MoonProt": set(data[data['MoonProtDB'] == True]['UniProtKB-AC']),
            "MultiTaskProtDB": set(data[data['MultiTaskProtDB'] == True]['UniProtKB-AC'])
        }

        # 1. Intersection of MoonProt and MultiTaskProtDB (Most Restrictive)
        intersection_set = sets_dict["MoonProt"].intersection(sets_dict["MultiTaskProtDB"])
        intersection_key = "humanMPs(MoonProt AND MultiTaskProtDB)"
        st.session_state['gene_lists'][intersection_key] = list(intersection_set)

        # 2. Union of MoonProt and MultiTaskProtDB
        union_two_set = sets_dict["MoonProt"].union(sets_dict["MultiTaskProtDB"])
        union_two_key = "humanMPs(MoonProt OR MultiTaskProtDB)"
        st.session_state['gene_lists'][union_two_key] = list(union_two_set)

        # 3. Union of all three databases (Least Restrictive)
        union_all_set = sets_dict["MoonDB"].union(sets_dict["MoonProt"]).union(sets_dict["MultiTaskProtDB"])
        union_all_key = "humanMPs_all"
        st.session_state['gene_lists'][union_all_key] = list(union_all_set)

        # Mark that default filters have been initialized
        st.session_state['default_filters_initialized'] = True

    # Create a form container for better organization and control
    with st.form(key='protein_filter_form'):
        # Create columns for horizontal layout of radio buttons
        col1, col2 = st.columns([3, 2])

        with col1:
            # Multiselect for database selection
            selected_sets = st.multiselect(
                "Select one or more sets:",
                options=["MoonDB", "MoonProt", "MultiTaskProtDB"],
                default=["MoonDB", "MoonProt", "MultiTaskProtDB"],
            )

        with col2:
            # Horizontal radio buttons using columns
            operation = st.radio(
                "Select the operation:",
                options=["Intersection", "Union"],  # Changed labels to be more intuitive
                index=1,
                key='filter_operation',
                horizontal=True  # This makes the radio buttons horizontal
            )

        # Form submit button
        submitted = st.form_submit_button("Apply Filters")

    if submitted and selected_sets:
        # Create the sets dictionary
        sets_dict = {
            "MoonDB": set(data[data['MoonDB'] == True]['UniProtKB-AC']),
            "MoonProt": set(data[data['MoonProtDB'] == True]['UniProtKB-AC']),
            "MultiTaskProtDB": set(data[data['MultiTaskProtDB'] == True]['UniProtKB-AC'])
        }

        # Get the selected sets' data
        selected_sets_data = [sets_dict[set_name] for set_name in selected_sets]

        # Perform set operation based on selection
        if operation == 'Intersection':
            filtered_set = set.intersection(*selected_sets_data)
            operation_str = " AND "
        else:  # Union operation
            filtered_set = set.union(*selected_sets_data)
            operation_str = " OR "

        # Filter the data
        filtered_data = data[data['UniProtKB-AC'].isin(filtered_set)].reset_index(drop=True)

        # Create a descriptive key for session state
        # Convert database names to title case and join them
        db_names = [name.title() for name in selected_sets]
        session_key = f"humanMPs({operation_str.join(db_names)})"
        #session_key2 = f"humanMPs({operation_str.join(db_names)})"

        # Store the filtered UniProtKB-AC list in session state
        st.session_state['gene_lists'][session_key] = filtered_data['UniProtKB-AC'].tolist()

        # Display results
        st.write(f"Number of filtered proteins: {filtered_data.shape[0]}")
        st.write(f"Results stored in session state with key: {session_key}")
        st.dataframe(filtered_data)


        return filtered_data

    elif submitted and not selected_sets:
        st.warning("Please select at least one database for filtering.")
        return data




def get_dataframes_from_session():
    dataframes = {}
    for key, value in st.session_state.items():
        if isinstance(value, pd.DataFrame):
            dataframes[key] = value
    return dataframes


def feature_analysis(filtered_data: pd.DataFrame, selected_df_name: str):
    """Analyze features related to Gene Ontology."""
    # st.subheader("Feature Analysis")

    # Filter columns related to Gene Ontology
    go_columns = [col for col in filtered_data.columns if "Gene Ontology" in col]
    column = st.selectbox("Select a feature", go_columns, key='feature_select')

    if column:
        if pd.api.types.is_numeric_dtype(filtered_data[column]):
            display_numeric_feature(filtered_data, column, selected_df_name)
        else:
            display_categorical_feature(filtered_data, column, selected_df_name)

@st.cache_data
def display_numeric_feature(data: pd.DataFrame, column: str, selected_df_name: str):
    """Display analysis for numeric features."""
    st.write(f"The selected column `{column}` is **continuous (numeric)**.")
    st.write("### Descriptive statistics:")
    st.write(data[column].describe())

    st.write(f"#### Histogram and boxplot of `{column}`:")
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    sns.histplot(data[column].dropna(), kde=True, bins=30, color='skyblue', ax=axes[0])
    axes[0].set_xlabel(column)
    axes[0].set_ylabel("Frequency")
    axes[0].set_title(f"Histogram of {column}")

    sns.boxplot(x=data[column].dropna(), color='lightgreen', ax=axes[1])
    axes[1].set_xlabel(column)
    axes[1].set_title(f"Boxplot of {column}")

    plt.tight_layout()
    st.pyplot(fig)


def display_categorical_feature(data: pd.DataFrame, column: str, selected_df_name: str):
    """Display analysis for categorical features with UniProtKB-AC tracking."""

    # Display basic information about the column
    st.write(f"The selected column `{column}` is **categorical** or treated as such.")
    unique_values = data[column].dropna().unique()

    # Count and plot top occurrences
    expanded_data = data[column].dropna().str.split('; ').explode()
    count_data = expanded_data.value_counts()

    if not count_data.empty:
        st.write(f"### Count of most common values in `{column}`:")
        count_df = count_data.head(20).reset_index()
        count_df.columns = [column, "Count"]

        # Display the dataframe with multi-row selection
        st.write("Select values from the table below to save associated UniProtKB-AC IDs for later analysis:")
        event = st.dataframe(
            count_df,
            use_container_width=True,
            hide_index=True,
            on_select="rerun",
            selection_mode="multi-row",
        )
        button = st.button("Save selection", key="save_selection")
        if event.selection and button:
            selected_rows = event.selection.rows  # Get selected row indices
            filtered_df = count_df.iloc[selected_rows]  # Filter dataframe
            selected_values = filtered_df[column].tolist()
            # Print in Streamlit
            #st.write(f"Selected values: {selected_values}")
            selected_uniprot_ids = data[
                data[column].apply(lambda x: any(val in str(x).split("; ") for val in selected_values))
            ]["UniProtKB-AC"].dropna().unique().tolist()

            #st.write(f"Associated UniProtKB-AC IDs: {selected_uniprot_ids}")
            session_key = f"{selected_df_name}-{column}:{selected_values}"
            st.session_state['gene_lists'][session_key] = selected_uniprot_ids
            st.success(f"Saved in session state with key `{session_key}`.")


    else:
        st.warning("No data available for visualization.")

    # Analyze Gene Ontology pairs if applicable
    if "Gene Ontology" in column:
        st.write(f"### Count of most common pairs in `{column}`:")
        go_terms = data[column].dropna().str.split('; ')
        go_terms = data[column].dropna().str.split('; ')
        # Split the GO terms into lists
        data['GO_terms_list'] = data[column].dropna().str.split('; ')

        # Function to generate GO term pairs
        def generate_go_pairs(go_terms):
            """
            Generate all unique pairs of GO terms from a list.

            Parameters:
            - go_terms (list): List of GO terms for a UniProt entry.

            Returns:
            - list of tuples: Each tuple contains a pair of GO terms.
            """
            # Ensure there are at least two GO terms to form a pair
            if not isinstance(go_terms, list):
                return []
            if len(go_terms) < 2:
                return []

            # Sort the GO terms to maintain consistency in pairs
            sorted_terms = sorted(go_terms)

            # Generate all unique combinations of two GO terms
            return list(itertools.combinations(sorted_terms, 2))

        # Apply the function to generate GO term pairs
        data['GO_term_pairs'] = data['GO_terms_list'].apply(generate_go_pairs)

        # Create the dictionary mapping UniProt to GO term pairs
        uniprot_to_go_pairs = pd.Series(data['GO_term_pairs'].values, index=data['UniProtKB-AC']).to_dict()
        go_pairs = data['GO_term_pairs'].explode().dropna().tolist()

        go_pair_counts = Counter(go_pairs)
        go_pair_counts_df = pd.DataFrame(go_pair_counts.items(), columns=["Pair", "Count"]).sort_values(by="Count",
                                                                                                        ascending=False).head(
            20)

        # Display the dataframe with multi-row selection
        st.write("Select rows from the table below to save associated UniProtKB-AC IDs for later analysis:")
        event = st.dataframe(
            go_pair_counts_df,
            use_container_width=True,
            hide_index=True,
            on_select="rerun",
            selection_mode="multi-row",
        )
        button2 = st.button("Save selection", key="save_selection2")
        if event.selection and button2:
            selected_rows = event.selection.rows
            selected_pairs = go_pair_counts_df.iloc[selected_rows]["Pair"].tolist()
            st.write(f"Selected pairs: {selected_pairs}")

            def find_uniprots_by_go_pair(uniprot_to_go_pairs, go_pair):
                return [uniprot for uniprot, pairs in uniprot_to_go_pairs.items() if go_pair in pairs]

            selected_uniprot_ids = []
            for pair in selected_pairs:
                selected_uniprot_ids.extend(find_uniprots_by_go_pair(uniprot_to_go_pairs, pair))

            selected_uniprot_ids = list(set(selected_uniprot_ids))  # Remove duplicates
            session_key = f"{selected_df_name}-{column}:{selected_pairs}"
            st.session_state['gene_lists'][session_key] = selected_uniprot_ids
            st.success(f"Saved in session state with key `{session_key}`.")



@st.cache_data
def sequence_analysis(filtered_data: pd.DataFrame):
    """Analyze protein sequences with extended physicochemical properties, with interactive plots and tabbed layout."""

    if 'Sequence' not in filtered_data.columns:
        st.error("The 'Sequence' column is missing in the DataFrame.")
        return

    st.write("## Protein Sequence Analysis")

    # 1. Sequence Length
    filtered_data['sequence_length'] = filtered_data['Sequence'].apply(lambda x: len(x) if pd.notnull(x) else 0)

    # --- Tabs for different plot types ---
    tabs = st.tabs([
        "Sequence Length",
        "Amino Acid Composition",
        "Physicochemical Properties",
        "Secondary Structure"
    ])


    with tabs[0]:
        st.write("### Protein Sequence Length:")
        st.write(filtered_data['sequence_length'].describe())

        fig = px.histogram(
            filtered_data,
            x='sequence_length',
            title="Distribution of Protein Sequence Lengths",
            labels={'sequence_length': 'Sequence Length', 'count': 'Frequency'},
            color_discrete_sequence=['salmon']
        )
        st.plotly_chart(fig)



    with tabs[1]:
        st.write("### Amino Acid Composition:")

        def get_amino_acids_percent(seq: str) -> dict:
            """Calculate amino acid composition percentage using ProtParam."""
            try:
                prot_param = ProteinAnalysis(seq)
                aa_percent = prot_param.get_amino_acids_percent()
                # Convert fractions to percentages
                aa_percent_scaled = {aa: percent * 100 for aa, percent in aa_percent.items()}
                return aa_percent_scaled
            except Exception as e:
                st.error(f"Error calculating amino acid composition: {e}")
                return {}

        aa_composition = filtered_data['Sequence'].dropna().apply(get_amino_acids_percent)
        aa_composition_df = pd.DataFrame(list(aa_composition)).fillna(0)
        mean_aa_composition = aa_composition_df.mean().sort_values(ascending=False)

        st.write("#### Average Amino Acid Composition (%):")
        st.dataframe(mean_aa_composition)
        fig_aa_bar = px.bar(
            x=mean_aa_composition.index,
            y=mean_aa_composition.values,
            title="Average Amino Acid Composition of Protein Sequences",
            labels={'x': "Amino Acid", 'y': "Percentage (%)"},
             color_discrete_sequence=px.colors.sequential.Magma
        )
        fig_aa_bar.update_xaxes(tickangle=90)
        st.plotly_chart(fig_aa_bar)


        st.write("### Amino Acid Frequency:")
        if not aa_composition_df.empty:
            fig_aa_heatmap = px.imshow(
                aa_composition_df.T,
                labels=dict(x="Proteins", y="Amino Acids", color="Percentage"),
                title="Heatmap of Amino Acid Composition",
                color_continuous_scale="YlGnBu"
            )
            st.plotly_chart(fig_aa_heatmap)
        else:
            st.write("No data to display Amino Acid Frequency Heatmap")


    with tabs[2]:
        st.write("### Physicochemical Properties of Protein Sequences:")

        def calculate_properties(seq: str) -> pd.Series:
                """Calculate various physicochemical properties of a protein sequence."""
                if pd.isnull(seq) or not seq:
                    return pd.Series({
                        "Molecular_Weight": None,
                        "Isoelectric_Point": None,
                        "Aromaticity": None,
                        "Instability_Index": None,
                        "Flexibility": None,
                        "Gravy": None,
                        "Molar_Extinction_Reduced": None,
                        "Molar_Extinction_Oxidized": None,
                        "Charge_pH7.0": None
                    })
                try:
                    prot_param = ProteinAnalysis(seq)
                    mw = prot_param.molecular_weight()
                    pI = prot_param.isoelectric_point()
                    aromaticity = prot_param.aromaticity()
                    instability = prot_param.instability_index()
                    flexibility = prot_param.flexibility()
                    gravy = prot_param.gravy()
                    epsilon_reduced, epsilon_oxidized = prot_param.molar_extinction_coefficient()
                    charge_pH7 = prot_param.charge_at_pH(7.0)

                    return pd.Series({
                        "Molecular_Weight": mw,
                        "Isoelectric_Point": pI,
                        "Aromaticity": aromaticity,
                        "Instability_Index": instability,
                        "Flexibility": flexibility,
                        "Gravy": gravy,
                        "Molar_Extinction_Reduced": epsilon_reduced,
                        "Molar_Extinction_Oxidized": epsilon_oxidized,
                        "Charge_pH7.0": charge_pH7
                    })
                except Exception as e:
                    #st.error(f"Error processing sequence: {e}")
                    return pd.Series({
                        "Molecular_Weight": None,
                        "Isoelectric_Point": None,
                        "Aromaticity": None,
                        "Instability_Index": None,
                        "Flexibility": None,
                        "Gravy": None,
                        "Molar_Extinction_Reduced": None,
                        "Molar_Extinction_Oxidized": None,
                        "Charge_pH7.0": None
                    })

        properties = filtered_data['Sequence'].apply(calculate_properties)
        filtered_data = pd.concat([filtered_data, properties], axis=1)


        # Helper function to safely plot histograms
        def safe_interactive_histplot(data, xlabel, title, color):
                data = data.dropna()
                if len(data) > 0:
                    fig = px.histogram(
                        x=data,
                        title=title,
                         labels={'x': xlabel, 'count': 'Frequency'},
                       color_discrete_sequence=[color]
                    )
                    st.plotly_chart(fig)
                else:
                    st.write(f"No data available for {title}.")


        # 3.1 Molecular Weight Statistics
        st.write("#### Molecular Weight Statistics:")
        st.write(filtered_data['Molecular_Weight'].describe())
        safe_interactive_histplot(
            filtered_data['Molecular_Weight'],
            xlabel="Molecular Weight (Da)",
            title="Distribution of Protein Molecular Weight",
            color='teal'
        )

        # 3.2 Isoelectric Point Statistics
        st.write("#### Isoelectric Point (pI) Statistics:")
        st.write(filtered_data['Isoelectric_Point'].describe())
        safe_interactive_histplot(
             filtered_data['Isoelectric_Point'],
             xlabel="Isoelectric Point (pI)",
            title="Distribution of Protein Isoelectric Point",
            color='orange'
        )

        # 3.3 Aromaticity Statistics
        st.write("#### Aromaticity Statistics:")
        st.write(filtered_data['Aromaticity'].describe())
        safe_interactive_histplot(
            filtered_data['Aromaticity'],
            xlabel="Aromaticity",
            title="Distribution of Protein Aromaticity",
            color='purple'
        )

        # 3.4 Instability Index Statistics
        st.write("#### Instability Index Statistics:")
        st.write(filtered_data['Instability_Index'].describe())
        safe_interactive_histplot(
            filtered_data['Instability_Index'],
            xlabel="Instability Index",
            title="Distribution of Protein Instability Index",
            color='darkred'
        )

        # Highlight unstable proteins (Instability Index > 40)
        unstable = filtered_data[filtered_data['Instability_Index'] > 40]
        st.write(f"**Number of Unstable Proteins (Instability Index > 40):** {unstable.shape[0]}")


        # 3.6 Gravy (Grand Average of Hydropathy) Statistics
        st.write("#### Gravy (Grand Average of Hydropathy) Statistics:")
        st.write(filtered_data['Gravy'].describe())
        safe_interactive_histplot(
            filtered_data['Gravy'],
            xlabel="Gravy (Hydropathy)",
            title="Distribution of Protein Gravy (Hydropathy)",
            color='gold'
        )

        # 3.7 Molar Extinction Coefficient Statistics
        st.write("#### Molar Extinction Coefficient (ε) Statistics:")
        st.write(filtered_data[['Molar_Extinction_Reduced', 'Molar_Extinction_Oxidized']].describe())

        # Plot Molar Extinction Coefficients
        epsilon_reduced = filtered_data['Molar_Extinction_Reduced'].dropna()
        epsilon_oxidized = filtered_data['Molar_Extinction_Oxidized'].dropna()

        if len(epsilon_reduced) > 0 or len(epsilon_oxidized) > 0:
                fig_epsilon = px.histogram(
                [epsilon_reduced,epsilon_oxidized],
                labels=dict(variable='Type',value="Molar Extinction Coefficient (ε)",count="Frequency"),
                title="Distribution of Molar Extinction Coefficient",
                color_discrete_sequence=['navy','magenta']
                )
                fig_epsilon.update_traces(opacity=0.75)
                st.plotly_chart(fig_epsilon)
        else:
            st.write("No data available to plot Molar Extinction Coefficient.")


        # 3.8 Charge at pH 7.0 Statistics
        st.write("#### Charge at pH 7.0 Statistics:")
        st.write(filtered_data['Charge_pH7.0'].describe())
        safe_interactive_histplot(
            filtered_data['Charge_pH7.0'],
            xlabel="Charge at pH 7.0",
            title="Distribution of Protein Charge at pH 7.0",
            color='brown'
        )


    with tabs[3]:
         # 4. Secondary Structure Analysis
        st.write("### Secondary Structure Analysis:")

        def calculate_secondary_structure(seq: str) -> pd.Series:
            """Calculate secondary structure fractions of a protein sequence."""
            if pd.isnull(seq) or not seq:
                return pd.Series({
                    "Fraction_Alpha_Helix": None,
                    "Fraction_Turn": None,
                    "Fraction_Beta_Sheet": None
                })
            try:
                prot_param = ProteinAnalysis(seq)
                ss_frac = prot_param.secondary_structure_fraction()
                # ProtParam returns (helix, turn, sheet)
                return pd.Series({
                    "Fraction_Alpha_Helix": ss_frac[0],
                    "Fraction_Turn": ss_frac[1],
                    "Fraction_Beta_Sheet": ss_frac[2]
                })
            except Exception as e:
                st.error(f"Error processing secondary structure for sequence: {e}")
                return pd.Series({
                    "Fraction_Alpha_Helix": None,
                    "Fraction_Turn": None,
                    "Fraction_Beta_Sheet": None
                })

        secondary_structure = filtered_data['Sequence'].apply(calculate_secondary_structure)
        filtered_data = pd.concat([filtered_data, secondary_structure], axis=1)

        # 4.1 Secondary Structure Statistics
        st.write("#### Secondary Structure Fractions:")
        st.write(filtered_data[['Fraction_Alpha_Helix', 'Fraction_Turn', 'Fraction_Beta_Sheet']].describe())

        # 4.2 Plot Average Secondary Structure Fractions
        ss_df = filtered_data[['Fraction_Alpha_Helix', 'Fraction_Turn', 'Fraction_Beta_Sheet']].dropna()
        if not ss_df.empty:
                fig_ss = px.box(
                    ss_df,
                    labels=dict(variable="Secondary Structure Type",value="Fraction"),
                    title="Distribution of Secondary Structure Fractions",
                    color_discrete_sequence=px.colors.qualitative.Vivid
                )
                st.plotly_chart(fig_ss)
        else:
            st.write("No data available to plot Secondary Structure Fractions.")


def safe_literal_eval(val):
    try:
        return ast.literal_eval(val)
    except (ValueError, SyntaxError):
        return val
def get_key_by_value(d, value):
    """Return the first key in dict d that has the given value."""
    for k, v in d.items():
        if v == value:
            return k
    return None

def interpro_analysis(data: pd.DataFrame, selected_df_name):
    """Analyze InterPro data."""
    st.write(
        "InterPro provides functional analysis of proteins by classifying them into families and predicting domains and important sites. To classify proteins in this way, InterPro uses predictive models, known as signatures, provided by several different databases (referred to as member databases) that make up the InterPro consortium:"
        " CATH, CDD, HAMAP, MobiDB Lite, Panther, Pfam, PIRSF, PRINTS, Prosite, SFLD, SMART, SUPERFAMILY AND NCBIfam. The [InterPro Consortium](https://interpro-documentation.readthedocs.io/en/latest/databases.html) section gives further information about the individual databases.")

    # Function to safely evaluate strings to Python literals

    abbr= {
        "IP": "InterPro",
        "TI": "ncbifam",
        "PT": "panther",
        "PS": "prosite",
        "NF": "ncbifam",
        "PR": "prints",
        "PI": "pirsf",
        "cd": "cdd",
        "SS": "ssf",
        "G3": "cathgene3d",
        "MF": "hamap",
        "SM": "smart",
        "SF": "sfld",
        "PF": "pfam",

    }
    # Convert the 'interpro_data' column from string to actual tuples

    st.write("### Enrichment of protein domains:")
    #read dict from domain_dict.json
    with st.expander("Statistical Methods Explanation"):
        st.markdown("#### Hypergeometric Test")
        st.markdown(
            "The hypergeometric test calculates the probability of observing at least $k$ proteins with a particular domain in a sample of size $n$, given that there are $K$ proteins with that domain in a population of size $N$."
        )

        st.latex(
            r"""P(X \geq k) = \sum_{i=k}^{\min(n,K)} \frac{{K \choose i}{{N-K} \choose {n-i}}}{{N \choose n}}"""
        )

        st.markdown("Where:")
        st.markdown("- $N$ = Total number of proteins in the background set")
        st.markdown("- $K$ = Number of proteins in the background set with the domain")
        st.markdown("- $n$ = Number of proteins in the query set")
        st.markdown("- $k$ = Number of proteins in the query set with the domain")
        st.markdown(
            "- $P(X \\geq k)$ = Probability of observing at least $k$ proteins with the domain by chance"
        )

        st.markdown("#### Benjamini-Hochberg Correction")
        st.markdown(
            "The Benjamini-Hochberg procedure controls the false discovery rate (FDR) when performing multiple hypothesis tests:"
        )

        st.markdown(
            "1. Sort all p-values in ascending order: $p_1 \\leq p_2 \\leq ... \\leq p_m$"
        )
        st.markdown(
            "2. For a given FDR level $\\alpha$ (typically 0.05), find the largest $k$ such that:"
        )

        st.latex(r"""p_k \leq \frac{k}{m} \cdot \alpha""")

        st.markdown(
            "3. Reject the null hypothesis (i.e., consider the enrichment significant) for all tests with p-values $p_1, p_2, ..., p_k$"
        )
        st.markdown("4. The adjusted p-value for each test is calculated as:")

        st.latex(
            r"""p_{adj,i} = \min\left(\min_{j \geq i}\left(\frac{m \cdot p_j}{j}\right), 1\right)"""
        )

        st.markdown(
            "This method controls the expected proportion of false positives among all rejected hypotheses."
        )
    #domain_dict = json.load(open("data/domain_dict.json", "r"))
    db_options = abbr.values()
    selected_db = st.selectbox(
        "Select InterPro member database to display the enrichment:",
        options=db_options,
        index=0
    )
    #delete all keys that starts with DP
    domain_dict = json.load(open("data/domain_dict.json", "r"))
    domain_dict = {k: v for k, v in domain_dict.items() if k.startswith(get_key_by_value(abbr, selected_db))}
    results_df=enrichment_analysis(dict=domain_dict, query_uniprots=df_to_uniprots(data), correction='fdr_bh')
    st.subheader(f"Enrichment Results for {selected_db} domains")

    # Add a column for negative log10 p-value (for plotting)
    #results_df["-log10(pval_corrected)"] = -results_df["pval_corrected"].apply(lambda x: 1e-300 if x <= 0 else x).apply(np.log10)
    results_df['weblink'] = results_df['term'].apply(
        lambda x: f"https://www.ebi.ac.uk/interpro/entry/{abbr[x[:2]]}/{x}" if x[1].isalpha() else f"https://www.ebi.ac.uk/interpro/entry/cathgene3d/G3DSA:{x}")
    results_df['database'] = results_df['term'].apply(
        lambda x: x[:2] if x[1].isalpha() else "CATH"
    )



    st.dataframe(
        results_df,
        column_config={
            "weblink": st.column_config.LinkColumn(
                "weblink",
                help="Link to InterPro entry"
            )
        },
        hide_index=True
    )

@st.cache_data
def rna_binding_analysis(data: pd.DataFrame, selected_df_name):
    st.write(
        "[RBPWorld](http://research.gzsys.org.cn/eurbpdb2/index.html) is an updated version of EuRBPDB, specifically designed to"
        " unveil the functions and disease associations of RNA-binding proteins (RBPs)"
        " with heightened efficacy. Within RBPWorld, an expansive collection of 1,393,686 RBPs"
        " across 445 species, including 3,303 human RBPs (hRBPs). ")
    st.write(f"There are {data['RBP type'].notna().sum()} RBPs annotated in RBPWorld in this dataset.")

    # Filter and display relevant RBP data
    datarbp = data[['UniProtKB-AC', 'Gene symbol', 'RBP type', 'No. RBPome']]
    datarbp = datarbp[datarbp['RBP type'].notna()]
    datarbp = datarbp.reset_index(drop=True)
    st.dataframe(datarbp)

    # Button to save RBPs to session state
    if st.button("Save to session state"):
        sessionkey = f"RBPs_{selected_df_name}"
        if sessionkey not in st.session_state:
            st.session_state['gene_lists'][sessionkey] = df_to_uniprots(datarbp)
            st.success(f"RBPs saved to session state with key: {sessionkey}")

    # RBP Type Distribution
    st.write("### RBP Type Distribution:")
    rbp_type_counts = data['RBP type'].value_counts().reset_index()
    rbp_type_counts.columns = ['RBP Type', 'Count']
    st.dataframe(rbp_type_counts)

    # Plot RBP Type Distribution
    st.write("#### RBP Type Distribution Plot")
    plt.figure(figsize=(10, 6))
    sns.barplot(data=rbp_type_counts, x='RBP Type', y='Count', palette='viridis')
    plt.xticks(rotation=45, ha='right')
    plt.xlabel('RBP Type')
    plt.ylabel('Number of RBPs')
    plt.title('Distribution of RBP Types in RBPWorld')
    plt.tight_layout()
    st.pyplot(plt)
    plt.clf()  # Clear the figure after plotting

    # Distribution of "No. RBPome"
    st.write("### Distribution of 'No. RBPome':")
    st.dataframe(data['No. RBPome'].describe())

    # Plot Distribution of "No. RBPome"
    st.write("#### 'No. RBPome' Distribution Plot")
    plt.figure(figsize=(10, 6))
    sns.histplot(data['No. RBPome'].dropna(), bins=30, kde=True, color='skyblue')
    plt.xlabel('Number of RBPome')
    plt.ylabel('Frequency')
    plt.title('Distribution of "No. RBPome"')
    plt.tight_layout()
    st.pyplot(plt)
    plt.clf()  # Clear the figure after plotting

@st.cache_data
def disorder_analysis(data: pd.DataFrame, selected_df_name: str):
    """
    Analyzes disorder-related columns in the given DataFrame and visualizes the data.

    Parameters:
    - data (pd.DataFrame): The input DataFrame containing disorder information.
    - selected_df_name (str): The name of the selected DataFrame for display purposes.
    """
    st.write(
        "Intrinsically disordered proteins (IDPs) or intrinsically disordered regions (IDRs) are proteins or regions"
        " of proteins that lack a fixed or ordered three-dimensional structure. IDPs are highly flexible and can"
        " adopt different conformations in different contexts. They play crucial roles in various biological processes,"
        " including signaling, transcription, and cell cycle regulation. Here, we analyze the disorder content of MPs in the selected dataset."
        "  [DisProt](https://www.disprot.org/) is a database of proteins with experimentally verified disordered regions, while AlphaDB and MobiDB Lite are"
        " databases that predict disordered regions in proteins.")
    # Select relevant columns for disorder analysis
    disorder_columns = [
        'UniProtKB-AC',
        'disprot_id',
        'disprot_disorder',
        'curated_disorder',  # Updated column
        'alphadb_disorder',
        'mobidblite_disorder'
    ]

    # Ensure all specified columns exist in the DataFrame
    missing_cols = set(disorder_columns) - set(data.columns)
    if missing_cols:
        st.error(f"The following required columns are missing from the data: {', '.join(missing_cols)}")
        return

    dataidp = data[disorder_columns].copy()
    dataidp['disprot_disorder'] = dataidp['disprot_disorder'].apply(
        lambda x: x / 100 if pd.notnull(x) else x
    )

    # Display the DataFrame
    st.subheader("Data Preview")
    st.dataframe(dataidp)

    # Calculate and display the number of NaNs per column


    # Identify numerical disorder columns for plotting
    numerical_cols = [
        'disprot_disorder',
        'curated_disorder',
        'alphadb_disorder',
        'mobidblite_disorder'
    ]
    idp_dict = json.load(open("data/IDP_ontology.json", "r"))
    # delete all keys that starts with DP
    idp_dict = {k: v for k, v in idp_dict.items() if not k.startswith("DP")}
    results_df = enrichment_analysis(dict=idp_dict, query_uniprots=df_to_uniprots(data), correction='fdr_bh')
    st.subheader("Enrichment Results")
    st.write("The table below shows the results of the enrichment analysis using IDP ontology terms.")
    st.dataframe(results_df)
    # Check if numerical columns exist
    disorder_numerical_cols = [col for col in numerical_cols if col in dataidp.columns]
    st.subheader("Disorder Types Box Plot")

    # Melt the DataFrame to long format for Plotly
    boxplot_data = dataidp.melt(
        id_vars=['UniProtKB-AC', 'disprot_id'],
        value_vars=disorder_numerical_cols,
        var_name='Disorder Type',
        value_name='Normalized Disorder Value'
    )

    # Remove rows with NaN in 'Normalized Disorder Value'
    boxplot_data = boxplot_data.dropna(subset=['Normalized Disorder Value'])

    # Create the interactive box plot using Plotly Express
    fig = px.box(
        boxplot_data,
        x='Disorder Type',
        y='Normalized Disorder Value',
        points='all',  # Show all points
        title='Disorder region percentage',
        hover_data=['UniProtKB-AC', 'disprot_id', 'Normalized Disorder Value'],
        labels={
            'Disorder Type': 'Disorder Type',
            'Normalized Disorder Value': 'Normalized Disorder Value (0-1)'
        },
        color='Disorder Type'  # Different colors for each disorder type
    )

    # Update layout for better readability
    fig.update_layout(
        boxmode='group',
        xaxis_title='Disorder Type',
        yaxis_title='Disorder percentage (0-1)',
        legend_title='Disorder Type',
        title_x=0.5  # Center the title
    )
    st.plotly_chart(fig, use_container_width=True)

    st.subheader("Disorder Range Analysis")


    # Plot distribution of ranges
    range_data = pd.DataFrame()
    for col in numerical_cols:
        dataidp[f'{col}_range'] = pd.cut(
            dataidp[col],
            bins=[0, 0.2, 0.4, 0.6, 0.8, 1],
            labels=['0-20%', '20-40%', '40-60%', '60-80%', '80-100%']
        )

    # Plot distribution of ranges
    range_data = []
    for col in numerical_cols:
        range_counts = dataidp[f'{col}_range'].value_counts()
        for range_name, count in range_counts.items():
            range_data.append({
                'Range': range_name,
                'Count': count,
                'Disorder Type': col
            })
    range_data = pd.DataFrame(range_data)

    fig_ranges = px.bar(
        range_data,
        x='Range',
        y='Count',
        color='Disorder Type',
        barmode='group',
        title='Distribution of Disorder Ranges',
        category_orders={'Range': ['0-20%', '20-40%', '40-60%', '60-80%', '80-100%']}
    )
    st.plotly_chart(fig_ranges, use_container_width=True)
    # Optionally, display correlation matrix
    st.subheader("Correlation Matrix")
    if len(disorder_numerical_cols) < 2:
        st.warning("At least two numerical disorder columns are required to display a correlation matrix.")
    else:
        corr = dataidp[disorder_numerical_cols].corr()

        # Create the interactive heatmap using Plotly
        fig = px.imshow(corr,
                        x=corr.columns,
                        y=corr.index,
                        color_continuous_scale='RdBu',
                        text_auto=True
                        )

        # Remove axis names
        fig.update_xaxes(title_text='')
        fig.update_yaxes(title_text='')

        # Make the chart bigger
        fig.update_layout(
            title='Correlation Matrix of Disorder Columns',
            height=600,  # Adjust height as needed
            width=800,  # Adjust width as needed
            margin=dict(l=20, r=20, t=40, b=20)  # Adjust margins as needed

        )

        st.plotly_chart(fig, use_container_width=True)

def plot_disease_classifications(df, selected_column,topn):
    # Create DataFrame with gene IDs and exploded classificatio
    expanded_df = df[["gene_symbol", selected_column]].copy()

    # Explode and then group by geneid to get sets of terms for each gene
    expanded_df[selected_column] = expanded_df[selected_column].str.split("|")
    expanded_df = expanded_df.explode(selected_column)

    # Group by geneid and create sets of terms
    grouped_terms = (
        expanded_df.groupby("gene_symbol")[selected_column]
        .agg(set)  # Convert to set to get unique terms
        .reset_index()
    )

    # Now explode these sets to count unique term frequencies
    exploded_unique = grouped_terms[selected_column].explode()
    value_counts = exploded_unique.value_counts().reset_index()
    value_counts.columns = ["Category", "Unique Genes Count"]

    # Create an interactive bar plot
    fig = px.bar(
        value_counts.head(topn),
        x="Category",
        y="Unique Genes Count",
        title=f"Top 20 {selected_column} Distribution",
        labels={
            "Unique Genes Count": "Number of associated MPs",
            "Category": selected_column,
        },
    )

    # Customize the layout
    fig.update_layout(
        xaxis_tickangle=-45,
        height=600,
        showlegend=False,
        title_x=0.5,
        title_font_size=20,
    )

    return fig, value_counts


def create_gene_disease_network(disease_df, classification_column):
    # Create a network
    net = Network(height="750px", width="100%", bgcolor="#ffffff", font_color="black")

    # Add disease nodes
    diseases = disease_df[classification_column].dropna().unique()
    for disease in diseases:
        net.add_node(disease, label=disease, color="#ff7f7f", title=disease)

    # Add gene nodes and edges
    for _, row in disease_df.dropna(
        subset=[classification_column, "gene_symbol"]
    ).iterrows():
        disease = row[classification_column]
        gene = row["gene_symbol"]

        # Add gene node if it doesn't exist
        net.add_node(gene, label=gene, color="#7f7fff", title=gene)

        # Add edge between disease and gene
        net.add_edge(disease, gene)

    net.show_buttons(filter_=[ "physics"])
    # Save the network
    net.save_graph("gene_disease_network.html")
    try:

        with open("gene_disease_network.html", "r", encoding="utf-8") as file:
            html_content = file.read()
        components.html(html_content, height=1600, scrolling=True)

        # Provide a download button for the HTML
        with open("gene_disease_network.html", "rb") as file:
            st.download_button(
                label="Download Network Visualization",
                data=file,
                file_name="gene_disease_network.html",
                mime="text/html",
            )
    except Exception as viz_error:
        st.error(f"Visualization Error: {viz_error}")



def disease(data: pd.DataFrame, selected_df_name: str):
    # --- Title & Introduction ---
    st.markdown("## Disease Association Analysis")

    # removed the expander
    st.markdown("### About OMIM/MIM")
    st.write(
        """
        **OMIM (Online Mendelian Inheritance in Man)** is a comprehensive, authoritative compendium of human genes and genetic phenotypes.
        The database is continuously updated with information on all known Mendelian disorders and the genes involved in their expression.
        """
    )
    # --- Check for mim_morbid_accession column ---
    if 'mim_morbid_accession' not in data.columns:
        st.error("The column 'mim_morbid_accession' was not found in the provided dataframe.")
        return
    disease_protein_count = data['mim_morbid_accession'].notna().sum()

    # --- Constants and Calculations ---
    total_proteins = 21792
    disease_genes_total = 5310
    enrichment_percentage = (disease_genes_total / total_proteins) * 100
    sample_size = len(data)
    dataset_disease_percentage = (disease_protein_count / sample_size) * 100 if sample_size > 0 else 0

    # --- Display Enrichment Info ---
    st.markdown("#### Enrichment Overview")
    col1, col2 = st.columns(2)
    with col1:
        st.metric("OMIM-associated proteins in Swissprot", f"{disease_genes_total} / {total_proteins}",
                  f"{enrichment_percentage:.2f}%")
    with col2:
        st.metric(f"OMIM-associated proteins in {selected_df_name}", f"{disease_protein_count} / {sample_size}",
                  f"{dataset_disease_percentage:.2f}%")

    # Calculate baseline probability
    baseline_probability = disease_genes_total / total_proteins

    # Choose the test to perform
    selected_test = st.selectbox(
        "Perform a statistical test",
        [
            "Chi-square Goodness-of-Fit Test",
            "Binomial Test",
            "Fisher's Exact Test",
        ],
    )

    if selected_test == "Chi-square Goodness-of-Fit Test":
        st.markdown("#### Chi-square Goodness-of-Fit Test")
        # Observed counts: [disease-associated, not disease-associated] in the sample
        observed = [disease_protein_count, sample_size - disease_protein_count]
        # Expected counts based on the baseline probability
        expected = [
            sample_size * baseline_probability,
            sample_size * (1 - baseline_probability),
        ]
        chi2_stat, p_value_chi2 = chisquare(f_obs=observed, f_exp=expected)
        st.write(f"**Chi-square statistic:** {chi2_stat:.4f}")
        st.write(f"**P-value:** {p_value_chi2:.4e}")
        if p_value_chi2 < 0.05:
            st.success(
                "The observed proportion significantly deviates from the expected proportion (p < 0.05)."
            )
        else:
            st.warning(
                "The observed proportion is not significantly different from the expected proportion (p ≥ 0.05)."
            )

    elif selected_test == "Binomial Test":
        st.markdown("#### Binomial Test")
        binom_result = binomtest(
            disease_protein_count,
            n=sample_size,
            p=baseline_probability,
            alternative="two-sided",
        )
        p_value_binomial = binom_result.pvalue
        st.write(f"**P-value:** {p_value_binomial:.4e}")
        if p_value_binomial < 0.05:
            st.success(
                "The observed proportion significantly deviates from the expected proportion (p < 0.05)."
            )
        else:
            st.warning(
                "The observed proportion is not significantly different from the expected proportion (p ≥ 0.05)."
            )


    elif selected_test == "Fisher's Exact Test":
        st.markdown("#### Fisher's Exact Test")
        # Construct the 2x2 contingency table
        # Rows: [Sample, Background - Sample]
        # Columns: [Disease-associated, Not disease-associated]
        a = disease_protein_count
        b = sample_size - disease_protein_count
        c = disease_genes_total - disease_protein_count
        d = (total_proteins - sample_size) - (
            disease_genes_total - disease_protein_count
        )
        table = [[a, b], [c, d]]
        st.write("Contingency Table:")
        st.write(f"[[{a}, {b}],")
        st.write(f" [{c}, {d}]]")
        oddsratio, p_value_fisher = fisher_exact(table, alternative="two-sided")
        st.write(f"**Odds Ratio:** {oddsratio:.4f}")
        st.write(f"**P-value:** {p_value_fisher:.4e}")
        if p_value_fisher < 0.05:
            st.success("The association is statistically significant (p < 0.05).")
        else:
            st.warning("The association is not statistically significant (p ≥ 0.05).")

    # --- Plot ---
    st.markdown("#### Visualization")

    # Prepare data for the bar chart
    labels = ['UniProtKB with Disease Association', f'{selected_df_name} with Disease Association']
    values = [enrichment_percentage, dataset_disease_percentage]

    # Create a bar chart using Plotly
    fig = go.Figure(data=[go.Bar(x=labels, y=values, marker_color=['#3498db', '#e74c3c'])])
    fig.update_layout(
        title_text=f'OMIM Disease Association Proportion',
        yaxis_title='Percentage',
        yaxis=dict(range=[0, max(values) + 5]),  # Add buffer to Y axis
        xaxis_title='Dataset'
    )
    st.plotly_chart(fig)

    st.write('## About DIsGeNET')
    st.write("DisGeNET is a comprehensive knowledge platform for exploring the genetic basis of human diseases. It integrates data from multiple sources, including curated datasets, genome-wide association studies (GWAS), animal models, and the scientific literature."
             " Curated search involves accessing manually reviewed and validated data from various specialized databases: **CLINGEN**, **CLINVAR**, **PSYGENET**, **ORPHANET**, **UNIPROT**, **MGD (Human)**, and **RGD (Human**"
             )
    disease = pd.read_csv("data/gene_disease_associations.tsv", sep="\t")
    associated=disease[disease["gene_symbol"].isin(data["Entrez"])]["geneid"].nunique()
    st.metric(label="Number of associated genes in DIsGeNET", value=f"{associated}/{len(data)}")
    disease=disease[disease["gene_symbol"].isin(data["Entrez"])]
    st.write('### Most common diseases')
    st.dataframe(disease['disease_name'].value_counts())
    classification_columns = [
        "diseaseUMLSCUI",
        "diseaseClasses_MSH",
        "diseaseClasses_UMLS_ST",
        "diseaseClasses_DO",
        "diseaseClasses_HPO",
    ]


    selected_classification = st.selectbox(
        "Select Disease Classification System:", classification_columns
    )

    # Add a number input for top N categories to show
    top_n = st.slider("Select number of top categories to display:", 5, 50, 20)

    # Create two columns for layout
    col1, col2 = st.columns([2, 1])

    # Generate and display the plot
    fig, value_counts = plot_disease_classifications(disease, selected_classification,top_n)

    with col1:
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.write("### All Categories")
        st.dataframe(value_counts, height=400)

    # Add some summary statistics
    st.write("### Summary Statistics")
    total_categories = len(value_counts)
    total_entries = value_counts["Unique Genes Count"].sum()

    # Create three columns for metrics
    metric_col1, metric_col2, metric_col3 = st.columns(3)

    with metric_col1:
        st.metric("Total Categories", total_categories)
    with metric_col2:
        st.metric("Total Entries", total_entries)
    with metric_col3:
        st.metric(
            "Average Entries per Category", f"{total_entries / total_categories:.2f}"
        )
    if st.button("Generate Network Visualization"):
        create_gene_disease_network(disease, selected_classification)



def enzymes(data: pd.DataFrame, selected_df_name: str):
    import plotly.graph_objects as go
    import plotly.express as px

    # EC class descriptions
    ec_descriptions = {
        "1": "Oxidoreductases: Transfer of H and O atoms or electrons between substances",
        "2": "Transferases: Transfer of functional groups between substances",
        "3": "Hydrolases: Formation of products by hydrolysis",
        "4": "Lyases: Non-hydrolytic addition or removal of groups",
        "5": "Isomerases: Intramolecular rearrangement",
        "6": "Ligases: Join molecules with ATP breakdown",
        "7": "Translocases: Movement of molecules across membranes",
    }

    # Define custom colors for consistency across plots
    colors = [
        "#ff9999",
        "#66b3ff",
        "#99ff99",
        "#ffcc99",
        "#ff99cc",
        "#99ccff",
        "#ffff99",
    ]

    enzymes_df = data[data["EC number"].notna()].copy()

    # Main metrics

    enzyme_count = enzymes_df.shape[0]
    total_count = data.shape[0]
    percentage = (enzyme_count / total_count) * 100

    st.metric(
        label=f"Enzymes in {selected_df_name}",
        value=f"{enzyme_count}/{total_count}"
    )

    # Process EC numbers
    enzymes_df["EC_exploded"] = enzymes_df["EC number"].str.split("; ")
    enzymes_exploded = enzymes_df.explode("EC_exploded")
    enzymes_exploded["EC_class"] = enzymes_exploded["EC_exploded"].str.split(".").str[0]

    protein_ec_classes = enzymes_exploded.groupby("UniProtKB-AC")[
        "EC_class"
    ].nunique()

    enzymes_exploded = enzymes_exploded["EC_exploded"]

    # EC Classes data preparation
    ec_classes = enzymes_exploded.str.split(".").str[0]
    class_counts = ec_classes.value_counts().sort_index()



    # Create legend labels
    legend_labels = [
        f"EC {i}: {ec_descriptions[i].split(':')[0]}" for i in class_counts.index
    ]

    st.write("### EC Classes Analysis", unsafe_allow_html=True)
    st.write(list(class_counts))
    print(list(class_counts))
    # Pie Chart using Plotly
    fig_pie = go.Figure(
        data=[
            go.Pie(
                labels=legend_labels,
                values=class_counts.values,
                hovertemplate="<b>%{label}</b><br>"
                + "Count: %{value}<br>"
                + "<b>Percentage: %{percent}</b><extra></extra>",
                textinfo="percent",
                marker=dict(colors=colors[: len(class_counts)]),
            )
        ]
    )

    fig_pie.update_layout(
        title="EC Classes Distribution (%)",
        width=1200,
        height=600,
        showlegend=True,
        legend=dict(yanchor="top", y=1, xanchor="left", x=1,font=dict(size=24)),
    )

    st.plotly_chart(fig_pie, use_container_width=True)

    st.write("### Enzymes with Multiple EC Classes")

    ec_class_distribution = protein_ec_classes.value_counts().reset_index()
    proteins_multi_classes = protein_ec_classes[protein_ec_classes > 1]
    ec_class_distribution.columns = ["Number of EC Classes", "Protein Count"]
    st.metric(
        label=f"Enzymes with Multiple EC Classes in {selected_df_name}",
        value=f"{len(proteins_multi_classes)}/{enzyme_count}",    )
    # Step 7: Create a pie chart using Plotly Express
    fig = px.pie(
        ec_class_distribution,
        names="Number of EC Classes",
        values="Protein Count",
        title="Distribution of Proteins by Number of EC Classes",
    )
    st.plotly_chart(fig, use_container_width=True)


    selected_proteins = data[data["UniProtKB-AC"].isin(proteins_multi_classes.index)]
#ADD Column from protein_ec_classes in data
    selected_proteins["Count_diff_EC_classes"] = selected_proteins["UniProtKB-AC"].map(protein_ec_classes)
    #reorder columns with   EC number first
    selected_proteins = selected_proteins[["UniProtKB-AC", "EC number","Count_diff_EC_classes"] + [col for col in selected_proteins.columns if col not in ["UniProtKB-AC", "EC number","Count_diff_EC_classes"]]]
    st.dataframe(selected_proteins)

    st.write("### EC Subclasses Analysis", unsafe_allow_html=True)
    ec_subclasses = enzymes_exploded.str.split(".").str[:2].str.join(".")
    subclass_counts = ec_subclasses.value_counts().head(10)  # Show top 10 subclasses

    # Interactive bar chart for subclasses
    fig_subclass = go.Figure()
    fig_subclass.add_trace(
        go.Bar(
            x=subclass_counts.index,
            y=subclass_counts.values,
            text=subclass_counts.values,
            textposition="auto",
            marker_color=colors[: len(subclass_counts)],
            hovertemplate="<b>%{x}</b><br>" + "Count: %{y}<extra></extra>",
        )
    )

    fig_subclass.update_layout(
        title="Top 10 EC Subclasses",
        xaxis_title="EC Subclass",
        yaxis_title="Count",
        width=1200,
        height=600,
        hoverlabel=dict(bgcolor="white"),
        xaxis={"categoryorder": "total descending"},
    )

    st.plotly_chart(fig_subclass, use_container_width=True)

    # Save enzyme gene list functionality
    st.write("## 💾 Save Enzyme Gene List", unsafe_allow_html=True)
    if st.button("Save Enzyme Gene List"):
        enzyme_genes = enzymes_df["UniProtKB-AC"].tolist()
        session_key = f"Enzymes_{selected_df_name}"
        if session_key not in st.session_state["gene_lists"]:
            st.session_state["gene_lists"][session_key] = enzyme_genes
            st.success(
                f"✅ Enzyme UniProtKB-AC list saved in session state with key: {session_key}"
            )
        else:
            st.warning("⚠️ Gene list already exists in session state.")

@st.cache_data
def load_hpa_rna_data(filename='rna_consensus_tissue.tsv'):
    """Loads specified HPA RNA data file, handling different column names."""
    # Adjust the path as needed for your project structure
    file_path = os.path.join("data", filename)
    if not os.path.exists(file_path):
        st.error(f"HPA data file not found: {file_path}")
        return None

    try:
        df = None
        if filename == 'rna_single_cell_type.tsv':
            # Columns for single-cell data: Gene, Gene name, Cell type, nTPM
            required_cols = ['Gene', 'Gene name', 'Cell type', 'nTPM']
            col_dtypes = {'Gene': str, 'Gene name': str, 'Cell type': str, 'nTPM': str} # Read nTPM as string first for robust parsing
            #st.info(f"Loading single-cell data ({filename}) with columns: {required_cols}")
            df = pd.read_csv(
                file_path,
                sep='\t',
                usecols=required_cols,
                dtype=col_dtypes
            )
            # Convert nTPM to numeric, coercing errors
            df['nTPM'] = pd.to_numeric(df['nTPM'], errors='coerce')
            df.dropna(subset=['nTPM'], inplace=True) # Remove rows where nTPM couldn't be parsed

        elif filename == 'rna_tissue_consensus.tsv':
             # Columns for bulk tissue data: Gene, Gene name, Tissue, nTPM
            required_cols = ['Gene', 'Gene name', 'Tissue', 'nTPM']
            col_dtypes = {'Gene': str, 'Gene name': str, 'Tissue': str, 'nTPM': str} # Read nTPM as string first
            #st.info(f"Loading bulk tissue data ({filename}) with columns: {required_cols}")
            df = pd.read_csv(
                file_path,
                sep='\t',
                usecols=required_cols,
                 dtype=col_dtypes
            )
             # Convert nTPM to numeric, coercing errors
            df['nTPM'] = pd.to_numeric(df['nTPM'], errors='coerce')
            df.dropna(subset=['nTPM'], inplace=True) # Remove rows where nTPM couldn't be parsed

        else:
            st.error(f"Unknown HPA data file type specified: {filename}")
            return None

        if df is not None and not df.empty:
            st.success(f"Successfully loaded and processed HPA data: {filename} ({len(df)} rows)")
            return df
        elif df is not None and df.empty:
             st.warning(f"Loaded HPA data file {filename}, but it resulted in an empty DataFrame after processing (check file content and parsing).")
             return None
        else:
             # This case should ideally not be reached if logic above is sound
             st.error(f"Failed to load data from {filename} for an unknown reason.")
             return None

    except ValueError as ve:
         # Catch errors specifically related to missing columns during read_csv
         st.error(f"Error loading HPA data from {file_path}. Potentially missing required columns: {ve}")
         return None
    except Exception as e:
        st.error(f"An unexpected error occurred loading HPA data from {file_path}: {e}")
        return None

@st.cache_data
def load_hpa_protein_data(file_path="data/normal_ihc_data.tsv"):
    """
    Loads HPA protein tissue consensus data and prints the column names.
    """
    try:
        # Load the data
        df = pd.read_csv(file_path, sep="\t")

        # Print the column names
        print("Columns in the DataFrame:")
        print(df.columns)

        return df
    except FileNotFoundError:
        print(f"File not found at {file_path}")
        return None
    except Exception as e:
        print(f"Error loading data: {e}")
        return None
def protein_expression(data: pd.DataFrame, selected_df_name: str):
    """
    Analyzes and visualizes HPA protein expression levels (IHC)
    for the selected protein list compared to all HPA proteins.
    """
    st.markdown("## Protein Expression Analysis (HPA - Immunohistochemistry)")
    st.markdown(
        """
        Comparing protein expression levels (Low, Medium, High) across tissues using
        Immunohistochemistry (IHC) data from the
        [Human Protein Atlas (HPA)](https://www.proteinatlas.org/about/download).
        The plot shows the fraction of proteins within each set ('All HPA Genes' vs. your selected list)
        that fall into each expression level category for a given tissue.
        This normalization accounts for the different total number of proteins in the two sets.
        """
    )

    # --- Load HPA Protein Data ---
    # Try the standard HPA filename first, fallback to the original name if needed
    hpa_df_full = load_hpa_protein_data(file_path="data/normal_tissue.tsv")
    if hpa_df_full is None:
        # Attempt fallback if primary load failed
        hpa_df_full = load_hpa_protein_data(file_path="data/normal_ihc_data.tsv") # Original name
        if hpa_df_full is None:
            st.warning("Cannot proceed with protein expression analysis: HPA protein data failed to load.")
            return # Stop if data loading failed

    # --- Verify required columns for analysis ---
    required_plot_cols = {'Gene', 'Tissue', 'Level'}
    if not required_plot_cols.issubset(hpa_df_full.columns):
        st.error(f"Loaded HPA Protein data is missing required columns for plotting: {required_plot_cols}. Found: {list(hpa_df_full.columns)}")
        return

    # --- Get Ensembl IDs from the selected dataset ---
    ensembl_ids = []
    id_source_col = None
    if 'Ensembl' in data.columns:
        ensembl_ids = data['Ensembl'].dropna().astype(str).unique()
        id_source_col = 'Ensembl'
        st.info("Using 'Ensembl' column for matching HPA protein data.")
    elif 'HPA' in data.columns:
        st.info("Using 'HPA' column for matching HPA protein data (extracting ENSG IDs).")
        id_source_col = 'HPA'
        def extract_ensg(hpa_string):
            if pd.isna(hpa_string): return None
            # Ensure input is string, split by ';', strip leading/trailing whitespace/chars
            parts = str(hpa_string).split(';')
            for part in parts:
                cleaned_part = part[1:]
                if cleaned_part.startswith('ENSG'):
                    # Basic check for valid ENSG format (ENSG + digits)
                    if len(cleaned_part) > 4 and cleaned_part[4:].isdigit():
                        return cleaned_part
            return None
        ensembl_ids = data['HPA'].apply(extract_ensg).dropna().astype(str).unique()
    else:
        st.error("Could not find a suitable column ('Ensembl' or 'HPA') in your dataset to link with HPA protein expression data.")
        return

    if len(ensembl_ids) == 0:
        st.warning(f"No valid Ensembl Gene IDs (starting with ENSG) found in the '{id_source_col}' column of the '{selected_df_name}' dataset after processing.")
        return
    else:
        st.write(f"Identified {len(ensembl_ids)} unique Ensembl IDs from column '{id_source_col}' in '{selected_df_name}'.")

    # --- Filter HPA Data for selected genes ---
    # Ensure HPA 'Gene' column is string for reliable matching
    hpa_df_full['Gene'] = hpa_df_full['Gene'].astype(str)
    hpa_selected_df = hpa_df_full[hpa_df_full['Gene'].isin(ensembl_ids)].copy()

    genes_found_hpa = hpa_selected_df['Gene'].nunique()
    if genes_found_hpa == 0:
         st.warning(f"No matching protein expression data found in HPA for the {len(ensembl_ids)} Ensembl IDs identified in '{selected_df_name}'. Analysis cannot proceed.")
         return
    else:
         st.success(f"Found HPA protein expression data (IHC) for {genes_found_hpa} out of {len(ensembl_ids)} genes from '{selected_df_name}'.")
         if genes_found_hpa < len(ensembl_ids):
              st.caption(f"{len(ensembl_ids) - genes_found_hpa} genes from your list were not found in the HPA protein dataset.")


    # --- Data Preparation for Plotting ---
    # Define the expression levels of interest
    valid_levels = ['Low', 'Medium', 'High']
    level_column = 'Level'
    tissue_column = 'Tissue'
    gene_column = 'Gene'

    # Filter both full and selected dataframes for valid levels and necessary columns
    df_all = hpa_df_full[[gene_column, tissue_column, level_column]].copy()
    df_all = df_all[df_all[level_column].isin(valid_levels)]
    df_all['Gene Type'] = 'All HPA Genes'


    df_selected = hpa_selected_df[[gene_column, tissue_column, level_column]].copy()
    df_selected = df_selected[df_selected[level_column].isin(valid_levels)]
    # No need to filter df_selected by level again, as it's a subset of df_all
    df_selected['Gene Type'] = selected_df_name # Use the dynamic name

    # Combine the dataframes
    combined_df = pd.concat([df_all, df_selected], ignore_index=True)

    # Get the count of unique proteins *present in the IHC data* for normalization
    # Important: Use the dataframes *after* filtering by valid levels if levels are sparse
    # However, HPA data usually covers most genes, so using nunique before level filtering is generally okay.
    # Let's calculate based on the data *before* filtering by level, to represent the total pool.
    unique_all_proteins = hpa_df_full[gene_column].nunique()
    unique_selected_proteins = hpa_selected_df[gene_column].nunique()

    if unique_selected_proteins == 0:
        st.error("Internal Error: Number of unique selected proteins with HPA data is zero. Cannot normalize.")
        return

    st.write(f"Normalization based on: {unique_all_proteins} unique genes in 'All HPA Genes' set (with IHC data), {unique_selected_proteins} unique genes in '{selected_df_name}' set (with IHC data).")


    # Group data to count occurrences of each level per tissue per gene type
    grouped_data = combined_df.groupby([tissue_column, level_column, 'Gene Type']).size().reset_index(name='Count')

    # Create a normalization factor map
    normalization_map = {'All HPA Genes': unique_all_proteins, selected_df_name: unique_selected_proteins}
    grouped_data['Total_Genes_In_Set'] = grouped_data['Gene Type'].map(normalization_map)

    # Calculate normalized counts (fraction of proteins in that set showing that level in that tissue)
    # Add a small epsilon to avoid division by zero, although Total_Genes_In_Set should be > 0 here
    grouped_data['Fraction'] = grouped_data['Count'] / (grouped_data['Total_Genes_In_Set'] + 1e-9)

    # --- Visualization ---
    st.subheader(f"Protein Expression Level Distribution: All HPA Genes vs. {selected_df_name}")

    # Determine number of tissues for facet wrapping
    num_tissues = grouped_data[tissue_column].nunique()
    facet_col_wrap = 3 if num_tissues > 6 else (2 if num_tissues > 3 else 1) # Adjust wrap based on number of tissues

    # Adjust height based on number of rows needed
    num_rows = -(-num_tissues // facet_col_wrap) # Ceiling division
    plot_height = max(600, num_rows * 300) # Minimum height 600px, scale with rows

    fig_protein = px.bar(
        grouped_data,
        x=level_column,
        y='Fraction', # Use the normalized value for the bar height
        color='Gene Type',
        barmode='group',
        facet_col=tissue_column,
        facet_col_wrap=facet_col_wrap,
        facet_row_spacing=0.04, # Adjust spacing as needed
        facet_col_spacing=0.03, # Adjust spacing as needed
        title=f'Protein Expression Levels by Tissue: All HPA Genes vs. {selected_df_name}',
        labels={
            'Fraction': 'Fraction of Proteins in Set', # Y-axis label
            level_column: 'Expression Level',
            tissue_column: 'Tissue',
            'Gene Type': 'Gene Set'
            },
        category_orders={
            level_column: valid_levels, # Ensure Low, Medium, High order
            'Gene Type': ['All HPA Genes', selected_df_name] # Control legend order
        },
        hover_data={
            'Gene Type': True,
            'Count': True,                   # Show raw count
            'Total_Genes_In_Set': ':,d',     # Show total unique proteins in the set
            'Fraction': ':.3f',             # Format fraction nicely
        },
        height=plot_height, # Dynamic height
        # width=1000 # Optional: set a fixed width if needed
    )

    # Improve layout
    fig_protein.update_yaxes(range=[0, 1]) # Set y-axis limit 0-1 for fractions
    fig_protein.update_xaxes(title_text='') # Remove redundant x-axis titles within facets
    # Add overall X axis title below the facets
    fig_protein.update_layout(
        margin=dict(l=40, r=40, t=80, b=80), # Adjust margins
        xaxis_title='Expression Level' # Add a general x-axis title at the bottom
    )
    # Rotate facet titles if they overlap (optional)
    # fig_protein.for_each_annotation(lambda a: a.update(textangle=0))
    # Make facet titles smaller if needed
    fig_protein.update_annotations(font_size=10)


    st.plotly_chart(fig_protein, use_container_width=True)

    # --- Display Data Table ---
    with st.expander("View Processed Protein Expression Data Table"):
         st.dataframe(grouped_data.round(3)) # Show rounded data

    st.markdown("---")
    st.markdown("Data Source: [Human Protein Atlas Normal Tissue IHC data](https://www.proteinatlas.org/about/download). Expression levels are categorized as Not detected, Low, Medium, High based on manual curation of IHC staining. This analysis focuses on Low, Medium, and High levels. Normalization calculates the fraction of proteins within each gene set exhibiting a specific expression level in a given tissue.")


def protein_expression2(data: pd.DataFrame, selected_df_name: str):
    """Analyzes HPA protein expression data for the selected protein list."""
    st.markdown("## Protein Expression Analysis (HPA)")
    st.markdown(
        """
        Comparing protein expression levels across tissues using consensus data from the
        [Human Protein Atlas (HPA)](https://www.proteinatlas.org/about/download).
        Expression is measured in normalized intensity.
        """
    )

    # --- Load HPA Data ---
    hpa_df = load_hpa_protein_data()
    if hpa_df is None:
        st.warning("Cannot proceed with protein expression analysis: HPA data failed to load.")
        return # Stop if data loading failed

    # --- Get Ensembl IDs from the selected dataset ---
    # Try to find an Ensembl ID column, fallback to HPA-like column extraction
    ensembl_ids = []
    if 'Ensembl' in data.columns:
        ensembl_ids = data['Ensembl'].dropna().unique()
        st.info("Using 'Ensembl' column for matching.")
    elif 'HPA' in data.columns:
        st.info("Using 'HPA' column for matching (extracting ENSG IDs).")
        def extract_ensg(hpa_string):
            if pd.isna(hpa_string): return None
            parts = str(hpa_string).split(';')
            for part in parts:
                cleaned_part = part[1:]
                if cleaned_part.startswith('ENSG'): return cleaned_part
            return None
        ensembl_ids = data['HPA'].apply(extract_ensg).dropna().unique()
    else:
        st.error("Could not find a suitable column ('Ensembl' or 'HPA') in the dataset to link with HPA protein expression data.")
        return

    if len(ensembl_ids) == 0:
        st.warning(f"No valid Ensembl Gene IDs (ENSG...) found in the '{selected_df_name}' dataset.")
        return

    # --- Filter HPA Data for selected genes ---
    hpa_moon_df = hpa_df[hpa_df['Gene'].isin(ensembl_ids)].copy() # Filter by Ensembl ID ('Gene' column in HPA tsv)

    genes_found = hpa_moon_df['Gene'].nunique()



def elm_analysis(data: pd.DataFrame, selected_df_name: str):
    """Analyzes ELM (Eukaryotic Linear Motif) instances in the selected protein list."""
    st.markdown("## ELM (Eukaryotic Linear Motif) Analysis")
    st.markdown(
        """
        Analyzing ELM instances in moonlighting proteins. This analysis uses data from the 
        [ELM database](http://elm.eu.org/) to identify linear motifs in the selected proteins
        and determine if they are enriched for specific ELM types.
        """
    )

    # Load ELM data
    elm_file = 'data/ELM_sapiens.tsv'

    try:
        elm_df = pd.read_csv(elm_file, sep='\t', comment='#')
        st.success(f"Successfully loaded ELM data with {elm_df.shape[0]} instances.")
    except Exception as e:
        st.error(f"Error loading ELM data: {e}")
        return

    # Get UniProt IDs from the selected dataset
    uniprot_ids = data['UniProtKB-AC'].dropna().unique().tolist()

    if len(uniprot_ids) == 0:
        st.warning(f"No UniProt IDs found in the '{selected_df_name}' dataset.")
        return

    # Filter ELM data for the selected proteins
    # The 'Primary_Acc' column contains the primary UniProt accession
    elm_moon_df = elm_df[elm_df['Primary_Acc'].isin(uniprot_ids)].copy()

    proteins_with_elm = elm_moon_df['Primary_Acc'].nunique()
    total_proteins = len(uniprot_ids)

    # Display summary statistics
    col1, col2 = st.columns(2)
    with col1:
        st.metric(
            label="Proteins with ELM instances",
            value=f"{proteins_with_elm}/{total_proteins}",
            delta=f"{proteins_with_elm/total_proteins:.1%}"
        )
    with col2:
        st.metric(
            label="Total ELM instances",
            value=elm_moon_df.shape[0]
        )

    # Display ELM types distribution
    st.subheader("ELM Types Distribution")
    elm_type_counts = elm_moon_df['ELMType'].value_counts().reset_index()
    elm_type_counts.columns = ['ELM Type', 'Count']

    # Create a bar chart
    fig = px.bar(
        elm_type_counts, 
        x='ELM Type', 
        y='Count',
        color='ELM Type',
        title="Distribution of ELM Types"
    )
    st.plotly_chart(fig)

    # Display ELM identifiers distribution (top 20)
    st.subheader("Top 20 ELM Identifiers")
    elm_id_counts = elm_moon_df['ELMIdentifier'].value_counts().reset_index()
    elm_id_counts.columns = ['ELM Identifier', 'Count']

    # Create a bar chart for top 20
    fig = px.bar(
        elm_id_counts.head(20), 
        x='ELM Identifier', 
        y='Count',
        color='Count',
        title="Top 20 ELM Identifiers"
    )
    st.plotly_chart(fig)

    # Create a dictionary mapping ELM types to UniProt IDs
    elm_type_dict = {}
    for elm_type in elm_moon_df['ELMType'].unique():
        elm_type_dict[elm_type] = elm_moon_df[elm_moon_df['ELMType'] == elm_type]['Primary_Acc'].unique().tolist()

    # Create a dictionary mapping ELM identifiers to UniProt IDs
    elm_id_dict = {}
    for elm_id in elm_moon_df['ELMIdentifier'].unique():
        elm_id_dict[elm_id] = elm_moon_df[elm_moon_df['ELMIdentifier'] == elm_id]['Primary_Acc'].unique().tolist()

    # Perform enrichment analysis for ELM types
    st.subheader("Enrichment Analysis for ELM Types")

    # Create a background dictionary of all human proteins with ELM instances
    all_human_elm_dict = {}
    for elm_type in elm_df['ELMType'].unique():
        all_human_elm_dict[elm_type] = elm_df[elm_df['ELMType'] == elm_type]['Primary_Acc'].unique().tolist()

    # Perform enrichment analysis
    results_df = enrichment_analysis(dict=all_human_elm_dict, query_uniprots=uniprot_ids, correction='fdr_bh',N=20000)

    # Display results
    st.dataframe(results_df)

    # Perform enrichment analysis for ELM identifiers
    st.subheader("Enrichment Analysis for ELM Identifiers")

    # Create a background dictionary of all human proteins with ELM identifiers
    all_human_elm_id_dict = {}
    for elm_id in elm_df['ELMIdentifier'].unique():
        all_human_elm_id_dict[elm_id] = elm_df[elm_df['ELMIdentifier'] == elm_id]['Primary_Acc'].unique().tolist()

    # Perform enrichment analysis
    results_id_df = enrichment_analysis(dict=all_human_elm_id_dict, query_uniprots=uniprot_ids, correction='fdr_bh',N=20000)

    # Display results (top 20 by p-value)
    st.dataframe(results_id_df.head(20))

    # Allow saving proteins with specific ELM types to session state
    st.subheader("Save Proteins with Specific ELM Types")

    elm_type_to_save = st.selectbox("Select ELM Type", elm_moon_df['ELMType'].unique())

    if st.button("Save to Session State"):
        proteins_with_elm_type = elm_moon_df[elm_moon_df['ELMType'] == elm_type_to_save]['Primary_Acc'].unique().tolist()
        session_key = f"{selected_df_name}_ELM_{elm_type_to_save}"

        if session_key not in st.session_state['gene_lists']:
            st.session_state['gene_lists'][session_key] = proteins_with_elm_type
            st.success(f"Saved {len(proteins_with_elm_type)} proteins with ELM type {elm_type_to_save} to session state with key: {session_key}")

def rna_expression(data: pd.DataFrame, selected_df_name: str):
    """Analyzes HPA RNA expression data (bulk or single-cell) for the selected protein list."""
    st.markdown("## RNA Expression Analysis (HPA)")
    st.markdown(
        """
        Comparing RNA expression levels using consensus data from the
        [Human Protein Atlas (HPA)](https://www.proteinatlas.org/about/download).
        Expression is measured in normalized Transcripts Per Million (nTPM).
        Select between bulk tissue data or single-cell type data.
        """
    )

    # --- User Selection: Bulk vs Single Cell ---
    analysis_type = st.radio(
        "Select HPA RNA Data Type:",
        ("Bulk Tissue", "Single Cell"),
        index=0,
        horizontal=True,
        key=f"rna_type_{selected_df_name}"
    )

    # --- Load Appropriate HPA Data ---
    hpa_data_file = None
    grouping_col = None
    hpa_df_all = None # Use a generic name for the loaded full dataset
    data_source_url = ""
    data_source_name = ""

    if analysis_type == "Bulk Tissue":
        hpa_data_file = 'rna_tissue_consensus.tsv'
        grouping_col = 'Tissue' # Grouping column for bulk data
        data_source_url = "https://www.proteinatlas.org/download/rna_consensus_tissue.tsv.zip"
        data_source_name = "HPA RNA Consensus Tissue data"

    elif analysis_type == "Single Cell":
        hpa_data_file = 'rna_single_cell_type.tsv'
        grouping_col = 'Cell type' # Grouping column for single-cell data
        data_source_url = "https://www.proteinatlas.org/download/rna_single_cell_type.tsv.zip"
        data_source_name = "HPA RNA Single Cell Type data"

    # Load data using the corrected function
    hpa_df_all = load_hpa_rna_data(hpa_data_file)

    if hpa_df_all is None:
        st.warning(f"Cannot proceed with {analysis_type} RNA expression analysis: HPA data ({hpa_data_file}) failed to load or is empty.")
        return # Stop if data loading failed

    # --- Get Ensembl IDs from the selected dataset ---
    # (This part remains the same, it extracts IDs from *your* input dataframe 'data')
    ensembl_ids = []
    id_column_used = None
    if 'Ensembl' in data.columns:
        ensembl_ids = data['Ensembl'].dropna().astype(str).unique() # Ensure string type
        id_column_used = 'Ensembl'


    if not id_column_used:
        st.error(f"Could not find a suitable column ('Ensembl' or 'HPA' containing ENSG IDs) in the '{selected_df_name}' dataset to link with HPA RNA expression data.")
        return

    # Filter out any potential non-ENSG IDs that might have slipped through
    ensembl_ids = [eid for eid in ensembl_ids if isinstance(eid, str) and eid.startswith('ENSG')]

    if len(ensembl_ids) == 0:
        st.warning(f"No valid Ensembl Gene IDs (starting with ENSG) found in the '{id_column_used}' column of the '{selected_df_name}' dataset after processing.")
        return

    st.write(f"Found {len(ensembl_ids)} unique valid Ensembl Gene IDs (ENSG...) in '{selected_df_name}'.")

    # --- Filter HPA Data for selected genes ---
    # Ensure 'Gene' column in HPA data is string for accurate matching
    hpa_df_all['Gene'] = hpa_df_all['Gene'].astype(str)
    hpa_moon_df = hpa_df_all[hpa_df_all['Gene'].isin(ensembl_ids)].copy()

    genes_found = hpa_moon_df['Gene'].nunique()
    genes_not_found = len(ensembl_ids) - genes_found

    if genes_found == 0:
         st.warning(f"No matching RNA expression data found in HPA {analysis_type} data for the {len(ensembl_ids)} valid genes identified in '{selected_df_name}'. Check if the Ensembl IDs exist in the selected HPA dataset.")
         # Optionally list IDs not found
         # st.write("IDs searched:", ensembl_ids)
         return
    else:
         st.success(f"Found HPA {analysis_type} RNA expression data for {genes_found} out of {len(ensembl_ids)} valid genes from '{selected_df_name}'.")
         if genes_not_found > 0:
             st.caption(f"({genes_not_found} genes from your list were not found in the selected HPA dataset)")


    # --- Calculate Aggregated Expression ---
    agg_func = st.radio("Select aggregation function:", ("Median", "Mean"), index=0, horizontal=True, key=f"rna_agg_{analysis_type}_{selected_df_name}")
    agg_method = 'median' if agg_func == "Median" else 'mean'
    agg_col_name_all = f"{agg_method}_nTPM_All_Genes"
    agg_col_name_moon = f"{agg_method}_nTPM_{selected_df_name}"

    # --- Perform Aggregation ---
    # Wrap in try-except in case grouping_col doesn't exist (shouldn't happen with loader fix)
    try:
        # Aggregate for all genes in the chosen HPA dataset
        grouped_all = hpa_df_all.groupby(grouping_col)["nTPM"].agg(agg_method).reset_index(name=agg_col_name_all)

        # Aggregate for the moonlighting genes subset within the chosen HPA dataset
        # Handle case where hpa_moon_df might be empty after filtering (although checked above)
        if not hpa_moon_df.empty:
            grouped_moon = hpa_moon_df.groupby(grouping_col)["nTPM"].agg(agg_method).reset_index(name=agg_col_name_moon)
        else:
             # Create an empty df with the right columns if no matching genes were found
             # This prevents merge errors later
             st.info("Creating empty aggregation for the selected gene set as no matches were found after filtering.")
             grouped_moon = pd.DataFrame(columns=[grouping_col, agg_col_name_moon])
             # Ensure grouping_col dtype matches grouped_all for merge compatibility
             grouped_moon[grouping_col] = grouped_moon[grouping_col].astype(grouped_all[grouping_col].dtype)


        # --- Merge and Visualize ---
        merged_expression = pd.merge(grouped_all, grouped_moon, on=grouping_col, how="left").fillna(0) # Use left merge to keep all tissues/cell types

        # Optional: Sort by expression of the selected set for better visualization
        # merged_expression = merged_expression.sort_values(by=agg_col_name_moon, ascending=False)

        st.subheader(f"{agg_func} nTPM Comparison by {grouping_col}")
        st.write(f"Comparing {agg_func.lower()} RNA expression levels across {grouping_col.lower()}s.")

        # Dynamically set labels based on grouping column
        plot_labels = {
            "value": f"{agg_func} nTPM",
            "variable": "Gene Set",
            grouping_col: grouping_col.replace("_", " ").title() # Nicer label formatting
        }

        fig_rna = px.bar(
            merged_expression,
            x=grouping_col, # Use the dynamic grouping column
            y=[agg_col_name_all, agg_col_name_moon],
            barmode='group',
            labels=plot_labels,
            title=f"{agg_func} nTPM by {grouping_col}: All HPA Genes vs. {selected_df_name}",
            height=600 if analysis_type == "Bulk Tissue" else 800 # Maybe more height for many cell types
        )
        fig_rna.update_layout(xaxis_tickangle=-90, legend_title_text='Gene Set')

        # Add specific layout adjustments if needed
        if analysis_type == "Single Cell" and not merged_expression.empty:
             try:
                 # Attempt to sort x-axis categories by the value of the user's selected gene set
                 sorted_categories = merged_expression.sort_values(by=agg_col_name_moon, ascending=False)[grouping_col].tolist()
                 fig_rna.update_xaxes(categoryorder='array', categoryarray=sorted_categories)
             except Exception as e:
                 st.caption(f"Could not sort x-axis automatically: {e}") # Fallback to default order

        st.plotly_chart(fig_rna, use_container_width=True)

        # --- Display Data Table ---
        with st.expander(f"View Aggregated {analysis_type} Expression Data Table"):
             # Show grouping column, and rounded expression values
             display_df = merged_expression[[grouping_col, agg_col_name_all, agg_col_name_moon]]
             st.dataframe(display_df.round(2), use_container_width=True)
        st.markdown("---")
        st.subheader(f"Single {grouping_col} nTPM Distribution")
        st.markdown(
            f"Select a specific **{grouping_col}** from the list below to see the distribution of nTPM values for individual genes within it."
        )

        # Get unique, sorted list of tissues/cell types present in the data
        # Use the merged_expression df's grouping col if available and populated, otherwise fallback to hpa_df_all
        if not merged_expression.empty and grouping_col in merged_expression.columns:
            grouping_options = sorted(merged_expression[grouping_col].unique())
        else:
            grouping_options = sorted(hpa_df_all[grouping_col].unique())

        # Add a placeholder option
        placeholder_option = f"--- Select a {grouping_col} ---"
        analysis_options = [placeholder_option] + grouping_options

        selected_entity = st.selectbox(
            f"Select {grouping_col} for distribution analysis:",
            options=analysis_options,
            index=None,  # Default to placeholder
            key=f"rna_entity_select_{analysis_type}_{selected_df_name}",  # Unique key
        )

        # --- Perform Analysis Only if a valid entity is selected ---
        if selected_entity != placeholder_option:
            st.write(f"Analyzing nTPM distribution within: **{selected_entity}**")

            # Filter the original dataframes for the selected entity
            hpa_all_filtered = hpa_df_all[
                hpa_df_all[grouping_col] == selected_entity
            ].copy()
            hpa_moon_filtered = hpa_moon_df[
                hpa_moon_df[grouping_col] == selected_entity
            ].copy()

            # Check if any data exists for the selected entity
            if hpa_all_filtered.empty:
                st.warning(
                    f"No HPA data found for the selected {grouping_col.lower()}: '{selected_entity}'. Cannot generate distribution plot."
                )
            else:
                # Prepare data for plotting histogram
                plot_data_all = hpa_all_filtered[
                    ["nTPM", "Gene name", "Gene"]
                ].copy()  # Include Gene ID
                plot_data_all["Gene Set"] = "All HPA Genes"

                if hpa_moon_filtered.empty:
                    st.info(
                        f"Note: None of the genes from '{selected_df_name}' had expression data recorded in '{selected_entity}'. Showing distribution for 'All HPA Genes' only."
                    )
                    plot_data_hist = plot_data_all
                    moon_count = 0
                else:
                    plot_data_moon = hpa_moon_filtered[
                        ["nTPM", "Gene name", "Gene"]
                    ].copy()  # Include Gene ID
                    plot_data_moon["Gene Set"] = f"{selected_df_name} Genes"
                    plot_data_hist = pd.concat(
                        [plot_data_all, plot_data_moon], ignore_index=True
                    )
                    moon_count = len(
                        hpa_moon_filtered
                    )  # Count of user's genes in this entity

                # Add option for log scale
                use_log_scale = st.checkbox(
                    "Use Log Scale for nTPM Axis",
                    value=True,
                    key=f"rna_dist_log_{analysis_type}_{selected_df_name}_{selected_entity}",  # Key specific to entity
                )

                # Handle log scale transformation
                plot_column = "nTPM_plot"
                if use_log_scale:
                    # Add 1 to handle nTPM = 0 before log transform. Use small epsilon if preferred.
                    plot_data_hist[plot_column] = np.log10(plot_data_hist["nTPM"] + 1)
                    x_axis_label = "log10(nTPM + 1)"
                else:
                    plot_data_hist[plot_column] = plot_data_hist["nTPM"]
                    x_axis_label = "nTPM"

                # --- Create Histogram ---
                total_all_count = len(hpa_all_filtered)
                plot_title = f"nTPM Distribution in {selected_entity}<br><sup>(All HPA Genes: {total_all_count}, {selected_df_name} Genes: {moon_count})</sup>"

                fig_dist = px.histogram(
                    plot_data_hist,
                    x=plot_column,
                    color="Gene Set",
                    marginal="rug",  # Show individual data points
                    barmode="overlay",  # Overlay histograms for comparison
                    opacity=0.65,  # Adjust transparency
                    labels={
                        plot_column: x_axis_label,
                        "Gene Set": "Gene Set",
                        "Gene name": "Gene Name",  # For hover
                        "Gene": "Ensembl ID",  # For hover
                    },
                    hover_data=[
                        "Gene name",
                        "Gene",
                        "nTPM",
                    ],  # Show useful info on hover
                    title=plot_title,
                    height=600,
                )
                fig_dist.update_layout(
                    xaxis_title=x_axis_label,
                    yaxis_title="Number of Genes",
                    legend_title_text="Gene Set",
                )

                st.plotly_chart(fig_dist, use_container_width=True)

                # --- Display Filtered Data Table (Single Entity) ---
                with st.expander(f"View Gene Expression Data for '{selected_entity}'"):
                    if hpa_moon_filtered.empty:
                        st.caption(
                            f"No genes from '{selected_df_name}' found in '{selected_entity}'."
                        )
                    else:
                        st.caption(
                            f"Showing nTPM values for {moon_count} genes from '{selected_df_name}' found in '{selected_entity}'."
                        )
                        st.dataframe(
                            hpa_moon_filtered[["Gene", "Gene name", "nTPM"]]
                            .sort_values("nTPM", ascending=False)
                            .round(2),
                            use_container_width=True,
                            # Optionally set height limit if needed
                            # height=300
                        )
    except KeyError as ke:
        st.error(f"Error during aggregation or merging: Missing expected column '{ke}'. This might indicate an issue with the loaded HPA data format for {analysis_type}.")
    except Exception as e:
         st.error(f"An unexpected error occurred during data aggregation or plotting: {e}")

    # --- Footer ---
    st.markdown("---")
    st.markdown(f"Data Source: [{data_source_name}]({data_source_url}). nTPM: Normalized Transcripts Per Million.")


# It's assumed that upsetplot is installed in the environment (e.g., via pip install upsetplot)


def cancer_driver(data: pd.DataFrame, selected_df_name: str):
    """
    Analyzes cancer driver-related data columns and visualizes the results.

    Parameters:
    - data (pd.DataFrame): Input DataFrame containing cancer-related columns
    - selected_df_name (str): Name of the selected dataset for display
    """
    st.markdown("## Cancer Driver Analysis")
    st.markdown(
        """
        Analysis of cancer driver annotations from multiple sources:
        - **Cancer Gene Census**: Catalogue of genes causally implicated in cancer
        - **Intogen**: Cancer driver genes identified through computational analysis
        - **OncoVar**: Cancer variants from ICGC and TCGA datasets
        """
    )

    # Import necessary for upsetplot
    import upsetplot
    from upsetplot import from_memberships

    cancer_sources = {
        "Cancer Gene Census": "CancerGeneCensus",
        "Intogen": "Intogen",
        "OncoVar ICGC": "OncoVar_ICGC",
        "OncoVar TCGA": "OncoVar_TCGA",
    }

    # Create a mutable copy for modifications if necessary, or ensure modifications are safe
    data_processed = data.copy()

    # Check if cancer source columns exist, otherwise, fill with False
    for col_name in cancer_sources.values():
        if col_name not in data_processed.columns:
            st.warning(
                f"Column '{col_name}' not found in the dataframe. It will be treated as all False for this analysis."
            )
            data_processed[col_name] = False
        else:
            # Ensure boolean type, fill NaNs with False
            data_processed[col_name] = (
                data_processed[col_name].fillna(False).astype(bool)
            )

    # Create metrics in columns
    cols = st.columns(len(cancer_sources))
    for i, (source_name, col_name) in enumerate(cancer_sources.items()):
        count = data_processed[col_name].sum()
        total = len(data_processed)
        percentage = (count / total) * 100 if total > 0 else 0
        cols[i].metric(
            label=source_name, value=f"{count}/{total}", delta=f"{percentage:.1f}%"
        )

    # Create upset plot for cancer drivers
    st.subheader("Overlap between Cancer Driver Sources")

    memberships = []
    for _, row in data_processed.iterrows():
        sources_present = []
        for source_name, col_name in cancer_sources.items():
            if row[col_name]:
                sources_present.append(source_name)
        if sources_present:  # Only add if protein is in at least one source
            memberships.append(
                tuple(sorted(sources_present))
            )  # Sort for consistent tuple representation

    if memberships:
        try:
            upset_data = from_memberships(memberships)
            if not upset_data.empty:
                fig_upset = plt.figure(figsize=(10, 6))
                upsetplot.plot(
                    upset_data, fig=fig_upset, subset_size="count"
                )  # Explicitly set subset_size
                st.pyplot(fig_upset)
                plt.clf()
            else:
                st.info(
                    "No data to display in Upset plot (empty after processing memberships)."
                )
        except Exception as e:
            st.error(f"Could not generate Upset plot: {e}")
            st.info(
                "Please ensure the 'upsetplot' library is installed and data is in the correct format. Explicitly setting subset_size='count' may help."
            )
    else:
        st.info(
            "No proteins found in any cancer driver sources to create an Upset plot."
        )

    with st.expander("About the Cancer Gene Census (CGC)"):
        st.markdown(
            """
            The Cancer Gene Census (CGC) is an ongoing effort to catalogue those genes which contain mutations that have been causally implicated in cancer and explain how dysfunction of these genes drives cancer. The content, the structure, and the curation process of the Cancer Gene Census was described and published in Nature Reviews Cancer.

            The census is not static, instead it is updated when new evidence comes to light. In particular we are grateful to Felix Mitelman and his colleagues in providing information on more genes involved in uncommon translocations in leukaemias and lymphomas. Currently, more than 1% of all human genes are implicated via mutation in cancer. Of these, approximately 90% contain somatic mutations in cancer, 20% bear germline mutations that predispose an individual to cancer and 10% show both somatic and germline mutations.

            **Census tiers**
            Genes in the Cancer Gene Census are divided into two groups, or tiers.
            """
        )

    # Analyze Role in Cancer distribution
    st.subheader("Role in Cancer Distribution")
    if "Role in Cancer" in data_processed.columns:
        role_counts = data_processed["Role in Cancer"].value_counts()
        if not role_counts.empty:
            fig_role = px.pie(
                values=role_counts.values,
                names=role_counts.index,
                title="Distribution of Cancer Roles",
                hole=0.3,
            )
            st.plotly_chart(fig_role)
        else:
            st.info("No data available for 'Role in Cancer'.")
    else:
        st.warning("Column 'Role in Cancer' not found.")

    with st.expander("Explanation of Cancer Gene Census Tiers"):
        st.markdown(
            """
            **Tier 1**
            To be classified into Tier 1, a gene must possess a documented activity relevant to cancer, along with evidence of mutations in cancer which change the activity of the gene product in a way that promotes oncogenic transformation. We also consider the existence of somatic mutation patterns across cancer samples gathered in COSMIC. For instance, tumour suppressor genes often show a broad range of inactivating mutations and dominant oncogenes usually demonstrate well defined hotspots of missense mutations. Genes involved in oncogenic fusions are included in Tier 1 when changes to their function caused by the fusion drives oncogenic transformation, or in cases when they provide regulatory elements to their partners (e.g. active promoter or dimerisation domain).

            **Tier 2**
            A new section of the Census, which consists of genes with strong indications of a role in cancer but with less extensive available evidence. These are generally more recent targets, where the body of evidence supporting their role is still emerging.
            """
        )

    # Analyze Tier distribution
    st.subheader("Cancer Gene Census Tier Distribution")
    if "Tier" in data_processed.columns:
        tier_counts = data_processed["Tier"].value_counts()
        if not tier_counts.empty:
            fig_tier = px.bar(
                x=tier_counts.index.astype(
                    str
                ),  # Ensure index is string for categorical axis
                y=tier_counts.values,
                title="Distribution of Cancer Gene Census Tiers",
                labels={"x": "Tier", "y": "Count"},
            )
            st.plotly_chart(fig_tier)
        else:
            st.info("No data available for 'Tier'.")
    else:
        st.warning("Column 'Tier' not found.")

    # Create detailed table of cancer drivers
    st.subheader("Cancer Driver Details")

    driver_cols_present = [
        col for col in cancer_sources.values() if col in data_processed.columns
    ]

    if not driver_cols_present:
        st.info("No cancer driver source columns found to identify cancer drivers.")
        return

    cancer_driver_mask = data_processed[driver_cols_present].any(axis=1)

    display_cols = []
    if "UniProtKB-AC" in data_processed.columns:
        display_cols.append("UniProtKB-AC")
    if "Gene symbol" in data_processed.columns:
        display_cols.append("Gene symbol")

    display_cols.extend(driver_cols_present)  # Add available driver columns

    if "Role in Cancer" in data_processed.columns:
        display_cols.append("Role in Cancer")
    if "Tier" in data_processed.columns:
        display_cols.append("Tier")

    # Ensure all display_cols actually exist in data_processed to avoid KeyErrors and ensure uniqueness
    final_display_cols = [
        col
        for col in list(dict.fromkeys(display_cols))
        if col in data_processed.columns
    ]

    if not final_display_cols:
        st.info("No relevant columns available to display for cancer drivers.")
        return

    cancer_drivers_df = data_processed.loc[cancer_driver_mask, final_display_cols]

    if not cancer_drivers_df.empty:
        st.dataframe(cancer_drivers_df)

        if st.button(
            "Save Cancer Drivers to Session State",
            key=f"save_cancer_drivers_{selected_df_name}",
        ):
            session_key = f"cancer_drivers_{selected_df_name}"
            if "gene_lists" not in st.session_state:
                st.session_state["gene_lists"] = {}

            if "UniProtKB-AC" in cancer_drivers_df.columns:
                st.session_state["gene_lists"][session_key] = (
                    cancer_drivers_df["UniProtKB-AC"].dropna().unique().tolist()
                )
                st.success(
                    f"Cancer drivers saved to session state with key: {session_key}"
                )
            else:
                st.error(
                    "Column 'UniProtKB-AC' not found in cancer drivers data. Cannot save."
                )
    else:
        st.info(
            "No cancer drivers found in the dataset based on the available source columns."
        )
# =======================
# Main Function
# =======================

def main():
    st.set_page_config(
        page_title="Moonlighting Proteins Analysis",
        layout="wide",
        initial_sidebar_state="expanded"
    )

    st.title("🌙 Moonlighting Proteins Analysis")

    # Load data
    df = load_data('data/moonhuman_mim.csv')  # Update with your actual file path
    if 'df' not in st.session_state:
        st.session_state['df'] = df
    if 'gene_lists' not in st.session_state:
        st.session_state['gene_lists'] = dict()
    st.write(
        "This interactive app allows for an exploratory data analysis of human moonlighting proteins. "
        "Each table provides functionalities for CSV export, specific value search, and column sorting. "
        "Below is a compendium of human moonlighting proteins obtained from the three main databases: "
        "[MoonProt](http://www.moonlightingproteins.org/), [MoonDB](http://moondb.hb.univ-amu.fr/), and [MultiTaskProtDB](http://wallace.uab.es/multitaskII). "
        "The MultiTaskProtDB dataset was obtained via email from one of the authors because the server has been down for months due to a cyber attack."
    )


    with st.expander("View Dataset Overview of Human Moonlighting Proteins (MPs)"):
        st.dataframe(df)

    # Venn Diagram
    st.subheader("Venn Diagram of Proteins Distribution Across Databases")
    st.write(
        "View the distribution of human moonlighting proteins across the three main databases: MoonDB, MoonProt, and MultiTaskProtDB."
    )
    plot_venn_diagram(df)

    # Filtering Section
    st.subheader("Filtering by membership in different databases")
    st.write("""
        Create a filtered dataset for further analysis by selecting the databases of interest and choosing whether you want the intersection or union of the results. 
        By default, there are already three predefined datasets:

        * **High-Consensus Dataset** (Intersection of MoonProt and MultiTaskProtDB - more restrictive) \t\t\t :green[humanMPs_highConsensus]
        * **Combined Literature Dataset** (Union of MoonProt and MultiTaskProtDB) \t\t\t :green[humanMPs_combinedLit]
        * **Comprehensive Dataset** (Union of all three databases - less restrictive) \t\t\t :green[humanMPs_comprehensive]
    """)

    filtered_data = filter_proteins(df)

    # Add genes to session state
    st.subheader("Select dataframe")
    st.write("Select the dataset to display analysis: ")

    df_names=sorted(get_lists_from_session())
    if'form_submitted' not in st.session_state:
        st.session_state['form_submitted'] = False
    if 'selected_df' not in st.session_state:
        st.session_state['selected_df'] = None

        # Define the form
    with st.form("data_selection_form"):
        st.selectbox(
            "Select a DataFrame for successive analysis:",
            df_names,
             index=df_names.index(st.session_state["selected_df"])
             if st.session_state["selected_df"]
             else 0,
            key='selected_df',
        )
        remove_unreviewed = st.checkbox("Remove unreviewed proteins")
        submitted = st.form_submit_button("Submit")
        if submitted:
            st.session_state['form_submitted'] = True
    selected_df_name = st.session_state['selected_df']
    if st.session_state['form_submitted']:
        selected_df_name = st.session_state['selected_df']
        st.write(f"Selected DataFrame: {selected_df_name}")
        selected_df_uniprots = st.session_state['gene_lists'][selected_df_name]
        data = df_from_uniprots(df, selected_df_uniprots)

        if remove_unreviewed:
            initial_count = data.shape[0]
            data = data[~data['Reviewed'].str.contains('unreviewed', case=False, na=False)]
            removed_count = initial_count - data.shape[0]
            st.success(f"Unreviewed proteins removed: {removed_count}. Total number of proteins: {data.shape[0]}.")
        else:
            st.info("Including unreviewed proteins.")

        # Display the DataFrame
        with st.expander("View Selected Data"):
            st.dataframe(data)

        # Display summary statistics
        total_proteins = data.shape[0]
        reviewed_count = data[data['Reviewed'] == 'reviewed'].shape[0]
        unreviewed_count = data['Reviewed'].str.contains('unreviewed', case=False, na=False).sum()

        st.write(
            f"**Total number of proteins:** {total_proteins}\n"
            f"**Number of reviewed proteins:** {reviewed_count}\n"
            f"**Number of non-reviewed proteins:** {unreviewed_count}"
        )



        tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9, tab10,tab11 = st.tabs(
            [
                "GO Terms",
                "Cancer Driver",
                "Protein Domains",
                "RNA Binding Proteins",
                "Enzymes",
                "Disorder",
                "Diseases",
                "RNA expression",
                "Protein Expression",
                "ELM Analysis",
                "Physicochemical Features"
            ]
        )

        with tab1:
            feature_analysis(data, selected_df_name)

        with tab2:
            cancer_driver(data, selected_df_name)

        with tab3:
            interpro_analysis(data, selected_df_name)

        with tab4:
            rna_binding_analysis(data, selected_df_name)

        with tab5:
            enzymes(data, selected_df_name)

        with tab6:
            disorder_analysis(data, selected_df_name)

        with tab7:
            disease(data, selected_df_name)

        with tab8:
            rna_expression(data, selected_df_name)

        with tab9:
            protein_expression(data, selected_df_name)

        with tab10:
            elm_analysis(data, selected_df_name)

        with tab11:
            sequence_analysis(data)



if __name__ == "__main__":
    main()
