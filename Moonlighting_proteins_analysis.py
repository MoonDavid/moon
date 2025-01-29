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
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import altair as alt
import json
import plotly.express as px
import plotly.graph_objects as go


# =======================
# Helper Functions
# =======================

@st.cache_data
def load_data(file_path: str) -> pd.DataFrame:
    """Load the dataset from a CSV file."""
    df = pd.read_csv(file_path)
    df = df.replace("", np.nan)
    return df


def enrichment_analysis(query_uniprots, dict, correction='fdr_bh'):
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

    N = len(all_uniprots)  # total background
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
        'pval_corrected': pvals_corrected
    })

    # Sort by corrected p-value ascending
    results_df.sort_values('pval_corrected', inplace=True)
    results_df.reset_index(drop=True, inplace=True)

    return results_df
def get_lists_from_session():
    """Retrieve keys from session_state that contain list objects."""
    lists = []
    for key, value in st.session_state.items():
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
        st.session_state[intersection_key] = list(intersection_set)
        #st.session_state["humanMPs(MoonProt AND MultiTaskProtDB)"] = data[data['UniProtKB-AC'].isin(intersection_set)].reset_index(drop=True)

        # 2. Union of MoonProt and MultiTaskProtDB
        union_two_set = sets_dict["MoonProt"].union(sets_dict["MultiTaskProtDB"])
        union_two_key = "humanMPs(MoonProt OR MultiTaskProtDB)"
        st.session_state[union_two_key] = list(union_two_set)
        #st.session_state["humanMPs(MoonProt OR MultiTaskProtDB)"] = data[data['UniProtKB-AC'].isin(union_two_set)].reset_index(drop=True)

        # 3. Union of all three databases (Least Restrictive)
        union_all_set = sets_dict["MoonDB"].union(sets_dict["MoonProt"]).union(sets_dict["MultiTaskProtDB"])
        union_all_key = "humanMPs_all"
        st.session_state[union_all_key] = list(union_all_set)
        #st.session_state["humanMPs(MoonDB OR MoonProt OR MultiTaskProtDB)"] = data[data['UniProtKB-AC'].isin(union_all_set)].reset_index(drop=True)

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
        st.session_state[session_key] = filtered_data['UniProtKB-AC'].tolist()
        #st.session_state[session_key2] = filtered_data

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
            st.write(f"Selected values: {selected_values}")
            selected_uniprot_ids = data[
                data[column].apply(lambda x: any(val in str(x).split("; ") for val in selected_values))
            ]["UniProtKB-AC"].dropna().unique().tolist()

            st.write(f"Associated UniProtKB-AC IDs: {selected_uniprot_ids}")
            session_key = f"{selected_df_name}-{column}:{selected_values}"
            st.session_state[session_key] = selected_uniprot_ids
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
            st.write(f"Associated UniProtKB-AC IDs: {selected_uniprot_ids}")
            session_key = f"{selected_df_name}-{column}:{selected_pairs}"
            st.session_state[session_key] = selected_uniprot_ids
            st.success(f"Saved in session state with key `{session_key}`.")



def sequence_analysis_deepseek(filtered_data: pd.DataFrame):
    """Analyze protein sequences with extended physicochemical properties using interactive Plotly plots and tabs."""

    if 'Sequence' in filtered_data.columns:
        # Initialize ProteinAnalysis objects
        st.write("## Protein Sequence Analysis")

        # Create tabs for different plots
        tab1, tab2, tab3, tab4, tab5 = st.tabs([
            "Sequence Length", "Amino Acid Composition", "Physicochemical Properties",
            "Molar Extinction Coefficient", "Secondary Structure"
        ])

        # 1. Sequence Length
        filtered_data['sequence_length'] = filtered_data['Sequence'].apply(lambda x: len(x) if pd.notnull(x) else 0)

        with tab1:
            st.write("### Protein Sequence Length:")
            st.write(filtered_data['sequence_length'].describe())

            fig = px.histogram(filtered_data, x='sequence_length', nbins=30, title="Distribution of Protein Sequence Lengths")
            fig.update_layout(xaxis_title="Sequence Length", yaxis_title="Frequency")
            st.plotly_chart(fig)

        # 2. Amino Acid Composition using ProtParam's get_amino_acids_percent()
        with tab2:
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

            fig = px.bar(mean_aa_composition, x=mean_aa_composition.index, y=mean_aa_composition.values,
                         labels={'x': 'Amino Acid', 'y': 'Percentage (%)'},
                         title="Average Amino Acid Composition of Protein Sequences")
            st.plotly_chart(fig)

            st.write("### Amino Acid Frequency:")
            fig = px.imshow(aa_composition_df.T, labels=dict(x="Proteins", y="Amino Acids", color="Percentage"),
                            title="Heatmap of Amino Acid Composition")
            st.plotly_chart(fig)

        # 3. Physicochemical Properties
        with tab3:
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
            def safe_histplot(data, xlabel, ylabel, title):
                data = data.dropna()
                if len(data) > 0:
                    fig = px.histogram(data, x=data, nbins=30, title=title)
                    fig.update_layout(xaxis_title=xlabel, yaxis_title=ylabel)
                    st.plotly_chart(fig)
                else:
                    st.write(f"No data available for {title}.")

            # 3.1 Molecular Weight Statistics
            st.write("#### Molecular Weight Statistics:")
            st.write(filtered_data['Molecular_Weight'].describe())
            safe_histplot(filtered_data['Molecular_Weight'], "Molecular Weight (Da)", "Frequency", "Distribution of Protein Molecular Weight")

            # 3.2 Isoelectric Point Statistics
            st.write("#### Isoelectric Point (pI) Statistics:")
            st.write(filtered_data['Isoelectric_Point'].describe())
            safe_histplot(filtered_data['Isoelectric_Point'], "Isoelectric Point (pI)", "Frequency", "Distribution of Protein Isoelectric Point")

            # 3.3 Aromaticity Statistics
            st.write("#### Aromaticity Statistics:")
            st.write(filtered_data['Aromaticity'].describe())
            safe_histplot(filtered_data['Aromaticity'], "Aromaticity", "Frequency", "Distribution of Protein Aromaticity")

            # 3.4 Instability Index Statistics
            st.write("#### Instability Index Statistics:")
            st.write(filtered_data['Instability_Index'].describe())
            safe_histplot(filtered_data['Instability_Index'], "Instability Index", "Frequency", "Distribution of Protein Instability Index")

            # Highlight unstable proteins (Instability Index > 40)
            unstable = filtered_data[filtered_data['Instability_Index'] > 40]
            st.write(f"**Number of Unstable Proteins (Instability Index > 40):** {unstable.shape[0]}")

            # 3.6 Gravy (Grand Average of Hydropathy) Statistics
            st.write("#### Gravy (Grand Average of Hydropathy) Statistics:")
            st.write(filtered_data['Gravy'].describe())
            safe_histplot(filtered_data['Gravy'], "Gravy (Hydropathy)", "Frequency", "Distribution of Protein Gravy (Hydropathy)")

            # 3.8 Charge at pH 7.0 Statistics
            st.write("#### Charge at pH 7.0 Statistics:")
            st.write(filtered_data['Charge_pH7.0'].describe())
            safe_histplot(filtered_data['Charge_pH7.0'], "Charge at pH 7.0", "Frequency", "Distribution of Protein Charge at pH 7.0")

        # 3.7 Molar Extinction Coefficient Statistics
        with tab4:
            st.write("#### Molar Extinction Coefficient (ε) Statistics:")
            st.write(filtered_data[['Molar_Extinction_Reduced', 'Molar_Extinction_Oxidized']].describe())

            epsilon_reduced = filtered_data['Molar_Extinction_Reduced'].dropna()
            epsilon_oxidized = filtered_data['Molar_Extinction_Oxidized'].dropna()

            if len(epsilon_reduced) > 0 or len(epsilon_oxidized) > 0:
                fig = px.histogram(epsilon_reduced, nbins=30, title="Distribution of Molar Extinction Coefficient (Reduced Cysteines)")
                fig.update_layout(xaxis_title="Molar Extinction Coefficient (ε)", yaxis_title="Frequency")
                st.plotly_chart(fig)

                fig = px.histogram(epsilon_oxidized, nbins=30, title="Distribution of Molar Extinction Coefficient (Disulfide Bridges)")
                fig.update_layout(xaxis_title="Molar Extinction Coefficient (ε)", yaxis_title="Frequency")
                st.plotly_chart(fig)
            else:
                st.write("No data available to plot Molar Extinction Coefficient.")

        # 4. Secondary Structure Analysis
        with tab5:
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
                fig = px.box(ss_df, title="Distribution of Secondary Structure Fractions")
                fig.update_layout(yaxis_title="Fraction", xaxis_title="Secondary Structure Type")
                st.plotly_chart(fig)
            else:
                st.write("No data available to plot Secondary Structure Fractions.")

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

def sequence_analysis1(filtered_data: pd.DataFrame):
    """Analyze protein sequences with extended physicochemical properties."""

    if 'Sequence' in filtered_data.columns:
        # Initialize ProteinAnalysis objects
        st.write("## Protein Sequence Analysis")

        # 1. Sequence Length
        filtered_data['sequence_length'] = filtered_data['Sequence'].apply(lambda x: len(x) if pd.notnull(x) else 0)

        st.write("### Protein Sequence Length:")
        st.write(filtered_data['sequence_length'].describe())

        plt.figure(figsize=(10, 6))
        sns.histplot(filtered_data['sequence_length'], kde=True, bins=30, color='salmon')
        plt.xlabel("Sequence Length")
        plt.ylabel("Frequency")
        plt.title("Distribution of Protein Sequence Lengths")
        st.pyplot(plt)

        # 2. Amino Acid Composition using ProtParam's get_amino_acids_percent()
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

        plt.figure(figsize=(14, 7))
        sns.barplot(x=mean_aa_composition.index, y=mean_aa_composition.values, palette="magma")
        plt.xticks(rotation=90)
        plt.xlabel("Amino Acid")
        plt.ylabel("Percentage (%)")
        plt.title("Average Amino Acid Composition of Protein Sequences")
        st.pyplot(plt)

        st.write("### Amino Acid Frequency:")
        plt.figure(figsize=(14, 7))
        sns.heatmap(aa_composition_df.T, cmap="YlGnBu", cbar=True)
        plt.xlabel("Proteins")
        plt.ylabel("Amino Acids")
        plt.title("Heatmap of Amino Acid Composition")
        st.pyplot(plt)

        # 3. Physicochemical Properties
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
        def safe_histplot(data, xlabel, ylabel, title, color):
            data = data.dropna()
            if len(data) > 1:
                sns.histplot(data, kde=True, bins=30, color=color)
                plt.xlabel(xlabel)
                plt.ylabel(ylabel)
                plt.title(title)
                st.pyplot(plt)
            elif len(data) == 1:
                st.write(f"Only one data point available for {title}. Skipping KDE.")
                sns.histplot(data, kde=False, bins=30, color=color)
                plt.xlabel(xlabel)
                plt.ylabel(ylabel)
                plt.title(title)
                st.pyplot(plt)
            else:
                st.write(f"No data available for {title}.")

        # 3.1 Molecular Weight Statistics
        st.write("#### Molecular Weight Statistics:")
        st.write(filtered_data['Molecular_Weight'].describe())

        plt.figure(figsize=(10, 6))
        safe_histplot(
            filtered_data['Molecular_Weight'],
            xlabel="Molecular Weight (Da)",
            ylabel="Frequency",
            title="Distribution of Protein Molecular Weight",
            color='teal'
        )

        # 3.2 Isoelectric Point Statistics
        st.write("#### Isoelectric Point (pI) Statistics:")
        st.write(filtered_data['Isoelectric_Point'].describe())

        plt.figure(figsize=(10, 6))
        safe_histplot(
            filtered_data['Isoelectric_Point'],
            xlabel="Isoelectric Point (pI)",
            ylabel="Frequency",
            title="Distribution of Protein Isoelectric Point",
            color='orange'
        )

        # 3.3 Aromaticity Statistics
        st.write("#### Aromaticity Statistics:")
        st.write(filtered_data['Aromaticity'].describe())

        plt.figure(figsize=(10, 6))
        safe_histplot(
            filtered_data['Aromaticity'],
            xlabel="Aromaticity",
            ylabel="Frequency",
            title="Distribution of Protein Aromaticity",
            color='purple'
        )

        # 3.4 Instability Index Statistics
        st.write("#### Instability Index Statistics:")
        st.write(filtered_data['Instability_Index'].describe())

        plt.figure(figsize=(10, 6))
        safe_histplot(
            filtered_data['Instability_Index'],
            xlabel="Instability Index",
            ylabel="Frequency",
            title="Distribution of Protein Instability Index",
            color='darkred'
        )

        # Highlight unstable proteins (Instability Index > 40)
        unstable = filtered_data[filtered_data['Instability_Index'] > 40]
        st.write(f"**Number of Unstable Proteins (Instability Index > 40):** {unstable.shape[0]}")




        # 3.6 Gravy (Grand Average of Hydropathy) Statistics
        st.write("#### Gravy (Grand Average of Hydropathy) Statistics:")
        st.write(filtered_data['Gravy'].describe())

        plt.figure(figsize=(10, 6))
        safe_histplot(
            filtered_data['Gravy'],
            xlabel="Gravy (Hydropathy)",
            ylabel="Frequency",
            title="Distribution of Protein Gravy (Hydropathy)",
            color='gold'
        )

        # 3.7 Molar Extinction Coefficient Statistics
        st.write("#### Molar Extinction Coefficient (ε) Statistics:")
        st.write(filtered_data[['Molar_Extinction_Reduced', 'Molar_Extinction_Oxidized']].describe())

        # Plot Molar Extinction Coefficients
        epsilon_reduced = filtered_data['Molar_Extinction_Reduced'].dropna()
        epsilon_oxidized = filtered_data['Molar_Extinction_Oxidized'].dropna()

        plt.figure(figsize=(10, 6))
        if len(epsilon_reduced) > 0:
            sns.histplot(epsilon_reduced, kde=True, bins=30, color='navy', label='Reduced Cysteines')
        else:
            st.write("No data available for Molar Extinction Coefficient (Reduced Cysteines).")

        if len(epsilon_oxidized) > 0:
            sns.histplot(epsilon_oxidized, kde=True, bins=30, color='magenta', label='Disulfide Bridges')
        else:
            st.write("No data available for Molar Extinction Coefficient (Disulfide Bridges).")

        if len(epsilon_reduced) > 0 or len(epsilon_oxidized) > 0:
            plt.xlabel("Molar Extinction Coefficient (ε)")
            plt.ylabel("Frequency")
            plt.title("Distribution of Molar Extinction Coefficient")
            plt.legend()
            st.pyplot(plt)
        else:
            st.write("No data available to plot Molar Extinction Coefficient.")

        # 3.8 Charge at pH 7.0 Statistics
        st.write("#### Charge at pH 7.0 Statistics:")
        st.write(filtered_data['Charge_pH7.0'].describe())

        plt.figure(figsize=(10, 6))
        safe_histplot(
            filtered_data['Charge_pH7.0'],
            xlabel="Charge at pH 7.0",
            ylabel="Frequency",
            title="Distribution of Protein Charge at pH 7.0",
            color='brown'
        )

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
        plt.figure(figsize=(10, 6))
        ss_df = filtered_data[['Fraction_Alpha_Helix', 'Fraction_Turn', 'Fraction_Beta_Sheet']].dropna()
        if not ss_df.empty:
            sns.boxplot(data=ss_df, palette="viridis")
            plt.ylabel("Fraction")
            plt.xlabel("Secondary Structure Type")
            plt.title("Distribution of Secondary Structure Fractions")
            st.pyplot(plt)
        else:
            st.write("No data available to plot Secondary Structure Fractions.")


def safe_literal_eval(val):
    try:
        return ast.literal_eval(val)
    except (ValueError, SyntaxError):
        return val
def interpro_analysis(data: pd.DataFrame, selected_df_name):
    """Analyze InterPro data."""
    st.write(
        "InterPro provides functional analysis of proteins by classifying them into families and predicting domains and important sites. To classify proteins in this way, InterPro uses predictive models, known as signatures, provided by several different databases (referred to as member databases) that make up the InterPro consortium:"
        " CATH, CDD, HAMAP, MobiDB Lite, Panther, Pfam, PIRSF, PRINTS, Prosite, SFLD, SMART, SUPERFAMILY AND NCBIfam. The [InterPro Consortium](https://interpro-documentation.readthedocs.io/en/latest/databases.html) section gives further information about the individual databases.")

    # Function to safely evaluate strings to Python literals

    abbr= {
        "TI": "ncbifam",
        "PT": "panther",
        "PS": "profile",
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
        "IP": "InterPro"
    }
    # Convert the 'interpro_data' column from string to actual tuples
    data.loc[:, 'interpro_data'] = data['interpro_data'].apply(safe_literal_eval)

    # Explode the 'interpro_data' column
    df_exploded = data.explode('interpro_data')
    # Count the occurrences of each value in the 'interpro_data' column
    count = df_exploded['interpro_data'].value_counts()
    count=count.reset_index(name='count')
    count['weblink']=count['interpro_data'].apply(lambda x: f"https://www.ebi.ac.uk/interpro/entry/{abbr[x[0][:2]]}/{x[0]}" if isinstance(x,tuple) else "")
    #st.write(f"Selected values: {selected_values}")
    st.write("Most common InterPro domains:")
    event = st.dataframe(
        count,
        column_config={
            "weblink": st.column_config.LinkColumn(
                "weblink",
                help="Link to InterPro entry"
            )
        },
        use_container_width=True,
        hide_index=True,
        on_select="rerun",
        selection_mode="multi-row",
    )
    if event.selection and st.button("Save selection"):
        selected_rows = event.selection.rows  # Get selected row indices
        selected_values = count.iloc[selected_rows]['interpro_data'].values.tolist()  # Filter dataframe
        uniprot_ids = df_exploded[df_exploded['interpro_data'].isin(selected_values)]['UniProtKB-AC'].values.tolist()
        session_key = f"{selected_df_name}-{selected_values}"
        st.write(f"Selected values: {selected_values}")
        st.write(f"Associated UniProtKB-AC IDs: {uniprot_ids}")
        if session_key not in st.session_state:
            st.session_state[session_key] = uniprot_ids
            st.success(f"Associated UniProtKB-AC IDs saved to session state with key: {session_key}")
    st.write("### Enrichment of InterPro domains:")
    #read dict from domain_dict.json
    domain_dict = json.load(open("domain_dict.json", "r"))
    #delete all keys that starts with DP
    domain_dict = {k: v for k, v in domain_dict.items() if not k.startswith("DP")}
    results_df=enrichment_analysis(dict=domain_dict, query_uniprots=df_to_uniprots(data), correction='fdr_bh')
    st.subheader("Enrichment Results")

    # Add a column for negative log10 p-value (for plotting)
    #results_df["-log10(pval_corrected)"] = -results_df["pval_corrected"].apply(lambda x: 1e-300 if x <= 0 else x).apply(np.log10)
    results_df['weblink'] = results_df['term'].apply(
        lambda x: f"https://www.ebi.ac.uk/interpro/entry/{abbr[x[:2]]}/{x}" if x[1].isalpha() else f"https://www.ebi.ac.uk/interpro/entry/cathgene3d/G3DSA:{x}")
    st.dataframe(results_df,column_config={
            "weblink": st.column_config.LinkColumn(
                "weblink",
                help="Link to InterPro entry"
            )},hide_index=True
        )


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
            st.session_state[sessionkey] = df_to_uniprots(datarbp)
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

def disorder_analysis(data: pd.DataFrame, selected_df_name: str):
    """
    Analyzes disorder-related columns in the given DataFrame and visualizes the data.

    Parameters:
    - data (pd.DataFrame): The input DataFrame containing disorder information.
    - selected_df_name (str): The name of the selected DataFrame for display purposes.
    """
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
    st.subheader("Missing Values Summary")
    nan_counts = dataidp.isna().sum()
    nan_df = nan_counts.reset_index()
    nan_df.columns = ['Column', 'Number of NaNs']
    st.table(nan_df)

    # Visualize the distribution of disorder columns
    st.subheader("Disorder regions percentage")

    # Identify numerical disorder columns for plotting
    numerical_cols = [
        'disprot_disorder',
        'curated_disorder',
        'alphadb_disorder',
        'mobidblite_disorder'
    ]
    idp_dict = json.load(open("IDP_ontology.json", "r"))
    # delete all keys that starts with DP
    idp_dict = {k: v for k, v in idp_dict.items() if not k.startswith("DP")}
    results_df = enrichment_analysis(dict=idp_dict, query_uniprots=df_to_uniprots(data), correction='fdr_bh')
    st.subheader("Enrichment Results")
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


# =======================
# Main Function
# =======================

def main():
    st.set_page_config(
        page_title="Moonlighting Proteins Analysis",
        layout="wide",
        initial_sidebar_state="expanded"
    )

    st.title("Moonlighting Proteins Analysis")

    # Load data
    df = load_data('moonhuman.csv')  # Update with your actual file path
    if 'df' not in st.session_state:
        st.session_state['df'] = df
    st.write(
        "Questa app interattiva permette un'analisi dati esplorativa delle proteine moonlighting umane. "
        "Ogni tabella offre funzionalità di esportazione in formato CSV, ricerca di valori specifici e ordinamento delle colonne. "
        "Qui di seguito è mostrato un compendio di proteine moonlighting umane ottenute dai 3 principali database: "
        "[MoonProt](http://www.moonlightingproteins.org/), [MoonDB](http://moondb.hb.univ-amu.fr/) e [MultiTaskProtDB](http://wallace.uab.es/multitaskII). "
        "Il dataset di MultiTaskProtDB è stato ottenuto tramite mail da uno degli autori perchè il server è down da mesi per attacco informatico.")
    # Display dataset overview
    st.subheader("Dataset Overview of Human Moonlighting Proteins (MPs)")
    st.dataframe(df)

    # Venn Diagram
    st.subheader("Venn Diagram of Proteins Distribution Across Databases")
    st.write(
        "Visualizza la distribuzione delle proteine moonlighting umane nei 3 principali database: MoonDB, MoonProt e MultiTaskProtDB")
    plot_venn_diagram(df)

    # Filtering Section
    st.subheader("Filtering by membership in different databases")
    st.write("""
        Crea un dataset filtrato per le analisi successive selezionando i database di interesse e se vuoi l'intersezione o l'unione dei risultati.
        Di default ce ne sono gia' 3:

        * Intersezione di MoonProt e MultiTaskProtDB (più restrittivo) \t\t\t :green[humanMPs(MoonProt AND MultiTaskProtDB)]
        * Unione di MoonProt e MultiTaskProtDB \t\t\t :green[humanMPs(MoonProt OR MultiTaskProtDB)]
        * Unione di tutti e 3 i database (meno restrittivo) \t\t\t :green[humanMPs_all]
    """)

    filtered_data = filter_proteins(df)

    # Add genes to session state
    st.subheader("Select dataframe")
    st.write("Seleziona il dataset filtrato da usare per le analisi successive:")

    df_names=sorted(get_lists_from_session())

    def on_change():
        st.session_state['index'] = df_names.index(st.session_state.selected_df)

    if 'index' not in st.session_state:
        st.session_state['index'] = None

        # Define the form
    with st.form("data_selection_form"):
        selected_df_name = st.selectbox(
            "Select a DataFrame for successive analysis:",
            df_names,
            key='selected_df',
            index=st.session_state['index']
        )

        remove_unreviewed = st.checkbox("Remove unreviewed proteins")

        # Submit button
        submitted = st.form_submit_button("Submit")

        if submitted:
            # Update session state with the selected DataFrame name
            st.session_state['index'] = df_names.index(selected_df_name) if selected_df_name in df_names else 0



    # Use the stored selections if available
    if selected_df_name in st.session_state:
        selected_df_uniprots = st.session_state[selected_df_name]
        data = df_from_uniprots(df, selected_df_uniprots)

        if remove_unreviewed:
            initial_count = data.shape[0]
            data = data[~data['Reviewed'].str.contains('unreviewed', case=False, na=False)]
            removed_count = initial_count - data.shape[0]
            st.success(f"Unreviewed proteins removed: {removed_count}. Total number of proteins: {data.shape[0]}.")
        else:
            st.info("Including unreviewed proteins.")

        # Display the DataFrame
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

        # Subsequent Analyses Sections
        st.subheader("GO terms")
        with st.expander("See analysis"):
            feature_analysis(data, selected_df_name)

        st.subheader("Protein physicochemical features ")
        # Automatically run analysis without a button
        with st.expander("See analysis Gemini"):
            sequence_analysis(data)

        st.subheader("RNA binding proteins")
        with st.expander("See analysis"):
            rna_binding_analysis(data, selected_df_name)

        st.subheader("Intrinsically disordered regions and proteins")
        with st.expander("See analysis"):
            disorder_analysis(data, selected_df_name)


if __name__ == "__main__":
    main()
