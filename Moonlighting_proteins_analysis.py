import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import itertools
from collections import Counter
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from gprofiler import GProfiler
import string
import requests
from urllib.parse import urlencode
import ast


# =======================
# Helper Functions
# =======================

@st.cache_data
def load_data(file_path: str) -> pd.DataFrame:
    """Load the dataset from a CSV file."""
    df = pd.read_csv(file_path)
    df = df.replace("", np.nan)
    return df


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

    plt.figure(figsize=(8, 8))
    venn = venn3([moondb_set, moonprot_set, multitask_set],
                 set_labels=('MoonDB', 'MoonProt', 'MultiTaskProtDB'))
    plt.title("Venn diagram showing the distribution of proteins across the 3 main databases")
    st.pyplot(plt)


def filter_proteins(data: pd.DataFrame) -> pd.DataFrame:
    """
    Filter proteins based on selected databases and operation using a form layout.
    Stores results in session state with a descriptive key indicating the operation.
    """
    st.subheader("Filtering by membership in different databases")

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
                key='filter_sets'
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
        else:  # OR operation
            filtered_set = set.union(*selected_sets_data)

        # Filter the data
        filtered_data = data[data['UniProtKB-AC'].isin(filtered_set)]

        # Create a descriptive key for session state
        # Convert database names to title case and join them
        db_names = [name.title() for name in selected_sets]
        operation_str = " AND " if operation == "Intersection" else " OR "
        session_key = f"humanMPs({operation_str.join(db_names)})_uniprotids"
        session_key2=f"humanMPs({operation_str.join(db_names)})"

        # Store the filtered UniProtKB-AC list in session state
        st.session_state[session_key] = filtered_data['UniProtKB-AC'].tolist()
        st.session_state[session_key2] = filtered_data

        # Display results
        st.write(f"Number of filtered proteins: {filtered_data.shape[0]}")
        st.write(f"Results stored in session state with key: {session_key}")
        st.dataframe(filtered_data)
        st.session_state.run_analysis=True

        return filtered_data

    elif submitted and not selected_sets:
        st.warning("Please select at least one database for filtering.")
        return data

    return data

def get_dataframes_from_session():
    dataframes = {}
    for key, value in st.session_state.items():
        if isinstance(value, pd.DataFrame):
            dataframes[key] = value
    return dataframes

def feature_analysis(filtered_data: pd.DataFrame):
    """Analyze features related to Gene Ontology."""
    st.subheader("Feature Analysis")

    # Filter columns related to Gene Ontology
    go_columns = [col for col in filtered_data.columns if "Gene Ontology" in col]
    column = st.selectbox("Select a feature", go_columns, key='feature_select')

    if column:
        if pd.api.types.is_numeric_dtype(filtered_data[column]):
            display_numeric_feature(filtered_data, column)
        else:
            display_categorical_feature(filtered_data, column)


def display_numeric_feature(data: pd.DataFrame, column: str):
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


def display_categorical_feature(data: pd.DataFrame, column: str):
    """Display analysis for categorical features with UniProtKB-AC tracking."""

    # Display basic information about the column
    st.write(f"The selected column `{column}` is **categorical** or treated as such.")
    unique_values = data[column].dropna().unique()
    st.write(f"### Example rows in `{column}`: ")
    st.write(unique_values[:10] if len(unique_values) > 10 else unique_values)

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
            session_key = f"{column} - selected_values"
            st.session_state[session_key] = selected_uniprot_ids
            st.success(f"Saved in session state with key `{session_key}`.")


    else:
        st.warning("No data available for visualization.")

    # Analyze Gene Ontology pairs if applicable
    if "Gene Ontology" in column:
        st.write(f"### Count of most common pairs in `{column}`:")
        go_terms = data[column].dropna().str.split('; ')
        go_pairs = go_terms.apply(lambda x: list(itertools.combinations(sorted(x), 2)))
        go_pairs = go_pairs.explode().dropna()

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

            # Save associated UniProtKB-AC values
            selected_uniprot_ids = []
            for pair in selected_pairs:
                element1, element2 = pair
                matched_rows = data[column].dropna().apply(
                    lambda x: element1 in x and element2 in x
                )
                selected_uniprot_ids.extend(data.loc[matched_rows, "UniProtKB-AC"].dropna().unique())

            selected_uniprot_ids = list(set(selected_uniprot_ids))  # Remove duplicates
            st.write(f"Associated UniProtKB-AC IDs: {selected_uniprot_ids}")
            session_key = f"{column} - {selected_pairs}"
            st.session_state[session_key] = selected_uniprot_ids
            st.success(f"Saved in session state with key `{session_key}`.")


def sequence_analysis(filtered_data: pd.DataFrame):
    """Analyze protein sequences."""
    st.header("Protein Sequence Analysis")

    if 'Sequence' in filtered_data.columns:
        st.subheader("Protein Sequence Exploration")

        if not filtered_data['Sequence'].empty:
            st.write("### Example Protein Sequence:")
            st.text(filtered_data['Sequence'].iloc[0])
        else:
            st.write("No sequence available in the filtered dataset.")

        # Calculate sequence lengths
        filtered_data['sequence_length'] = filtered_data['Sequence'].apply(lambda x: len(x) if pd.notnull(x) else 0)

        st.write("### Protein Sequence Length:")
        st.write(filtered_data['sequence_length'].describe())

        plt.figure(figsize=(10, 6))
        sns.histplot(filtered_data['sequence_length'], kde=True, bins=30, color='salmon')
        plt.xlabel("Sequence Length")
        plt.ylabel("Frequency")
        plt.title("Distribution of Protein Sequence Lengths")
        st.pyplot(plt)

        # Amino Acid Composition
        st.write("### Amino Acid Composition:")

        def get_aa_composition(seq: str) -> dict:
            """Calculate amino acid composition percentage."""
            aa_counts = Counter(seq)
            total = sum(aa_counts.values())
            return {aa: (count / total * 100) for aa, count in aa_counts.items()}

        aa_composition = filtered_data['Sequence'].dropna().apply(get_aa_composition)
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

        # Molecular Weight and Isoelectric Point
        st.write("### Physicochemical Properties of Protein Sequences:")

        def calculate_properties(seq: str) -> pd.Series:
            """Calculate molecular weight and isoelectric point of a protein sequence."""
            if pd.isnull(seq) or not seq:
                return pd.Series({"Molecular_Weight": None, "Isoelectric_Point": None})
            try:
                prot_param = ProteinAnalysis(seq)
                mw = prot_param.molecular_weight()
                pI = prot_param.isoelectric_point()
                return pd.Series({"Molecular_Weight": mw, "Isoelectric_Point": pI})
            except:
                return pd.Series({"Molecular_Weight": None, "Isoelectric_Point": None})

        properties = filtered_data['Sequence'].apply(calculate_properties)
        filtered_data = pd.concat([filtered_data, properties], axis=1)

        st.write("#### Molecular Weight Statistics:")
        st.write(filtered_data['Molecular_Weight'].describe())

        st.write("#### Isoelectric Point (pI) Statistics:")
        st.write(filtered_data['Isoelectric_Point'].describe())

        plt.figure(figsize=(10, 6))
        sns.histplot(filtered_data['Molecular_Weight'].dropna(), kde=True, bins=30, color='teal')
        plt.xlabel("Molecular Weight (Da)")
        plt.ylabel("Frequency")
        plt.title("Distribution of Protein Molecular Weight")
        st.pyplot(plt)

        plt.figure(figsize=(10, 6))
        sns.histplot(filtered_data['Isoelectric_Point'].dropna(), kde=True, bins=30, color='orange')
        plt.xlabel("Isoelectric Point (pI)")
        plt.ylabel("Frequency")
        plt.title("Distribution of Protein Isoelectric Point")
        st.pyplot(plt)
    else:
        st.warning("The 'Sequence' column is not present in the dataset.")


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
    df = load_data('moonhumannew.csv')  # Update with your actual file path

    # Display dataset overview
    st.subheader("Dataset Overview of Human Moonlighting Proteins")
    st.dataframe(df)

    # Venn Diagram
    st.subheader("Venn Diagram of Proteins Distribution Across Databases")
    plot_venn_diagram(df)

    # Filtering Section
    filtered_data = filter_proteins(df)
    # Add genes to session state
    st.subheader("DataFrames in Session State")
    dataframes_dict = get_dataframes_from_session()
    if not dataframes_dict:
        st.info("No DataFrames found in session state.")
    else:
        df_names = list(dataframes_dict.keys())
        selected_df_name = st.selectbox("Select a DataFrame for successive analysis:", df_names)
        data = dataframes_dict[selected_df_name]

        st.subheader(f"Displaying: {selected_df_name}")
        st.dataframe(data)
    if 'run_analysis' not in st.session_state:
        st.session_state.run_analysis = False
    if st.session_state.run_analysis:
        # Feature Analysis
        feature_analysis(data)

        # Sequence Analysis
        sequence_analysis(data)

        # =======================
        # InterPro Data Analysis Section
        st.subheader("InterPro Data Analysis")

        # Function to safely evaluate strings to Python literals
        def safe_literal_eval(val):
            try:
                return ast.literal_eval(val)
            except (ValueError, SyntaxError):
                return val

        # Convert the 'interpro_data' column from string to actual tuples
        data['interpro_data'] = data['interpro_data'].apply(safe_literal_eval)
        # Convert the 'interpro_data' column from string to actual tuples

        # Convert each list in the 'interpro_data' column to a set to remove duplicates within each row
        data['interpro_data'] = data['interpro_data'].apply(lambda x: list(set(x)) if isinstance(x, list) else x)

        # Explode the 'interpro_data' column
        df_exploded = data.explode('interpro_data')

        # Count the occurrences of each value in the 'interpro_data' column
        most_frequent_accessions = df_exploded['interpro_data'].value_counts()

        # Display the most frequent values
        st.write("### Most Frequent InterPro Data Values:")
        st.dataframe(most_frequent_accessions)


if __name__ == "__main__":
    main()


