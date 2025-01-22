import streamlit as st
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import io


def get_lists_from_session():
    """Retrieve keys from session_state that contain list objects."""
    lists = []
    for key, value in st.session_state.items():
        if isinstance(value, list):
            lists.append(key)
    return lists


def df_from_uniprots(df, uniprots):
    """
    Filter the DataFrame based on a list of UniProt IDs.

    Args:
        df (pd.DataFrame): The DataFrame to filter.
        uniprots (list): List of UniProt IDs to filter by.

    Returns:
        pd.DataFrame: Filtered DataFrame.
    """
    return df[df['UniProtKB-AC'].isin(uniprots)].reset_index(drop=True)


def main():
    # Configure the Streamlit page
    st.set_page_config(
        page_title="Save Tables and Gene Lists",
        page_icon="📈",
        layout="wide"
    )

    # Sidebar Configuration
    st.sidebar.header("Save Tables and Gene Lists")

    # Main Header
    st.subheader("Save Tables and Gene Lists")

    available_uniprot_lists = get_lists_from_session()
    if not available_uniprot_lists:
        st.warning("No UniProt ID lists available in the session state to save as identifiers.")
    data_type = st.selectbox(
        "Select the type of data you want to save:",
        ("UniProtID list", "Table", "FASTA")
    )
    diz = {'tab': '\t', 'whitespace': ' ', 'newline': '\n', ',': ',', ';': ';', '|': '|'}

    if data_type == "UniProtID list":
        if available_uniprot_lists:
            # Identifier List Selection
            selected_list = st.selectbox(
                "Select the UniProt ID list to save:",
                available_uniprot_lists
            )

            # Delimiter Selection
            delimiter = st.selectbox(
                "Select delimiter for the TXT file:",
                (",", ";", "|", "newline", "whitespace", "tab")
            )
        else:
            st.info("No available UniProt ID lists to select.")

    elif data_type == "Table":
        selected_table = st.selectbox(
            "Select the table to save:",
            available_uniprot_lists
        )
    elif data_type == "FASTA":
        selected_fasta = st.selectbox(
            "Select the table with sequences to save in fasta format:",
            available_uniprot_lists
        )

    # Handle form submission
    if data_type == "UniProtID list":
        if available_uniprot_lists:
            identifier_list = st.session_state[selected_list]
            if identifier_list:
                # Join the list with the selected delimiter
                content = diz[delimiter].join(identifier_list)
                # Provide a download button for the TXT file
                st.download_button(
                    label="Download TXT",
                    data=content,
                    file_name=f"{selected_list}.txt",
                    mime="text/plain"
                )
            else:
                st.error(f"The list '{selected_list}' is empty.")
        else:
            st.error("No UniProt ID lists available to save.")

    elif data_type == "Table":
        # Check if the required data exists in session_state
        if 'humanMPs_all' in st.session_state and selected_table in st.session_state:
            df = df_from_uniprots(
                st.session_state['df'],
                st.session_state[selected_table]
            )
            if not df.empty:
                # Convert DataFrame to CSV
                csv = df.to_csv(index=False)
                # Provide a download button for the CSV file
                st.download_button(
                    label="Download CSV",
                    data=csv,
                    file_name=f"{selected_table}.csv",
                    mime="text/csv"
                )
                # st.success("Successfully prepared the 'filtered_table.csv' file for download.")
            else:
                st.error("The filtered table is empty. Please check your UniProt IDs.")
        else:
            st.error("Required data ('humanMPs_all' and 'uniprots') not found in session state.")
    elif data_type == "FASTA":
        if 'df' in st.session_state and selected_fasta in st.session_state:
            df = df_from_uniprots(
                st.session_state['df'],
                st.session_state[selected_fasta]
            )
            if not df.empty and 'Sequence' in df.columns:

                fasta_sequences = []
                for index, row in df.iterrows():

                    try:
                        seq = Seq(row['Sequence'])
                        record = SeqRecord(seq, id=row['UniProtKB-AC'], description=row['Gene Names']+' '+row['Protein names'])
                        fasta_sequences.append(record)
                    except:
                        st.warning(f"Sequence not valid, skipping row {index} with {row['UniProtKB-AC']}")

                if fasta_sequences:
                    with io.StringIO() as handle:
                        SeqIO.write(fasta_sequences, handle, "fasta")
                        fasta_string = handle.getvalue()
                    st.download_button(
                        label="Download FASTA",
                        data=fasta_string,
                        file_name=f"{selected_fasta}.fasta",
                        mime="text/fasta"
                    )
                else:
                    st.error(f"No valid sequence found for '{selected_fasta}'.")
            else:
                st.error(
                    "The filtered table is empty or doesn't have a 'Sequence' column. Please check your UniProt IDs.")
        else:
            st.error("Required data ('df' and the selected table) not found in session state.")


if __name__ == "__main__":
    main()