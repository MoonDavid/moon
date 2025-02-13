import streamlit as st
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import io


def get_lists_from_session():
    """Retrieve keys from session_state that contain list objects."""
    return [key for key, value in st.session_state['gene_lists'].items() if isinstance(value, list)]


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


def configure_page():
    st.set_page_config(
        page_title="Save Tables and Gene Lists",
        page_icon="📈",
        layout="wide"
    )
    st.sidebar.header("Save Tables and Gene Lists")
    st.subheader("Save Tables and Gene Lists")


def handle_uniprot_id_list(diz, available_uniprot_lists):
    if not available_uniprot_lists:
        st.info("No available UniProt ID lists to select.")
        return

    selected_list = st.selectbox(
        "Select the UniProt ID list to save:",
        available_uniprot_lists
    )
    delimiter = st.selectbox(
        "Select delimiter for the TXT file:",
        (",", ";", "|", "newline", "whitespace", "tab")
    )

    identifier_list = st.session_state['gene_lists'].get(selected_list, [])
    if not identifier_list:
        st.error(f"The list '{selected_list}' is empty.")
        return

    content = diz[delimiter].join(identifier_list)
    st.download_button(
        label="Download TXT",
        data=content,
        file_name=f"{selected_list}.txt",
        mime="text/plain"
    )


def handle_table(available_uniprot_lists):
    selected_table = st.selectbox("Select the table to save:", available_uniprot_lists)

    if 'humanMPs_all' in st.session_state['gene_lists']:
        df = df_from_uniprots(
            st.session_state['df'],
            st.session_state['gene_lists'][selected_table]
        )

        if not df.empty:
            csv = df.to_csv(index=False)
            st.download_button(
                label="Download CSV",
                data=csv,
                file_name=f"{selected_table}.csv",
                mime="text/csv"
            )
        else:
            st.error("The filtered table is empty. Please check your UniProt IDs.")
    else:
        st.error("Required data ('humanMPs_all' and 'uniprots') not found in session state.")


def handle_fasta(available_uniprot_lists):
    selected_fasta = st.selectbox(
        "Select the table with sequences to save in FASTA format:",
        available_uniprot_lists
    )

    if 'df' in st.session_state:
        df = df_from_uniprots(
            st.session_state['df'],
            st.session_state['gene_lists'][selected_fasta]
        )

        if not df.empty and 'Sequence' in df.columns:
            include_description = st.checkbox(
                "Include description in FASTA headers", value=True)
            filter_selenocysteine = st.checkbox(
                "Exclude sequences with selenocysteine (U)", value=False)

            fasta_sequences = []
            for index, row in df.iterrows():
                try:
                    seq = Seq(row['Sequence'])

                    if "U" in seq:
                        st.warning(
                            f"Sequence contains selenocysteine (U) at row {index} with UniProt ID: {row['UniProtKB-AC']}")
                        if filter_selenocysteine:
                            continue

                    description = (f"{row['Gene Names']} {row['Protein names']}"
                                   if include_description else "")
                    record = SeqRecord(seq, id=row['UniProtKB-AC'], description=description)
                    fasta_sequences.append(record)
                except Exception as e:
                    st.warning(
                        f"Sequence not valid, skipping row {index} with UniProt ID: {row['UniProtKB-AC']}. Error: {e}")

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
            st.error("The filtered table is empty or doesn't have a 'Sequence' column.")
    else:
        st.error("Required data ('df' and the selected table) not found in session state.")


def main():
    configure_page()
    if 'gene_lists' not in st.session_state:
        st.warning(
            "No pre-loaded MPs list found in session state. Please load first main page Moonlighting Proteins Analysis to add some."
        )
        st.stop()
    available_uniprot_lists = get_lists_from_session()


    data_type = st.selectbox(
        "Select the type of data you want to save:",
        ("UniProtID list", "Table", "FASTA")
    )

    diz = {'tab': '\t', 'whitespace': ' ', 'newline': '\n', ',': ',', ';': ';', '|': '|'}

    if data_type == "UniProtID list":
        handle_uniprot_id_list(diz, available_uniprot_lists)
    elif data_type == "Table":
        handle_table(available_uniprot_lists)
    elif data_type == "FASTA":
        handle_fasta(available_uniprot_lists)


if __name__ == "__main__":
    main()
