import numpy as np
import plotly.express as px
import streamlit as st
from gprofiler import GProfiler


def run_gprofiler_profile():
    """Run Gene Ontology enrichment analysis using g:Profiler."""
    st.header("Functional Profiling with g:Profiler")
    st.info("""
    **What is gProfiler?**
    gProfiler is a powerful web tool for functional enrichment analysis and identifier conversion for genes and proteins.
    Use the sidebar to switch between:
    - **Profile**: Analyze your gene list for enrichment(ORA, Over Rappresentation Analysis) using g:Profiler
    - **Convert**: Transform gene/protein IDs between different naming conventions and databases
    """)

    st.write("""
    This application allows you to perform Gene Ontology enrichment analysis using 
    [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost).
    """)

    # GO Namespace Options
    go_options = {
        "GO:BP (Biological Process)": "GO:BP",
        "GO:MF (Molecular Function)": "GO:MF",
        "GO:CC (Cellular Component)": "GO:CC",
        "Reactome Pathways": "REAC",
        "KEGG Pathways": "KEGG",
        "WikiPathways": "WP",
        "Transcription Factor Targets (TRANSFAC)": "TF",
        "Human Protein Atlas": "HPA",
        "CORUM": "CORUM",
        "Human Phenotype Ontology": "HP",
        "miRTarBase": "MIRNA"
    }

    default_selected = [
        "GO:BP (Biological Process)",
        "GO:MF (Molecular Function)",
        "GO:CC (Cellular Component)"
    ]

    selected_options = []

    with st.form(key='enrichment_form_gProfiler'):
        st.write("### Select the Gene List and GO Namespaces:")

        # Select the list of genes from session state
        available_gene_lists = {key: value for key, value in st.session_state['gene_lists'].items() if isinstance(value, list)}
        if available_gene_lists:
            # Sort the gene lists to show the oldest first based on insertion order
            sorted_gene_list_keys = list(available_gene_lists.keys())
            selected_gene_list = st.selectbox(
                "Select the gene list:",
                options=sorted_gene_list_keys,
                format_func=lambda x: x if x else "No list available"
            )
        else:
            st.warning("No gene list available in session state.")
            return

        col1, col2 = st.columns(2)

        # Display checkboxes for GO namespaces
        for idx, (label, code) in enumerate(go_options.items()):
            if idx % 2 == 0:
                with col1:
                    if st.checkbox(label, value=(label in default_selected), key=f"gprofiler_checkbox_{idx}"):
                        selected_options.append(label)
            else:
                with col2:
                    if st.checkbox(label, value=(label in default_selected), key=f"gprofiler_checkbox_{idx}"):
                        selected_options.append(label)

        st.markdown("---")

        # Top-N Selection
        top_n = st.selectbox(
            "Select the number of top GO terms to visualize:",
            options=[10, 20, 40],
            index=0  # Default to Top 10
        )

        submit_button = st.form_submit_button(label='Run Functional Enrichment')

    if submit_button:
        genes = available_gene_lists.get(selected_gene_list, [])

        if not genes:
            st.warning("Please select a valid list of gene identifiers.")
            return

        go_sources = [go_options[option] for option in selected_options]
        gp = GProfiler(return_dataframe=True)

        with st.spinner("Running GO enrichment. Please wait..."):
            try:
                results = gp.profile(
                    organism='hsapiens',
                    query=genes,
                    sources=go_sources if go_sources else None
                )
                # Store results in session state for later use
                st.session_state['enrichment_results'] = results
                # Store the current top_n selection
                st.session_state['last_top_n'] = top_n
            except Exception as e:
                st.error(f"Error during enrichment analysis: {e}")
                st.write(genes)
                st.session_state['enrichment_results'] = None

        results = st.session_state.get('enrichment_results')

        if results is not None:
            if results.empty:
                st.warning("No enrichment results found for the provided input.")
            else:
                st.success("Enrichment analysis completed successfully!")
                st.dataframe(results)

                csv_results = results.to_csv(index=False).encode('utf-8')
                st.download_button(
                    label="Download Complete Results as CSV",
                    data=csv_results,
                    file_name="gprofiler_enrichment.csv",
                    mime="text/csv"
                )

                # Generate and display the plot based on top_n
                top_results = results.head(top_n)  # Use head() since results are already sorted by p-value

                # Apply logarithmic transformation to p-values
                # Adding a small constant to avoid log(0)
                top_results = top_results.copy()  # Avoid SettingWithCopyWarning
                top_results['log_p_value'] = -np.log10(top_results['p_value'] + 1e-300)

                fig = px.bar(
                    top_results,
                    x='log_p_value',
                    y='name',
                    orientation='h',
                    color='log_p_value',  # Use log_p_value for color mapping
                    color_continuous_scale='Viridis',  # Choose a color scale that represents magnitude effectively
                    labels={
                        'name': 'GO Term',
                        'log_p_value': '-Log10(P-Value)'
                    },
                    title=f'Top {top_n} Enriched GO Terms by P-Value',
                    hover_data={
                        'description': True,
                        'log_p_value': True,
                        'p_value': True
                    }
                )
                fig.update_layout(
                    yaxis={'categoryorder': 'total ascending'},
                    coloraxis_colorbar=dict(
                        title="-Log10(P-Value)",
                        tickvals=[top_results['log_p_value'].min(), top_results['log_p_value'].max()],
                        ticktext=[f"{top_results['log_p_value'].min():.2f}", f"{top_results['log_p_value'].max():.2f}"]
                    )
                )
                st.plotly_chart(fig, use_container_width=True)
        else:
            st.warning("No enrichment results available to display.")

    # Optionally, display the plot if results are already present in session state
    elif 'enrichment_results' in st.session_state and st.session_state['enrichment_results'] is not None:
        results = st.session_state['enrichment_results']
        if not results.empty:
            # Retrieve the last selected top_n from the form or default to 10
            top_n = st.session_state.get('last_top_n', 10)

            top_results = results.head(top_n)  # Use head() since results are already sorted by p-value

            # Apply logarithmic transformation to p-values
            top_results = top_results.copy()  # Avoid SettingWithCopyWarning
            top_results['log_p_value'] = -np.log10(top_results['p_value'] + 1e-300)

            fig = px.bar(
                top_results,
                x='log_p_value',
                y='name',
                orientation='h',
                color='log_p_value',  # Use log_p_value for color mapping
                color_continuous_scale='Viridis',  # Choose a color scale that represents magnitude effectively
                labels={
                    'log_p_value': '-Log10(P-Value)',
                    'name': 'GO Term'
                },
                title=f'Top {top_n} Enriched GO Terms by P-Value',
                hover_data={
                    'description': True,
                    'log_p_value': True,
                    'p_value': True
                }
            )
            fig.update_layout(
                yaxis={'categoryorder': 'total ascending'},
                coloraxis_colorbar=dict(
                    title="-Log10(P-Value)",
                    tickvals=[top_results['log_p_value'].min(), top_results['log_p_value'].max()],
                    ticktext=[f"{top_results['log_p_value'].min():.2f}", f"{top_results['log_p_value'].max():.2f}"]
                )
            )
            st.plotly_chart(fig, use_container_width=True)


def run_gprofiler_convert():
    """Run gene identifier conversion using g:Profiler."""
    st.header("Gene Identifier Conversion with g:Profiler")

    st.info("""
    **What is gProfiler?**
    gProfiler is a powerful web tool for functional enrichment analysis and identifier conversion for genes and proteins.
    Use the sidebar to switch between:
    - **Profile**: Analyze your gene list for enrichment(ORA, Over Rappresentation Analysis) using g:Profiler
    - **Convert**: Transform gene/protein IDs between different naming conventions and databases
    """)
    st.write("""
    This application allows you to convert gene identifiers using 
    [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost).
    """)

    with st.form(key='convert_form'):
        st.write("### Select Conversion Parameters:")

        # Select the list of genes from session state
        available_gene_lists = {key: value for key, value in st.session_state['gene_lists'].items() if
                                isinstance(value, list)}
        # copy


        if available_gene_lists:
            sorted_gene_list_keys = list(available_gene_lists.keys())
            selected_gene_list = st.selectbox(
                "Select the gene list to convert:",
                options=sorted_gene_list_keys,
                format_func=lambda x: x if x else "No list available"
            )
        else:
            st.warning("No gene list available in session state.")
            return

        # Target Namespace Selection
        target_namespace_options = {
            "BioGRID": "BIOGRID",
            "ChEMBL": "CHEMBL",
            "Codelink": "CODELINK_CODELINK",
            "DBASS3 Accession": "DBASS3_ACC",
            "DBASS5 Accession": "DBASS5_ACC",
            "EMBL": "EMBL",
            "Ensembl LRG Gene": "ENS_LRG_GENE",
            "Ensembl LRG Transcript": "ENS_LRG_TRANSCRIPT",
            "Ensembl Gene ID": "ENSG",
            "Ensembl Protein ID": "ENSP",
            "Ensembl Transcript ID": "ENST",
            "Entrez Gene": "ENTREZGENE",
            "Entrez Gene Accession": "ENTREZGENE_ACC",
            "Gene Ontology": "GO",
            "GeneCards": "GENECARDS",
            "GeneCards Accession": "GENECARDS_ACC",
            "GO Slim GOA": "GOSLIM_GOA",
            "HGNC": "HGNC",
            "HGNC Accession": "HGNC_ACC",
            "HGNC Transcript Name": "HGNC_TRANS_NAME",
            "Human Protein Atlas": "HPA",
            "Human Protein Atlas Accession": "HPA_ACC",
            "Illumina HumanRef 8 v3": "ILLUMINA_HUMANREF_8_V3",
            "Illumina HumanWG 6 v3": "ILLUMINA_HUMANWG_6_V3",
            "MEROPS": "MEROPS",
            "MIM Gene": "MIM_GENE",
            "MIM Gene Accession": "MIM_GENE_ACC",
            "MIM Morbid": "MIM_MORBID",
            "MIM Morbid Accession": "MIM_MORBID_ACC",
            "RNAcentral": "RNACENTRAL",
            "RFAM": "RFAM",
            "RFAM Accession": "RFAM_ACC",
            "RFAM Transcript Name": "RFAM_TRANS_NAME",
            "Reactome": "REACTOME",
            "Reactome Gene": "REACTOME_GENE",
            "Reactome Transcript": "REACTOME_TRANSCRIPT",
            "RefSeq mRNA": "REFSEQ_MRNA",
            "RefSeq mRNA Accession": "REFSEQ_MRNA_ACC",
            "RefSeq mRNA Predicted": "REFSEQ_MRNA_PREDICTED",
            "RefSeq mRNA Predicted Accession": "REFSEQ_MRNA_PREDICTED_ACC",
            "RefSeq ncRNA": "REFSEQ_NCRNA",
            "RefSeq ncRNA Accession": "REFSEQ_NCRNA_ACC",
            "RefSeq ncRNA Predicted": "REFSEQ_NCRNA_PREDICTED",
            "RefSeq ncRNA Predicted Accession": "REFSEQ_NCRNA_PREDICTED_ACC",
            "RefSeq Peptide": "REFSEQ_PEPTIDE",
            "RefSeq Peptide Predicted": "REFSEQ_PEPTIDE_PREDICTED",
            "RefSeq Peptide Predicted Accession": "REFSEQ_PEPTIDE_PREDICTED_ACC",
            "UCSC": "UCSC",
            "UniParc": "UNIPARC",
            "UniProt Gene Name": "UNIPROT_GN",
            "UniProt Gene Name Accession": "UNIPROT_GN_ACC",
            "UniProt Isoform": "UNIPROT_ISOFORM",
            "UniProtKB/Swiss-Prot": "UNIPROTSWISSPROT",
            "UniProtKB/Swiss-Prot Accession": "UNIPROTSWISSPROT_ACC",
            "UniProtKB/TrEMBL": "UNIPROTSPTREMBL",
            "UniProtKB/TrEMBL Accession": "UNIPROTSPTREMBL_ACC",
            "WikiGene": "WIKIGENE",
            "WikiGene Accession": "WIKIGENE_ACC",
            "miRBase": "MIRBASE",
        }

        default_target = "ENTREZGENE"

        target_namespace = st.selectbox(
            "Select the target namespace for conversion:",
            options=list(target_namespace_options.keys()),
            index=list(target_namespace_options.keys()).index("Entrez Gene")
        )

        submit_button = st.form_submit_button(label='Run Identifier Conversion')

    if submit_button:
        genes = available_gene_lists.get(selected_gene_list, [])

        if not genes:
            st.warning("Please select a valid list of gene identifiers.")
            return

        target_ns_code = target_namespace_options.get(target_namespace, "ENTREZGENE")

        gp = GProfiler(return_dataframe=True)

        with st.spinner("Running identifier conversion. Please wait..."):
            try:
                results = gp.convert(
                    organism='hsapiens',
                    query=genes,
                    target_namespace=target_ns_code
                )
                # Store results in session state for later use
                st.session_state['conversion_results'] = results
            except Exception as e:
                st.error(f"Error during identifier conversion: {e}")
                st.session_state['conversion_results'] = None

        results = st.session_state.get('conversion_results')

        if results is not None:
            if results.empty:
                st.warning("No conversion results found for the provided input.")
            else:
                st.success("Identifier conversion completed successfully!")
                st.dataframe(results)

                csv_results = results.to_csv(index=False).encode('utf-8')
                st.download_button(
                    label="Download Conversion Results as CSV",
                    data=csv_results,
                    file_name="gprofiler_conversion.csv",
                    mime="text/csv"
                )
        else:
            st.warning("No conversion results available to display.")

    # Optionally, display the conversion results if already present in session state
    elif 'conversion_results' in st.session_state and st.session_state['conversion_results'] is not None:
        results = st.session_state['conversion_results']
        if not results.empty:
            st.dataframe(results)
            csv_results = results.to_csv(index=False).encode('utf-8')
            st.download_button(
                label="Download Conversion Results as CSV",
                data=csv_results,
                file_name="gprofiler_conversion.csv",
                mime="text/csv"
            )


def main():
    st.set_page_config(page_title="gProfiler Enrichment Tool", layout="wide")

    st.sidebar.header("gProfiler Enrichment Tool")
    if "gene_lists" not in st.session_state:
        st.warning(
            "No pre-loaded MPs list found in session state. Please load first main page Moonlighting Proteins Analysis to add some."
        )
        st.stop()
    # Add a radio button in the sidebar to select between 'Profile' and 'Convert'
    service = st.sidebar.radio(
        "Select gProfiler Service",
        options=["Profile", "Convert"],
        help="Choose whether to perform functional enrichment analysis or convert gene identifiers."
    )

    if service == "Profile":
        run_gprofiler_profile()
    elif service == "Convert":
        run_gprofiler_convert()


if __name__ == "__main__":
    main()
