import streamlit as st
import requests
import pandas as pd
from urllib.parse import urlencode

@st.cache_data
def run_panther_analysis(
        gene_list: str,
        organism: str,
        annot_dataset: str,
        enrichment_test: str,
        correction: str,
        ref_list: str = None,
        ref_organism: str = None
) -> (pd.DataFrame, str):
    """
    Submit a GET request to PANTHER to run overrepresentation analysis.
    """
    PANTHER_ENDPOINT = "http://pantherdb.org/services/oai/pantherdb/enrich/overrep"

    params = {
        "geneInputList": gene_list,
        "organism": organism,
        "annotDataSet": annot_dataset,
        "enrichmentTestType": enrichment_test,
        "correction": correction,
    }

    if ref_list:
        params["refInputList"] = ref_list
        if ref_organism:
            params["refOrganism"] = ref_organism

    query_string = urlencode(params)
    full_url = f"{PANTHER_ENDPOINT}?{query_string}"

    try:
        response = requests.get(full_url)
        response.raise_for_status()
        data_json = response.json()
    except requests.exceptions.RequestException as e:
        st.error(f"PANTHER request failed: {e}")
        return pd.DataFrame(), full_url
    except ValueError:
        st.error("Unable to parse JSON response from PANTHER.")
        return pd.DataFrame(), full_url

    if "results" in data_json and "result" in data_json["results"]:
        df = pd.DataFrame(data_json["results"]["result"])
        return df, full_url
    else:
        st.warning("No results returned from PANTHER.")
        return pd.DataFrame(), full_url


def panther_overrepresentation_analysis():
    """
    Run PANTHER overrepresentation analysis within a Streamlit form.
    """
    st.header("PANTHER Overrepresentation (GO Enrichment)")

    available_gene_lists = {key: value for key, value in st.session_state['gene_lists'].items() if isinstance(value, list)}
    if not available_gene_lists:
        st.warning("No gene list available in session state.")
        return

    # Create a form for input parameters
    with st.form(key='panther_enrichment_form'):
        # Gene List Selection
        sorted_gene_list_keys = list(available_gene_lists.keys())
        selected_gene_list = st.selectbox(
            "Select the gene list:",
            options=sorted_gene_list_keys,
            format_func=lambda x: x if x else "No list available"
        )

        # Organism Taxon ID
        taxon_id = st.text_input(
            "Organism Taxon ID (e.g., 9606 for Human)",
            "9606",
            key='panther_taxon_id'
        )
        annotation_datasets = {
            "Biological Process": "GO:0008150",
            "Molecular Function": "GO:0003674",
            "Cellular Component": "GO:0005575",
            "PANTHER GO Slim Molecular Function": "ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF",
            "PANTHER GO Slim Biological Process": "ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP",
            "PANTHER GO Slim Cellular Location": "ANNOT_TYPE_ID_PANTHER_GO_SLIM_CC",
            "Protein Class": "ANNOT_TYPE_ID_PANTHER_PC",
            "PANTHER Pathway": "ANNOT_TYPE_ID_PANTHER_PATHWAY",
            "Reactome Pathway": "ANNOT_TYPE_ID_REACTOME_PATHWAY"
        }

        # Annotation Dataset Selection
        annot_dataset_choice = st.selectbox(
            "Annotation Dataset",
            annotation_datasets.keys(),
            index=0,
            key='panther_annot_dataset'
        )
        annot_dataset_id = annotation_datasets[annot_dataset_choice]
        # Enrichment Test Type
        enrichment_test_type = st.selectbox(
            "Enrichment Test Type",
            ["FISHER", "BINOMIAL"],
            key='panther_enrichment_test'
        )

        # Multiple Testing Correction Method
        correction_method = st.selectbox(
            "Multiple Testing Correction",
            ["FDR", "BONFERRONI", "NONE"],
            key='panther_correction'
        )

        # Custom Reference List
        custom_ref_list = st.text_area(
            "Custom Reference List (comma-separated)",
            "",
            key='panther_ref_list'
        ).strip()

        # Reference Organism Taxon ID (conditionally displayed)
        ref_taxon = None
        if custom_ref_list:
            ref_taxon = st.text_input(
                "Reference Organism Taxon ID",
                taxon_id,
                key='panther_ref_taxon'
            )

        # Submit Button
        submit_button = st.form_submit_button(label="Run PANTHER Enrichment")

    # Handle form submission
    if submit_button or "panther_results" in st.session_state:

        with st.spinner("Submitting request to PANTHER. Please wait..."):
            # Prepare the gene list as a comma-separated string
            gene_list_str = ', '.join(available_gene_lists[selected_gene_list])

            results_df,full_url = run_panther_analysis(
                gene_list=gene_list_str,
                organism=taxon_id,
                annot_dataset=annot_dataset_id,
                enrichment_test=enrichment_test_type,
                correction=correction_method,
                ref_list=custom_ref_list if custom_ref_list else None,
                ref_organism=ref_taxon if custom_ref_list else None
            )

        if results_df.empty:
            st.warning("No enrichment results returned from PANTHER or an error occurred.")
            st.write(full_url)
        else:
            st.success("PANTHER Overrepresentation Analysis Completed!")
            if "panther_results" not in st.session_state:
                st.session_state["panther_results"] = results_df
            st.dataframe(st.session_state["panther_results"])

            # Provide a download button for the results
            csv_results = st.session_state["panther_results"].to_csv(index=False).encode('utf-8')
            st.download_button(
                label="Download Full Results as CSV",
                data=csv_results,
                file_name="panther_enrichment_results.csv",
                mime="text/csv"
            )


def main():
    st.set_page_config(page_title="PANTHER Enrichment Tool", page_icon="📈", layout="wide")
    st.sidebar.header("PANTHER Enrichment Tool")

    panther_overrepresentation_analysis()


if __name__ == "__main__":
    main()
