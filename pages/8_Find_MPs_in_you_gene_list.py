import streamlit as st
import pandas as pd
from collections import defaultdict
from mygene import MyGeneInfo
from gprofiler import GProfiler
from detect_delimiter import detect
from venny4py.venny4py import *
if "gene_lists" not in st.session_state:
    st.warning(
        "No pre-loaded MPs list found in session state. Please load first main page Moonlighting Proteins Analysis to add some."
    )
    st.stop()
# Retrieve session_state items (make sure these are set somewhere in your app)
DF = st.session_state.get('df', None)
HUMANMPs_all = st.session_state['gene_lists'].get('humanMPs_all', None)


def df_from_uniprots(df, uniprots):
    """Create a DataFrame from a list of UniProt IDs."""
    return df[df['UniProtKB-AC'].isin(uniprots)].reset_index(drop=True)


def create_dataframe_from_exploded_dict(converted):
    """
    Creates a Pandas DataFrame from an exploded dictionary.

    Args:
        converted (dict): A dictionary where keys represent query names and
                           values are lists of UniProt IDs.

    Returns:
        pandas.DataFrame: A DataFrame with 'query' and 'uniprot' columns.
    """
    data = []
    for query, uniprot_ids in converted.items():
        for uniprot_id in uniprot_ids:
            data.append({'query': query, 'uniprot': uniprot_id})
    return pd.DataFrame(data)

@st.cache_data
def convert_identifiers_mygene(identifiers, to_id='uniprot'):
    """Convert identifiers using MyGene.info."""
    mg = MyGeneInfo()
    results = mg.querymany(
        identifiers,
        scopes=['ensembl.gene', 'entrezgene', 'uniprot', 'symbol'],
        fields=to_id,
        species='human'
    )
    converted = defaultdict(list)
    for res in results:
        if 'query' in res and to_id in res:
            # If the returned field is a dictionary (as expected), extract Swiss-Prot/TrEMBL
            if isinstance(res[to_id], dict):
                if 'Swiss-Prot' in res[to_id] and res[to_id]['Swiss-Prot']:
                    if isinstance(res[to_id]['Swiss-Prot'], list):
                        converted[res['query']].extend(res[to_id]['Swiss-Prot'])
                    else:
                        converted[res['query']].append(res[to_id]['Swiss-Prot'])

                if 'TrEMBL' in res[to_id] and res[to_id]['TrEMBL']:
                    if isinstance(res[to_id]['TrEMBL'], list):
                        converted[res['query']].extend(res[to_id]['TrEMBL'])
                    else:
                        converted[res['query']].append(res[to_id]['TrEMBL'])

            # If the returned field is a list or string, simply add it
            else:
                converted[res['query']].append(res[to_id])
    return create_dataframe_from_exploded_dict(converted)

@st.cache_data
def convert_identifiers_gprofiler(identifiers, target_namespace='UNIPROT_GN_ACC'):
    """Convert identifiers using g:Profiler."""
    gp = GProfiler(return_dataframe=True)
    df = gp.convert(organism='hsapiens', query=identifiers, target_namespace=target_namespace)
    df = df[['incoming', 'converted']]
    df.columns = ['query', 'uniprot']
    return df


def main():
    # Set up the page configuration
    st.set_page_config(page_title="Moonlighting Proteins Checker", layout="wide")

    # App Title and Description
    st.title("🧬  Moonlighting Proteins Checker")
    st.markdown("## Gene & Protein Identifier Converter and Checker")
    st.markdown(
        """
        This application takes a list of gene/protein identifiers (e.g., Entrez, Ensembl, UniProt) and checks which ones 
        are annotated as moonlighting proteins. It leverages **MyGene.info** and **g:Profiler** to perform identifier conversions.
        """
    )

    # Sidebar with instructions and about info
    st.info(
            """
            1. **Enter** your gene/protein identifiers in the text area below.  
            2. **Select** your input delimiter or choose auto-detection.  
            3. **Click** on the **Convert & Check** button.  
            4. **Review** the conversion results and see which proteins are marked as moonlighting and in which databases.
            """
        )


    # User input for identifiers
    user_input = st.text_area(
        "Enter gene/protein identifiers:",
        height=150,
        placeholder="e.g., BRCA1, TP53, EGFR, ..."
    )

    # Delimiter selection
    separator_hint = st.selectbox(
        "Specify or let the app auto-detect your input separator:",
        ["Auto-detect commas/spaces/newlines", "Commas (,)", "Spaces", "Newlines"]
    )

    if st.button("Convert & Check"):
        raw_input = user_input.strip()
        if not raw_input:
            st.error("Please enter at least one gene/protein identifier.")
            return

        # Process the input based on the selected delimiter option
        if separator_hint == "Commas (,)":
            identifiers = [x.strip() for x in raw_input.split(",") if x.strip()]
        elif separator_hint == "Spaces":
            identifiers = [x.strip() for x in raw_input.split() if x.strip()]
        elif separator_hint == "Newlines":
            identifiers = [x.strip() for x in raw_input.split("\n") if x.strip()]
        else:
            # Auto-detect delimiter
            delimiter = detect(raw_input)
            identifiers = [x.strip() for x in raw_input.split(delimiter) if x.strip()]

        st.write(identifiers)
        # Conversion using MyGene.info
        with st.spinner("Converting identifiers using MyGene.info..."):
            mygene_converted = convert_identifiers_mygene(identifiers)

        # Conversion using g:Profiler
        with st.spinner("Converting identifiers using g:Profiler..."):
            gprofiler_converted = convert_identifiers_gprofiler(identifiers)
        with st.expander("Conversion results"):
            if separator_hint != "Auto-detect commas/spaces/newlines":
                st.info(f"**Auto-detected delimiter:** `{delimiter}`")
        # Display the conversion results side-by-side
            st.markdown("### Identifier Conversion Results")
            col1, col2 = st.columns(2)
            with col1:
                st.subheader("MyGene.info Conversion")
                if not mygene_converted.empty:
                    st.dataframe(mygene_converted)
                else:
                    st.warning("No conversions returned from MyGene.info.")
            with col2:
                st.subheader("g:Profiler Conversion")
                if not gprofiler_converted.empty:
                    st.dataframe(gprofiler_converted)
                else:
                    st.warning("No conversions returned from g:Profiler.")

        # Combine UniProt IDs from both methods


            mygene_uniprot_ids = set(mygene_converted['uniprot'].unique()) if not mygene_converted.empty else set()
            gprofiler_uniprot_ids = set(gprofiler_converted['uniprot'].unique()) if not gprofiler_converted.empty else set()
            merged_uniprot_ids = mygene_uniprot_ids.union(gprofiler_uniprot_ids)

        # Display intersections with the known moonlighting proteins (HUMANMPs_all)
        if HUMANMPs_all is not None:
            human_mps_set = set(HUMANMPs_all)
            mygene_intersection = mygene_uniprot_ids.intersection(human_mps_set)
            gprofiler_intersection = gprofiler_uniprot_ids.intersection(human_mps_set)
            combined_intersection = merged_uniprot_ids.intersection(human_mps_set)

            st.markdown("### MPs found in gene list")
            col1, col2 = st.columns(2)
            with col1:
                st.write(f"**MyGene.info Intersection:** {len(mygene_intersection)} moonlighting protein(s)")
                st.write(mygene_intersection)
            with col2:
                st.write(f"**g:Profiler Intersection:** {len(gprofiler_intersection)} moonlighting protein(s)")
                st.write(gprofiler_intersection)
            # Create columns for metrics
            col1, col2, col3 = st.columns(3)
            percentage_mps = (len(combined_intersection) / len(identifiers)) * 100 if len(identifiers) > 0 else 0
            with col1:
                st.metric(label="Number of MPs", value=len(combined_intersection))

            with col2:
                st.metric(label="Number of Non-MPs", value=len(identifiers) - len(combined_intersection))

            with col3:
                st.metric(label="Percentage of MPs", value=f"{percentage_mps:.2f}%")
            moondb_set = set(DF[DF['MoonDB'] == True]['UniProtKB-AC'])
            moonprot_set = set(DF[DF['MoonProtDB'] == True]['UniProtKB-AC'])
            multitask_set = set(DF[DF['MultiTaskProtDB'] == True]['UniProtKB-AC'])
            genelistuser=combined_intersection | set(range(1,len(identifiers)-len(combined_intersection)+1))
            dict_sets={
                "MoonDB":moondb_set,
                "MoonProt":moonprot_set,
                "MultiTaskProtDB":multitask_set,
                "ListaInput":genelistuser
            }
            fig, ax = plt.subplots(figsize=(8, 8))  # Explicitly create a figure object
            venny4py(sets=dict_sets,asax=ax)  # Pass the `ax` to `venny4py`

            # Display the figure using Streamlit
            st.pyplot(fig, use_container_width=True)

            # Optionally, show detailed table of all mapped proteins if available
            if DF is not None:

                combined_result_df = df_from_uniprots(DF, merged_uniprot_ids)
                with st.expander("View data about all moonlighting proteins mapped"):
                    st.dataframe(combined_result_df)
        else:
            st.error("The HUMANMPs_all list is not loaded or available.")


if __name__ == "__main__":
    main()
