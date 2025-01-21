import streamlit as st
from Bio import Entrez
import spacy
import pandas as pd
import pickle
from urllib.parse import quote

from Bio.Entrez import api_key

# Set up BioPython/Entrez
Entrez.email = "davide.gotta@gmail.com"
# Entrez.api_key = "YOUR_NCBI_API_KEY"  # Uncomment and set if you have an API key

# Cache the spaCy model to avoid reloading on every run
@st.cache_resource
def load_spacy_model():
    try:
        nlp = spacy.load("en_ner_bionlp13cg_md")
        return nlp
    except OSError:
        st.error("spaCy model 'en_ner_bionlp13cg_md' not found. Please install it using:\n\npython -m spacy download en_ner_bionlp13cg_md")
        st.stop()

# Function to search PubMed
def search_pubmed(search_term, retmax):
    handle = Entrez.esearch(db="pubmed", term=search_term, retmax=retmax)
    search_results = Entrez.read(handle)
    handle.close()
    id_list = search_results.get("IdList", [])
    return id_list

# Function to fetch PubMed abstracts
def fetch_abstracts(id_list):
    abstracts = []
    articles = []
    if id_list:
        fetch_handle = Entrez.efetch(db="pubmed", id=",".join(id_list), retmode="xml")
        records = Entrez.read(fetch_handle)
        fetch_handle.close()

        for pubmed_article in records.get("PubmedArticle", []):
            # Extract MedlineCitation and Article sections
            medline = pubmed_article.get("MedlineCitation", {})
            article = medline.get("Article", {})

            # Extract the article title
            title = article.get("ArticleTitle", "")

            # Extract the abstract text
            abstract_data = article.get("Abstract", {}).get("AbstractText", [])

            # Extract DOI from PubmedData -> ArticleIdList
            pubmed_data = pubmed_article.get("PubmedData", {})
            article_ids = pubmed_data.get("ArticleIdList", [])
            doi = ""
            for id_obj in article_ids:
                # id_obj can be a string or a dictionary-like object with attributes
                if isinstance(id_obj, dict):
                    id_type = id_obj.get("IdType", "").lower()
                    if id_type == "doi":
                        doi = id_obj.get("#text", "")
                        break
                elif isinstance(id_obj, str):
                    # If id_obj is a string, sometimes the DOI can be identified by its format
                    if id_obj.startswith("10."):  # DOI typically starts with '10.'
                        doi = id_obj
                        break

            # Concatenate abstract paragraphs if present
            if abstract_data:
                full_abstract = " ".join(abstract_data)
                articles.append({
                    "title": title,
                    "doi": doi,
                    "abstract": full_abstract
                })
            else:
                # Handle articles without an abstract
                articles.append({
                    "title": title,
                    "doi": doi,
                    "abstract": ""
                })
    return articles

# Function to perform NER on abstracts
def perform_ner(articles, nlp):
    diz = {}
    for article in articles:
        doi = article["doi"]
        if article['abstract']:
            doc = nlp(article['abstract'])
            genes = set([ent.text for ent in doc.ents if ent.label_ in ["GENE_OR_GENE_PRODUCT", "PROTEIN"]])
            diz[doi] = genes
        else:
            diz[doi] = set()
    return diz

# Function to create clickable DOI links
def create_doi_link(doi):
    if doi:
        url = f"https://doi.org/{doi}"
        return url
    else:
        return "N/A"
api_key=None
# Streamlit App Layout
def main():
    st.set_page_config(page_title="Literature mining", layout="wide")
    st.title("PubMed Abstract Search and Name Entity Recognition(NER)")
    st.markdown("""

    ### What does this app do?

    - **Search PubMed**: Enter your search terms to find relevant scientific articles from the PubMed database.
    - **Fetch Abstracts**: Retrieves the abstracts of the articles matching your search criteria.
    - **Perform NER**: Uses Natural Language Processing to identify and extract gene or protein names mentioned in the abstracts.
    - **View and Download Results**: Displays the results in an interactive table with clickable DOI links and allows you to download the data as a CSV file for further analysis.
""")

    st.header("Search Parameters")

    search_term = st.text_input(
        "Search Terms",
        value="moonlighting protein OR multifunctional protein",
        help="Enter keywords or phrases related to your research topic."
    )

    # Replace number_input with slider for retmax
    retmax = st.slider(
        "Maximum Results (retmax)",
        min_value=1,
        max_value=10000,
        value=500,
        step=1,
        help="Specify the maximum number of articles to retrieve (up to 10,000)."
    )

    # Optional input for Entrez API Key
    api_key = st.text_input(
        "Entrez API Key (optional)",
        type="password",
        help="If you have an NCBI API key, enter it here to increase your request limits."
    )

    if api_key:
        Entrez.api_key = api_key

    if st.button("Search"):
        with st.spinner("Searching PubMed..."):
            id_list = search_pubmed(search_term, retmax)
            st.write(f"Found {len(id_list)} articles for the search term: **{search_term}**")

        if id_list:
            with st.spinner("Fetching abstracts..."):
                articles = fetch_abstracts(id_list)
                st.success(f"Fetched {len(articles)} articles.")

            with st.spinner("Performing Named Entity Recognition (NER)..."):
                nlp = load_spacy_model()
                diz = perform_ner(articles, nlp)
                st.success("NER completed.")

            # Prepare DataFrame
            data = []
            for article in articles:
                doi_link = create_doi_link(article["doi"])
                entities = ", ".join(diz.get(article["doi"], []))
                data.append({
                    "DOI": doi_link,
                    "Title": article["title"],
                    "Entities": entities
                })

            df = pd.DataFrame(data)

            # Display DataFrame with DOI as clickable links
            st.markdown("### Search Results")
            st.dataframe(df,column_config={
            "DOI": st.column_config.LinkColumn(
                "Article",
                help="Link to scientific article"
             )})


            # Allow users to download the results
            csv = df.to_csv(index=False)
            st.download_button(
                label="Download Results as CSV",
                data=csv,
                file_name='pubmed_search_results.csv',
                mime='text/csv',
            )
        else:
            st.warning("No articles found for the given search terms.")

if __name__ == "__main__":
    main()
