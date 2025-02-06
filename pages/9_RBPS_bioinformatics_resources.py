import streamlit as st
import streamlit.components.v1 as components

# Set the page configuration
st.set_page_config(page_title="RNA-Binding Protein Databases & Tools", layout="wide")

# Title
st.title("RNA-Binding Protein Databases & Tools")
st.markdown("""
This app provides information on various databases and tools related to RNA-binding proteins.
Click on the expanders below to see details and (when available) an embedded view of the website.
""")

#############################
# DATABASES Section
#############################

st.header("Databases")


# A helper function to create an expander for each database
def database_expander(name, url, doi, description):
    with st.expander(f"**{name}**"):
        st.markdown(f"**URL:** [Visit {name}]({url})")
        st.markdown(f"**DOI:** {doi}")
        st.markdown(f"**Description:** {description}")
        # Embed the webpage in an iframe if possible
        try:
            components.iframe(url, height=600, scrolling=True)
        except Exception as e:
            st.error(f"Could not embed {name} due to: {e}")


# POSTAR3
database_expander(
    name="POSTAR3",
    url="http://111.198.139.65/",
    doi="10.1093/nar/gkab702",
    description=(
        "POSTAR3 is an updated, comprehensive platform focused on post‐transcriptional regulation. "
        "It integrates an extensive collection of RBP binding sites and regulatory features derived from multiple high‐throughput techniques such as CLIP‐seq, Ribo‐seq, structure‐seq, and degradome‐seq. "
        "Covering data from seven species and multiple experimental technologies, POSTAR3 helps users explore the interplay between RBP binding, RNA secondary structure, and genomic variants."
    )
)

# ATtRACT
database_expander(
    name="ATtRACT",
    url="https://attract.cnic.es/",
    doi="10.1093/database/baw035",
    description=(
        "ATtRACT is a manually curated database that focuses on RNA‐binding proteins and their associated binding motifs. "
        "It compiles experimentally validated data from several sources, providing consensus motifs and sequence profiles for many RBPs. "
        "In addition, it offers tools to scan RNA sequences for potential binding sites—making it a practical resource for researchers aiming to predict RNA–RBP interactions."
    )
)

# ENCORI/StarBase
database_expander(
    name="ENCORI/StarBase",
    url="https://rnasysu.com/encori/",
    doi="",
    description=(
        "Originally launched as StarBase, the ENCORI platform decodes RNA interactomes by integrating multiple types of high‐throughput data "
        "(such as CLIP‐seq, PAR‐CLIP, iCLIP, CLASH, and degradome‐seq). It provides comprehensive interaction maps covering miRNA–mRNA, miRNA–lncRNA, protein–RNA, and other RNA–RNA networks. "
        "With features for pan-cancer analysis and functional prediction, ENCORI/StarBase is an invaluable tool for exploring complex post‐transcriptional regulatory networks."
    )
)

st.subheader("Additional Databases")
with st.expander("RBPDB"):
    st.markdown(
        "**RBPDB** compiles experimental observations of RNA-binding sites (both in vitro and in vivo) from the literature for several metazoan species. It includes binding motifs in formats such as position weight matrices and provides a scanning tool for identifying potential RBP sites.")
    st.markdown("**Source:** [PMC.NCBI.NLM.NIH.GOV](https://www.ncbi.nlm.nih.gov/)")

with st.expander("RBP2GO"):
    st.markdown(
        "**RBP2GO** is a pan-species database that aggregates proteome-wide datasets for RNA-binding proteins. It features gene ontology annotations, network analyses, and functional insights, making it useful for both exploratory studies and hypothesis-driven research.")
    st.markdown("**Source:** [RBP2GO.DKFZ.DE](https://rbp2go.dkfz.de/)")

with st.expander("EuRBPDB / RBPWorld"):
    st.markdown(
        "**EuRBPDB / RBPWorld** is one of the most extensive resources available, cataloging hundreds of thousands of RBPs from a wide range of eukaryotic species. It offers rich annotations including RNA-binding domains, expression profiles, and even disease associations.")
    st.markdown("**Source:** [EURBPDB.GZSYS.ORG.CN](http://eurbpdb.gzsys.org.cn/)")

with st.expander("oRNAment"):
    st.markdown(
        "**oRNAment** provides a catalog of putative RBP binding site instances in both coding and non-coding RNAs across various species. It is particularly useful for computational prediction of RNA-binding sites based on available motif data.")
    st.markdown("**Source:** [MOTIFMAP-RNA.ICS.UCI.EDU](https://motifmap-rna.ics.uci.edu/)")

with st.expander("CLIPdb"):
    st.markdown(
        "**CLIPdb** collects binding sites identified by different CLIP-seq technologies. It facilitates the analysis of RBP–RNA interactions by allowing users to browse and download binding site data from multiple species.")
    st.markdown("**Source:** [BMCGENOMICS.BIOMEDCENTRAL.COM](https://bmcgenomics.biomedcentral.com/)")

with st.expander("RBPbase"):
    st.markdown(
        "**RBPbase** integrates high-throughput RNA interactome capture studies and provides a table-based interface for exploring RBP annotations across multiple organisms. It is designed to quickly retrieve and analyze data on RBPs, including their RNA-binding profiles.")

#############################
# TOOLS Section
#############################

st.header("Tools")


def tool_expander(name, url, doi, description):
    with st.expander(f"**{name}**"):
        st.markdown(f"**URL:** [Visit {name}]({url})")
        if doi:
            st.markdown(f"**PAPER:** [Read {name} paper] ({doi})")
        st.markdown(f"**Description:** {description}")
        try:
            components.iframe(url, height=600, scrolling=True)
        except Exception as e:
            st.error(f"Could not embed {name} due to: {e}")


# RBPsuite
tool_expander(
    name="RBPsuite",
    url="http://www.csbio.sjtu.edu.cn/bioinf/RBPsuite/",
    doi="",
    description=(
        "RBPsuite contains deep learning-based methods for predicting RBP binding sites on RNAs. "
        "It offers two main approaches: iDeepS for linear RNAs and iDeepC for circular RNAs. "
        "These methods are designed to capture the differences in binding preferences of RBPs to different RNA types."
    )
)

# Cis-BP-RNA
tool_expander(
    name="Cis-BP-RNA",
    url="http://cisbp-rna.ccbr.utoronto.ca/index.php",
    doi="",
    description=(
        "Cis-BP-RNA is the Catalog of Inferred Sequence Binding Proteins of RNA. "
        "It is a library of RNA binding protein motifs and specificities, organized for ease of searching, browsing, and downloading. "
        "It also provides built-in web tools for scanning sequences and predicting binding motifs."
    )
)

# EnrichRBP
tool_expander(
    name="EnrichRBP",
    url="https://airbp.aibio-lab.com/app/api/home/index/",
    doi="",
    description=(
        "EnrichRBP is an automated and interpretable platform for the prediction and analysis of RBP binding sites across circRNAs, linear RNAs, "
        "and RNAs in various cellular contexts. It supports state-of-the-art models and provides comprehensive result visualization and model interpretability."
    )
)

# PST-PRNA
tool_expander(
    name="PST-PRNA",
    url="http://www.zpliulab.cn/PSTPRNA",
    doi="10.1093/bioinformatics/btac078",
    description=(
        "PST-PRNA is a novel model for predicting RNA-binding sites based on protein surface topography (PST). "
        "It uses 3D structural information and deep learning to learn latent structural features of protein surfaces and predict binding sites."
    )
)

#############################
# Footer
#############################

st.markdown("---")
st.markdown("Developed with Streamlit.")

