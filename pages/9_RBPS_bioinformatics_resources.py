import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import matplotlib.pyplot as plt
import altair as alt
import requests


PAGE_TITLE = "RNA-Binding Protein bioinformatics resources"
PAGE_DESCRIPTION = """
This page provides information on various databases and tools related to RNA-binding proteins.
Use the tabs below to switch between **Databases** and **Tools**.
"""


#############################
# Helper Function
#############################

def resource_expander(name, url=None, server=None, code=None, doi=None, description="", iframe=True):
    """
    Creates an expander displaying resource information with options to embed or link to external pages.

    Parameters:
    - name (str): Name of the resource.
    - url (str): URL to the main database or resource page (optional).
    - server (str): URL to the web server for the resource (optional).
    - code (str): URL to the code repository for the resource (optional).
    - doi (str): DOI link to the resource's associated paper or publication (optional).
    - description (str): A brief description of the resource.
    - iframe (bool): Whether to attempt embedding the URL in an iframe (default is True).
    """
    st.markdown(f"### **{name}**")
    if url:
        st.markdown(f"[Visit {name} database]({url})")
    if server:
        st.markdown(f"[Visit {name} server]({server})")
    elif not url:
        st.warning(f"No web server available for {name}")
    if code:
        st.markdown(f"[View {name} code repo]({code})")
    if doi:
        st.markdown(f"[Read {name} paper]({doi})")

    # Add resource description
    if description:
        st.markdown(f"**Description:** {description}")

    # Attempt to embed the URL in an iframe
    if iframe:
        with st.expander(f"Embed {name}"):
            try:
                components.iframe(url, height=600, scrolling=True)
            except Exception as e:
                st.error(f"Could not embed {name} due to: {e}")
    st.divider()


def initialize_page() -> None:
    """Initialize the Streamlit page configuration and header."""
    st.set_page_config(page_title=PAGE_TITLE, layout="wide")
    st.title(PAGE_TITLE)
    st.markdown(PAGE_DESCRIPTION)

#############################
# Create Tabs
#############################
def main():
    initialize_page()
    tab_databases, tab_tools = st.tabs(["Databases", "Tools"])

    #############################
    # Databases Tab
    #############################
    with tab_databases:
        st.subheader("Databases")

        # POSTAR3
        resource_expander(
            name="POSTAR3",
            url="http://111.198.139.65/",
            doi="https://10.1093/nar/gkab702",
            description=(
                "POSTAR3 is an updated, comprehensive platform focused on post‐transcriptional regulation. "
                "It integrates an extensive collection of RBP binding sites and regulatory features derived from multiple high‐throughput techniques such as CLIP‐seq, Ribo‐seq, structure‐seq, and degradome‐seq. "
                "Covering data from seven species and multiple experimental technologies, POSTAR3 helps users explore the interplay between RBP binding, RNA secondary structure, and genomic variants."
            )
        )

        # ATtRACT
        resource_expander(
            name="ATtRACT",
            url="https://attract.cnic.es/",
            doi="https://10.1093/database/baw035",
            description=(
                "ATtRACT is a manually curated database that focuses on RNA‐binding proteins and their associated binding motifs. "
                "It compiles experimentally validated data from several sources, providing consensus motifs and sequence profiles for many RBPs. "
                "In addition, it offers tools to scan RNA sequences for potential binding sites—making it a practical resource for researchers aiming to predict RNA–RBP interactions."
            )
            ,iframe=False
        )

        # ENCORI/StarBase
        resource_expander(
            name="ENCORI/StarBase",
            url="https://rnasysu.com/encori/",
            doi="https://doi.org/10.1093/nar/gkt1248",
            description=(
                "Originally launched as StarBase, the ENCORI platform decodes RNA interactomes by integrating multiple types of high‐throughput data "
                "(such as CLIP‐seq, PAR‐CLIP, iCLIP, CLASH, and degradome‐seq). It provides comprehensive interaction maps covering miRNA–mRNA, miRNA–lncRNA, protein–RNA, and other RNA–RNA networks. "
                "With features for pan-cancer analysis and functional prediction, ENCORI/StarBase is an invaluable tool for exploring complex post‐transcriptional regulatory networks."
            )
        )

        # RBPDB
        resource_expander(
            name="RBPDB",
            url="http://rbpdb.ccbr.utoronto.ca/",
            doi="https://doi.org/10.1093/nar/gkq1069",
            description=(
                "RBPDB is a comprehensive database of RNA-binding protein (RBP) binding sites and their associated functional annotations. "
                "It includes data from 100+ species, including human, mouse, rat, and zebrafish. "
                "It is a comprehensive resource for researchers interested in the development of novel RBPs and their interactions with other genomic entities."
            )
        )

        # RBP2GO
        resource_expander(
            name="RBP2GO",
            url="https://rbp2go.dkfz.de/",
            doi="https://doi.org/10.1093/nar/gkaa1040",
            description=(
                "RBP2GO is a comprehensive database of RNA-binding protein (RBP) binding sites and their associated functional annotations. "
                "It includes data from 100+ species, including human, mouse, rat, and zebrafish."
            )
        )

        # EuRBPDB / RBPWorld
        resource_expander(
            name="EuRBPDB / RBPWorld",
            url="http://eurbpdb.gzsys.org.cn/",
            doi="https://10.1093/nar/gkab702",
            description=(
                "EuRBPDB is a comprehensive database of RNA-binding protein (RBP) binding sites and their associated functional annotations."
            )
        )

        # oRNAment
        resource_expander(
            name="oRNAment",
            url="https://motifmap-rna.ics.uci.edu/",
            doi="https://doi.org/10.1093/nar/gkz986",
            description=(
                "oRNAment is a catalog of putative RBP binding site instances in both coding and non-coding RNAs across various species. "
                "It is particularly useful for computational prediction of RNA-binding sites based on available motif data."
            )
        )

        # RBPbase
        resource_expander(
            name="RBPbase",
            url="https://apps.embl.de/rbpbase/",
            doi="https://10.1093/nar/gkab702",
            description=(
                "RBPbase integrates high-throughput RNA interactome capture studies and provides a table-based interface for exploring RBP annotations "
                "across multiple organisms. It is designed to quickly retrieve and analyze data on RBPs, including their RNA-binding profiles."
            )
        )

        resource_expander(
            name="RBP Image Database",
            url="https://rnabiology.ircm.qc.ca/RBPImage/",
            doi="https://10.1093/nar/gkac971",
            description=(
                "The RBP Image Database is a comprehensive resource for studying the subcellular localization properties of 301 RNA-binding "
                "proteins (RBPs) in human HepG2 and HeLa cell lines. It features systematic immunofluorescence studies with ∼250,000 microscopy "
                "images and annotations detailing localization features based on organelle and subcellular markers. This database supports the "
                "rapid querying of RBPs through a user-friendly interface and facilitates insights into the functional regulation of RBPs in "
                "specific subcellular regions. The database provides a visual and descriptive framework for researchers exploring the roles of "
                "RNA-binding proteins in RNA metabolism."
            )
        )

    #############################
    # Tools Tab
    #############################
    with tab_tools:

        st.header("Protein-based")
        resource_expander(
            name='MucLiPred',
            code='https://github.com/sethzhangjs/MucLiPred',
            doi='https://doi.org/10.1021/acs.jcim.3c01471',
            description="MucLiPred: Multi-Level Contrastive Learning for Predicting Nucleic Acid Binding Residues of Proteins."
                "\nProtein-molecule interactions play a crucial role in various biological functions, with their accurate prediction "
                "being pivotal for drug discovery. MucLiPred introduces a dual contrastive learning mechanism to improve prediction "
                "performance for such interactions, offering significant advances in identifying molecule-binding residues.",
            iframe=False
        )
        resource_expander(
            name="Pprint2",
            server="https://webs.iiitd.edu.in/raghava/pprint2/index.php",
            code="https://github.com/raghavagps/pprint2",
            doi="10.1093/bib/bbac538",
            description=(
                "Pprint2 is an improved version of Pprint designed for predicting RNA-interacting residues in a protein. "
                "This tool utilizes convolutional neural networks (CNN) and evolutionary profile information to achieve an AUC of 0.82 "
                "with a Matthews correlation coefficient of 0.49 on the validation dataset. The tool outperforms other existing methods "
                "and is available as both a standalone software and a web-based server."
            )
        )

        resource_expander(
            name="iDRNA-ITF",
            server="http://bliulab.net/iDRNA-ITF/",
            doi="https://10.1093/bib/bbac236",
            description=(
                "iDNA-ITF is a novel sequence-based method that identifies DNA- and RNA-binding residues in proteins using an innovative induction and transfer framework. "
                "It uniquely incorporates functional properties of residues through a specialized feature extraction network, followed by separate DNA and RNA prediction networks. "
                "The tool demonstrates state-of-the-art performance across multiple test sets, making it valuable for studying protein functions and supporting drug design efforts."
            )
        )
        resource_expander(
            name="NucBind",
            url="http://yanglab.nankai.edu.cn/NucBind",
            doi="https://doi.org/10.1093/bioinformatics/bty756",
            description=(
                "NucBind is a cutting-edge method designed to predict nucleic acids-binding residues in proteins. It combines an "
                "ab-initio method (SVMnuc) and a template-based method (COACH-D). SVMnuc leverages improved features extracted from three "
                "complementary sequence profiles, while COACH-D utilizes homologous template information from a nucleic acids-binding library. "
                "Through a combination of these complementary approaches, NucBind achieves superior performance compared to other state-of-the-art "
                "methods. Despite higher accuracy, some cross-prediction errors between DNA and RNA-binding residues remain a challenge."
            )
        )

        resource_expander(
            name="DRNApred",
            server="http://biomine.cs.vcu.edu/servers/DRNApred/",
            doi="https://10.1093/nar/gkx059",
            description=(
                "DRNApred is a high-throughput sequence-based method that accurately predicts and discriminates between DNA- and RNA-binding residues in proteins. "
                "Using a novel two-layered architecture and specialized regression that penalizes cross-predictions, it significantly reduces false positives between DNA and RNA binding predictions. "
                "The tool outperforms existing predictors by maintaining high accuracy while minimizing cross-predictions, making it particularly valuable for whole-genome scale applications and identification of novel nucleic acid binding proteins."
            )
        )

        resource_expander(
            name="DeepBind",
            server="http://www.rnainter.org/DeepBind/",
            doi="https://doi.org/10.1038/nbt.3300",
            description="""
            DeepBind is a deep learning-based method for predicting RNA and DNA binding sites in proteins""",
            iframe=False

        )
        resource_expander(
            name="RNAProt",
            code="https://github.com/BackofenLab/RNAProt",
            doi="https://doi.org/10.1093/gigascience/giab054",
            description=(
                "RNAProt is an advanced computational framework for predicting RBP (RNA-binding protein) binding sites using recurrent neural networks. "
                "It offers state-of-the-art predictive performance, superior run-time efficiency, and flexibility with additional features like structure information and user-defined features, "
                "making it broader in capability than existing tools. RNAProt is easy to install, comes with comprehensive documentation, and provides informative visualizations "
                "to aid in understanding binding preferences. [GigaScience, Volume 10, Issue 8, August 2021](https://doi.org/10.1093/gigascience/giab054)."
            ),
            iframe=False
        )

        st.header("RNA-based")
        resource_expander(
            name="iDeepE",
            url="",
            doi="",
            description=(
                "iDeepE is a deep learning framework that integrates local and global sequence features using CNNs for accurate prediction of RNA-protein binding sites. "
                "By combining local and global CNNs, iDeepE improves prediction accuracy, balances computational efficiency, and outperforms state-of-the-art methods (e.g., DeepBind, GraphProt). "
                "Tested on large-scale datasets, iDeepE successfully identifies known motifs and enhances predictions by leveraging both local sequence motifs and global sequence contexts."
            )
        )

        resource_expander(
            name="Reformer",
            server="https://doi.org/10.1016/j.patter.2024.101150",
            doi="https://10.1016/j.patter.2024.101150",
            description=(
                "Reformer is a deep learning model for RNA-protein binding site prediction at single-base resolution. "
                "Trained on 225 eCLIP-seq datasets across 155 RNA-binding proteins, it identifies subtle motifs beyond traditional methods, revealing the effects of mutations on RNA regulation. "
                "Reformer offers high accuracy in binding affinity predictions and aids in prioritizing mutations for disease research and RNA-based therapeutics."
            )
        )

        resource_expander(
            name="iDeep",
            server="http://www.csbio.sjtu.edu.cn/bioinf/iDeep",
            doi="https://doi.org/10.1186/s12859-017-1561-8",
            description=(
                "iDeep is a deep learning-based framework that predicts RNA-protein binding sites and motifs using a hybrid "
                "convolutional neural network and deep belief network. By integrating multi-source data at a high abstraction level, "
                "iDeep achieves enhanced prediction performance and interprets binding motifs effectively, outperforming state-of-the-art predictors."
            )
        )

        resource_expander(
            name="RBPmap",
            server="http://rbpmap.technion.ac.il/",
            doi="10.1093/nar/gkab702",
            description=(
                "RBPmap is a web-based tool for predicting RNA-binding protein binding sites on RNA sequences. "
                "It uses a deep learning-based model to predict binding sites on RNA sequences. "
                "It is designed to capture the differences in binding preferences of RBPs to different RNA types."
            )
        )

        # RBPsuite
        resource_expander(
            name="RBPsuite",
            server="http://www.csbio.sjtu.edu.cn/bioinf/RBPsuite/",
            doi="",
            description=(
                "RBPsuite contains deep learning-based methods for predicting RBP binding sites on RNAs. "
                "It offers two main approaches: iDeepS for linear RNAs and iDeepC for circular RNAs. "
                "These methods are designed to capture the differences in binding preferences of RBPs to different RNA types."
            )
        )

        # Cis-BP-RNA
        resource_expander(
            name="Cis-BP-RNA",
            server="http://cisbp-rna.ccbr.utoronto.ca/index.php",
            doi="",
            description=(
                "Cis-BP-RNA is the Catalog of Inferred Sequence Binding Proteins of RNA. "
                "It is a library of RNA binding protein motifs and specificities, organized for ease of searching, browsing, and downloading. "
                "It also provides built-in web tools for scanning sequences and predicting binding motifs."
            )
        )

        # EnrichRBP
        resource_expander(
            name="EnrichRBP",
            server="https://airbp.aibio-lab.com/app/api/home/index/",
            doi="",
            description=(
                "EnrichRBP is an automated and interpretable platform for the prediction and analysis of RBP binding sites across circRNAs, linear RNAs, "
                "and RNAs in various cellular contexts. It supports state-of-the-art models and provides comprehensive result visualization and model interpretability."
            )
        )

        # PST-PRNA
        resource_expander(
            name="PST-PRNA",
            server="http://www.zpliulab.cn/PSTPRNA",
            doi="https://10.1093/bioinformatics/btac078",
            description=(
                "PST-PRNA is a novel model for predicting RNA-binding sites based on protein surface topography (PST). "
                "It uses 3D structural information and deep learning to learn latent structural features of protein surfaces and predict binding sites."
            )
        )
        st.header("Tools that use both RNA and Protein")
        resource_expander(
            name="RPISeq",
            server="http://pridb.gdcb.iastate.edu/RPISeq/",
            doi="",
            description=(
                "RPISeq predicts the probability of RNA-protein interactions based on sequence (RNA-protein interaction prediction). "
                "This tool takes both RNA and protein sequences as inputs and uses Random Forest (RF) and Support Vector Machine (SVM) classifiers. "
                "Users can input RNA and protein sequences to predict their likelihood of interaction."
            )
        )
        resource_expander(
            name="catRAPID omics v2.0",
            server="http://service.tartaglialab.com/page/catrapid_omics2_group",
            doi="https://doi.org/10.1093/nar/gkab393",
            description=(
                "catRAPID omics v2.0 is a web server designed for the prediction of protein–RNA interaction propensities at the "
                "transcriptome- and RNA-binding proteome-level in 8 model organisms. The server allows predictions for multiple input "
                "protein or RNA sequences and computes interaction scores using updated precompiled libraries. Major features include "
                "predicting interactions between custom protein and RNA sets, analyzing long linear RNAs and circular RNAs, and identifying "
                "RNA-binding motifs in predicted RNA targets of proteins. The server also highlights predicted binding sites and checks "
                "for conservation in orthologous protein–RNA pairs. This updated version enhances the prediction accuracy and expands the "
                "capability to handle a wider range of biologically relevant scenarios."
            )
        )

    st.divider()
    st.header("RNA Binding Region(RBR) Prediction Methods Benchamrk")



    st.subheader("Metrics used")
    # Create two columns for MCC and AUC
    col1, col2 = st.columns(2)

    # Column for Matthews Correlation Coefficient (MCC)
    with col1:
        st.header("Matthews Correlation Coefficient (MCC)")
        st.write("""
        The Matthews Correlation Coefficient is a balanced measure that takes into account true and false positives and negatives.
        It is especially useful for imbalanced datasets.
        """)
        st.latex(r'''
        \text{MCC} = \frac{(TP \times TN) - (FP \times FN)}{\sqrt{(TP+FP)(TP+FN)(TN+FP)(TN+FN)}}
        ''')
        st.write("""
        **Interpretation:**
        - **+1**: Perfect prediction.
        - **0**: No better than random guessing.
        - **-1**: Total disagreement between prediction and observation.
        """)

    # Column for Area Under the ROC Curve (AUC)
    with col2:
        st.header("Area Under the ROC Curve (AUC)")
        st.write("""
        The AUC measures the ability of a classifier to distinguish between classes.
        It is the area under the ROC curve, which plots the True Positive Rate (TPR) against the False Positive Rate (FPR)
        at various threshold settings.
        """)
        st.write("""
        **Interpretation of AUC:**
        - **1.0:** Perfect classifier.
        - **0.5:** No better than random guessing.
        - AUC represents the probability that a randomly chosen positive instance is ranked higher than a randomly chosen negative instance.
        """)
    data = {
        "Method name": [
            "Pprint", "RNABindR-plus", "DRNApred", "HybridNAP", "NucBind",
            "ProNA2020", "NCBRPred", "iDRNA-ITF", "Pprint2", "MucLiPred"
        ],
        "Year published": [2008, 2014, 2017, 2019, 2019, 2020, 2021, 2022, 2023, 2024],
        "AUC": [0.63, 0.73, 0.52, 0.59, 0.74, 0.68, 0.69, 0.77, 0.82, 0.84],
        "MCC": [0.08, 0.15, 0.01, 0.05, 0.14, 0.10, 0.14, 0.19, 0.49, 0.43],
    }

    # Create a DataFrame
    df = pd.DataFrame(data)

    # Streamlit display
    st.title("RBR Prediction Methods")
    st.dataframe(df)

    # Button to generate the plot
    if st.button("Generate Interactive Plot"):
        # Create an Altair chart
        st.scatter_chart(
            df,
            x="AUC",
            y="MCC",
            color="Method name",
        )

if __name__ == "__main__":
    main()
