import streamlit as st
import streamlit.components.v1 as components
def main():
    st.set_page_config(page_title="Other tools", page_icon="📈", layout="wide")
    st.sidebar.header("Other tools")
    st.subheader("DMRpred")
    st.write(
        "[DMRpred](http://biomine.cs.vcu.edu/servers/DMRpred/) is a web server that predicts disordered moonlighting regions (DMRs) within protein sequences. The DMRpred server is designed to predict disordered moonlighting regions (DMRs) within protein sequences. For each residue in the input protein sequence, the server generates a numerical score that quantifies its propensity to form a disordered moonlighting region. A higher propensity score indicates a greater likelihood that the residue is part of a disordered moonlighting region. Additionally, the server provides binary annotations, determining whether each residue is predicted to be involved in a disordered moonlighting function or not."
    )

    # Embed the DMRpred website using an iframe
    st.subheader("MEME")
    st.write("[MEME Suite](https://meme-suite.org/meme/) is a comprehensive collection of bioinformatics tools designed for the discovery and analysis of sequence motifs in DNA, RNA, and protein sequences. "
    )
    components.iframe(
        src="https://meme-suite.org/meme/",
        width=1000,
        height=600,
        scrolling=True,
    )

if __name__ == "__main__":
    main()
