import streamlit as st
import streamlit.components.v1 as components

def initialize():
    st.set_page_config(
        page_title="Expression analysis",
        page_icon="📈",
        layout="wide"
    )
    st.sidebar.header("Expression analysis")
    st.subheader("Expression analysis")
def main():
    initialize()
    

if __name__ == "__main__":
    main()