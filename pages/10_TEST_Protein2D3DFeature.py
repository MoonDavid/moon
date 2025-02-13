import streamlit as st
import streamlit.components.v1 as components
from stmol import showmol, render_pdb, render_pdb_resn, render_pdb_resi
import py3Dmol
import json
from Bio import SeqIO
import random

# Page title and description
PAGE_TITLE = "Protein Viewer"
PAGE_DESCRIPTION = 'This app allows you to visualize protein sequences and structures. The "Primary Sequence" tab displays the primary sequence of a protein with RNA binding propensity. The "Protein Structure" tab displays the 3D structure of a protein with various visualization options.'
def initialize_page() -> None:
    """Initialize the Streamlit page configuration and header."""
    st.set_page_config(page_title=PAGE_TITLE, layout="wide")
    st.title(PAGE_TITLE)
    st.markdown(PAGE_DESCRIPTION)
def parse_protein_data(file_path):
    """
    Parse protein data file and return a dictionary with protein information
    Returns: dict with structure {protein_id: {
        'sequence': str,
        'rna_propensity': list[float],
        'rna_binary': list[int]
    }}
    """
    proteins = {}
    current_id = None

    try:
        with open(file_path, "r") as f:
            lines = f.readlines()

        i = 0
        while i < len(lines):
            line = lines[i].strip()

            if line.startswith(">"):
                current_id = line[1:]  # Remove '>' character
                proteins[current_id] = {
                    "sequence": "",
                    "rna_propensity": [],
                    "rna_binary": [],
                }
                i += 1

                # Get sequence
                if i < len(lines):
                    proteins[current_id]["sequence"] = lines[i].strip()
                    i += 1

                # Get RNA Propensity
                if i < len(lines) and lines[i].startswith("RNA_Propensity:"):
                    propensity_values = (
                        lines[i].replace("RNA_Propensity:", "").strip().split()
                    )
                    proteins[current_id]["rna_propensity"] = [
                        float(x) for x in propensity_values
                    ]
                    i += 1

                # Get RNA Binary
                if i < len(lines) and lines[i].startswith("RNA_Binary:"):
                    binary_values = lines[i].replace("RNA_Binary:", "").strip().split()
                    proteins[current_id]["rna_binary"] = [int(x) for x in binary_values]
                    i += 1
            else:
                i += 1

    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return {}
    except Exception as e:
        print(f"Error parsing file: {e}")
        return {}

    return proteins





def create_feature_viewer(sequence, rna_probabilities):
    # Create HTML with feature viewer
    assert len(sequence) == len(rna_probabilities), "Length of sequence and probabilities should be the same"
    data_demo = [{"x": i * 2, "y": rna_probabilities[i-1]} for i in range(1, len(sequence) + 1)]

    # Convert the Python list to a JSON string
    data_json = json.dumps(data_demo)
    html_code = f"""
    <!DOCTYPE html>
    <html>
      <head>
        <link rel="stylesheet" type="text/css" href="https://cdn.jsdelivr.net/gh/calipho-sib/feature-viewer@v1.1.0/dist/feature-viewer.bundle.css">
        <script src="https://cdn.jsdelivr.net/gh/calipho-sib/feature-viewer@v1.1.0/dist/feature-viewer.bundle.js"></script>
        <style>
          body {{ margin: 0; padding: 10px; }}
        </style>
      </head>
      <body>
        <div id="fv1"></div>
        <script>
          var proteinSeq = "{sequence}";

          var ft = new FeatureViewer.createFeature(proteinSeq, "#fv1", {{
            showAxis: true,
            showSequence: true,
            brushActive: true,
            toolbar: true,
            bubbleHelp: true,
            zoomMax: 50
          }});

          // Create RNA binding probability data
          var rnaBindingData = [];
          var probabilities = {data_json};

          ft.addFeature({{
            data: probabilities,
            name: "RNA Binding Prob",
            className: "rnabinding",
            color: "#008B8D",
            type: "line",
            height: "5"
          }});
        </script>
      </body>
    </html>
    """
    return html_code



def main():
    # Title of the app
    st.title("Protein Feature Viewer")

    # Protein ID input

    # Two tabs for Primary Sequence and Protein Structure
    tab_seq, tab_structure = st.tabs(["Primary Sequence", "Protein Structure"])

    with tab_seq:
        st.subheader("Primary Sequence Feature Viewer")

        proteins = parse_protein_data("data/batch1.txt")

        # Create selectbox with protein IDs
        protein_ids = list(proteins.keys())
        selected_protein = st.selectbox("Select protein:", protein_ids)

        if selected_protein:
            protein_data = proteins[selected_protein]
            sequence = protein_data["sequence"]
            rna_probabilities = protein_data["rna_propensity"]

            # Create and display the feature viewer
            html_code = create_feature_viewer(sequence, rna_probabilities)
            components.html(html_code, height=500, scrolling=True)

        # [Your existing HTML code here]

    with tab_structure:
        st.subheader("3D Protein Structure Viewer")

        # Create columns for controls and viewer
        col1, col2 = st.columns([1, 2])

        with col1:
            # Create a form for all the controls
            with st.form("structure_controls"):
                st.subheader("Visualization Controls")

                # Style selector
                style = st.selectbox(
                    "Select Style",['cartoon','line','cross','stick','sphere']
                )

                # Color scheme
                color_scheme = st.selectbox(
                    "Color Scheme", ["spectrum", "hydrophobicity", "chain", "residue", "secondary structure"]
                )

                # Background color
                bg_color = st.color_picker("Background Color", "#FFFFFF")

                # Surface options
                show_surface = st.checkbox("Show Surface")
                if show_surface:
                    surface_opacity = st.slider("Surface Opacity", 0.0, 1.0, 0.5)
                    surface_color = st.color_picker("Surface Color", "#FFFFFF")

                # Residue highlighting options
                highlight_type = st.radio(
                    "Highlight by:", ["None", "Residue Name", "Residue Position"]
                )

                if highlight_type == "Residue Name":
                    residue_names = st.text_input(
                        "Enter residue names (comma-separated)", "ALA,ARG,LYS"
                    )
                    highlight_color = st.color_picker("Highlight Color", "#FF0000")

                elif highlight_type == "Residue Position":
                    residue_positions = st.text_input(
                        "Enter residue positions (e.g., '42-44,48,49')", "42-44,48,49"
                    )
                    highlight_color = st.color_picker("Highlight Color", "#FF0000")

                # Label options
                show_labels = st.checkbox("Show Labels")
                if show_labels:
                    label_type = st.multiselect(
                        "Label Information",
                        ["Residue Name", "Residue Number", "Atom Name", "Chain"],
                    )

                # Submit button
                submit_button = st.form_submit_button("Update Visualization")

        with col2:
            protein_id = st.text_input("Enter Protein ID or PDB Code", "1A2C")

            if submit_button:
                # Initialize viewer
                viewer = render_pdb(protein_id)

                # Set style
                if color_scheme == "spectrum":
                    style_options = {"color": color_scheme}
                else:
                    style_options = {"colorscheme": color_scheme}
                viewer.setStyle({style: style_options})


                # Set background
                viewer.setBackgroundColor(bg_color)

                # Add surface if selected
                if show_surface:
                    viewer.addSurface(
                        py3Dmol.VDW, {"opacity": surface_opacity, "color": surface_color}
                    )

                # Handle highlighting
                if highlight_type == "Residue Name":
                    residue_list = [r.strip() for r in residue_names.split(",")]
                    if residue_list:
                        viewer = render_pdb_resn(viewer, residue_list)
                        viewer.setStyle(
                            {"resn": residue_list}, {"stick": {"color": highlight_color}}
                        )

                elif highlight_type == "Residue Position":
                    position_list = [p.strip() for p in residue_positions.split(",")]
                    if position_list:
                        viewer = render_pdb_resi(viewer, position_list)
                        viewer.setStyle(
                            {"resi": position_list}, {"stick": {"color": highlight_color}}
                        )

                # Add labels if selected
                if show_labels:
                    label_dict = {}
                    if "Residue Name" in label_type:
                        label_dict["resn"] = True
                    if "Residue Number" in label_type:
                        label_dict["resi"] = True
                    if "Atom Name" in label_type:
                        label_dict["atom"] = True
                    if "Chain" in label_type:
                        label_dict["chain"] = True

                   # viewer.addResLabels(label_dict)


                st.write(dir(viewer))
                # Show the molecule
                showmol(viewer, height=500, width=800)
            else:
                st.info("Please submit the form to visualize the structure")

if __name__ == "__main__":
    initialize_page()
    main()