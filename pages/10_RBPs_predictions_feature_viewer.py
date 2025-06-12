import streamlit as st
import streamlit.components.v1 as components
from stmol import showmol, render_pdb, render_pdb_resn, render_pdb_resi
import py3Dmol
import json
from Bio import SeqIO
import random
import pandas as pd

# Page title and description
PAGE_TITLE = "Protein Viewer"
PAGE_DESCRIPTION = 'This app allows you to visualize protein sequences and structures. The "Primary Sequence" tab displays the primary sequence of a protein with RNA binding propensity. The "Protein Structure" tab displays the 3D structure of a protein with various visualization options.'
def initialize_page() -> None:
    """Initialize the Streamlit page configuration and header."""
    st.set_page_config(page_title=PAGE_TITLE, layout="wide")
    st.title(PAGE_TITLE)
    st.markdown(PAGE_DESCRIPTION)

def get_alphafold_url(uniprot_id):
    return f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb'

# Fetch AlphaFold2 structure
def fetch_structure(uniprot_id):
    url = get_alphafold_url(uniprot_id)
    response = requests.get(url)
    response.raise_for_status()
    return response.text
def parse_pstprna(file_path, threshold=0.5):
    """
    Parse RNA binding propensity data from PST-PRNA format files.

    Args:
        file_path: Path to the result file
        threshold: Score threshold for binary classification (default: 0.5)

    Returns:
        Dictionary containing sequence, RNA binding propensity, and binary binding sites
    """
    # Three-letter to one-letter amino acid code mapping
    aa_map = {
        "Ala": "A",
        "Arg": "R",
        "Asn": "N",
        "Asp": "D",
        "Cys": "C",
        "Glu": "E",
        "Gln": "Q",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Leu": "L",
        "Lys": "K",
        "Met": "M",
        "Phe": "F",
        "Pro": "P",
        "Ser": "S",
        "Thr": "T",
        "Trp": "W",
        "Tyr": "Y",
        "Val": "V",
    }

    positions = []
    residues = []
    states = []
    scores = []

    with open(file_path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 4:  # Ensure line has all required fields
                positions.append(int(parts[0]))
                residues.append(parts[1])  # Amino acid name
                states.append(parts[2])  # State ('s' or 'i')
                scores.append(float(parts[3]))  # Binding probability

    # Convert residues to one-letter codes and build sequence
    sequence = "".join(aa_map.get(res, "X") for res in residues)

    # Create binary classification based on threshold
    binary = [
        1 if (states[i] == "s" and scores[i] > threshold) else 0
        for i in range(len(scores))
    ]

    return {"sequence": sequence, "RNA_propensity": scores, "RNA_binary": binary}
def parse_iDRNA_TF(file_path):
    proteins = {}
    current_protein = None

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()

            # Detect new protein ID
            if line.startswith(">"):
                current_protein = line[1:].strip()
                proteins[current_protein] = {
                    "sequence": "",  # Will join amino acids later
                    "DNA_propensity": [],
                    "DNA_binary": [],
                    "RNA_propensity": [],
                    "RNA_binary": [],
                }
                continue

            # Skip headers or empty lines
            if line.startswith("Amino") or not line:
                continue

            # Extract data
            fields = line.split()
            amino_acid = fields[0]
            prob_dna = float(fields[1])
            bin_dna = int(fields[2])
            prob_rna = float(fields[3])
            bin_rna = int(fields[4])

            # Append values to respective lists
            proteins[current_protein]["sequence"] += amino_acid  # Build the sequence
            proteins[current_protein]["DNA_propensity"].append(prob_dna)
            proteins[current_protein]["DNA_binary"].append(bin_dna)
            proteins[current_protein]["RNA_propensity"].append(prob_rna)
            proteins[current_protein]["RNA_binary"].append(bin_rna)

    return proteins




def parse_hybridRNAbind(file_path):
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
                    "RNA_propensity": [],
                    "RNA_binary": [],
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
                    proteins[current_id]["RNA_propensity"] = [
                        float(x) for x in propensity_values
                    ]
                    i += 1

                # Get RNA Binary
                if i < len(lines) and lines[i].startswith("RNA_Binary:"):
                    binary_values = lines[i].replace("RNA_Binary:", "").strip().split()
                    proteins[current_id]["RNA_binary"] = [int(x) for x in binary_values]
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



import requests
def create_feature_viewer2(
    uniprot, sequence, rna_binding, interpro_list, mro_info=None, max_points=500
):
    """
    Create a feature viewer for protein sequence with RNA binding propensity and other features.
    Uses data sampling for long sequences to improve performance.
    """

    # Helper function to sample data points for long sequences
    def sample_data_points(values, start=1, max_points=500):
        if len(values) <= max_points:
            return [
                {"x": i, "y": values[i - start]}
                for i in range(start, len(values) + start)
            ]
        st.warning(
            f"Warning: The sequence is too long ({len(values)} points). Data has been sampled for better performance."
        )

        # Sample at regular intervals
        step = max(1, len(values) // max_points)
        sampled_points = []

        for i in range(start, len(values) + start, step):
            # Calculate average value for this segment to preserve pattern
            if i + step <= len(values) + start:
                segment = values[i - start : i + step - start]
                avg_value = sum(segment) / len(segment)
                sampled_points.append({"x": i, "y": avg_value})
            else:
                sampled_points.append({"x": i, "y": values[i - start]})

        return sampled_points

    # Helper function to create binary feature sites
    def create_binary_sites(binary_values, start=1):
        return [
            {"x": i, "y": i, "description": f"Binding site"}
            for i in range(start, len(binary_values) + start)
            if binary_values[i - start] == 1
        ]

    # Generate data for RNA binding probabilities and binary sites
    features_data = {}

    # Process hybridRNAbind data with sampling
    features_data["hybridRNAbind_prob"] = json.dumps(
        sample_data_points(
            rna_binding["hybridRNAbind"]["RNA_propensity"], max_points=max_points
        )
    )
    features_data["hybridRNAbind_binary"] = json.dumps(
        create_binary_sites(rna_binding["hybridRNAbind"]["RNA_binary"])
    )

    # Process iDNA-TF data with sampling
    features_data["iDRNA_TF_prob"] = json.dumps(
        sample_data_points(
            rna_binding["iDRNA-TF"]["RNA_propensity"], max_points=max_points
        )
    )
    features_data["iDRNA_TF_binary"] = json.dumps(
        create_binary_sites(rna_binding["iDRNA-TF"]["RNA_binary"])
    )

    # Fetch protein features from UniProt
    url = f"https://www.ebi.ac.uk/proteins/api/features/{uniprot}"
    try:
        features = requests.get(url).json()
    except:
        features = {"features": []}  # Handle API errors gracefully

    # Process domain and site features
    domain_list = []
    for feature in features.get("features", []):
        if feature["category"] == "DOMAINS_AND_SITES":
            domain_list.append(
                {
                    "x": int(feature.get("begin", 0)),
                    "y": int(feature.get("end", 0)),
                    "description": feature.get("description", ""),
                }
            )
    features_data["domains"] = json.dumps(domain_list)

    # Process PTM features
    ptm_list = []
    for feature in features.get("features", []):
        if feature["category"] == "PTM":
            ptm_list.append(
                {
                    "x": int(feature.get("begin", 0)),
                    "y": int(feature.get("end", 0)),
                    "description": feature.get("description", ""),
                }
            )
    features_data["ptm"] = json.dumps(ptm_list)

    # Process InterPro features
    features_data["interpro"] = json.dumps(interpro_list)

    # Process MRO information if available
    mro_features = []
    if mro_info:
        for role in mro_info:
            mro_features.append(
                {
                    "x": role.get("start", 0),
                    "y": role.get("end", 0),
                    "description": role.get("description", "Moonlighting role"),
                }
            )
    features_data["mro"] = json.dumps(mro_features)

    # Generate feature viewer configuration
    feature_configs = [
        {
            "data": "hybridRNAbind_prob",
            "name": "HybridRNAbind Prob",
            "className": "rnabinding",
            "color": "#008B8D",
            "type": "line",
            "height": "5",
        },
        {
            "data": "iDRNA_TF_prob",
            "name": "iDRNA-TF Prob",
            "className": "rnabinding",
            "color": "#0066CC",
            "type": "line",
            "height": "5",
        },
        {
            "data": "hybridRNAbind_binary",
            "name": "HybridRNAbind Sites",
            "className": "rnasites",
            "color": "#008000",
            "type": "rect",
        },
        {
            "data": "iDRNA_TF_binary",
            "name": "iDRNA-TF Sites",
            "className": "rnasites",
            "color": "#00AA00",
            "type": "rect",
        },
        {
            "data": "domains",
            "name": "Domains and Sites",
            "className": "domains",
            "color": "#FFA500",
            "type": "rect",
        },
        {
            "data": "ptm",
            "name": "PTM",
            "className": "ptm",
            "color": "#FF0000",
            "type": "rect",
        },
        {
            "data": "interpro",
            "name": "InterPro",
            "className": "interpro",
            "color": "#9932CC",
            "type": "rect",
        },
    ]

    # Add MRO features if available
    if mro_info:
        feature_configs.append(
            {
                "data": "mro",
                "name": "Moonlighting Roles",
                "className": "mro",
                "color": "#FF6347",
                "type": "rect",
            }
        )

    # Generate JavaScript feature configurations
    js_features = []
    for config in feature_configs:
        js_features.append(f"""
          ft.addFeature({{
            data: {config["data"]},
            name: "{config["name"]}",
            className: "{config["className"]}",
            color: "{config["color"]}",
            type: "{config["type"]}"{"," if "height" in config else ""}{
            f' height: "{config["height"]}"' if "height" in config else ""
        }
          }});
        """)

    # Add sampling note for long sequences
    sampling_note = ""
    if len(sequence) > max_points:
        sampling_note = f"""
        <div style="background-color: #e8f4f8; padding: 8px; margin-bottom: 10px; border-radius: 4px; font-size: 12px;">
            Note: This protein has {len(sequence)} amino acids. Data has been sampled for better performance.
        </div>
        """

    # Generate the HTML
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
        {sampling_note}
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

          // Data variables
          {"; ".join([f"var {key} = {value}" for key, value in features_data.items()])}

          // Add features
          {"".join(js_features)}
        </script>
      </body>
    </html>
    """

    return html_code
def create_feature_viewer(uniprot, sequence, rna_binding, interpro_list, mro_info=None):
    """
    Create a feature viewer for protein sequence with RNA binding propensity and other features.

    Args:
        uniprot: UniProt ID for the protein
        sequence: Protein sequence
        rna_binding: Dictionary containing RNA binding data from different methods
        interpro_list: List of InterPro features
        mro_info: Optional dictionary containing moonlighting roles information

    Returns:
        HTML code for the feature viewer
    """

    # Helper function to create data points
    def create_data_points(values, start=1):
        return [
            {"x": i, "y": values[i - start]} for i in range(start, len(values) + start)
        ]

    # Helper function to create binary feature sites
    def create_binary_sites(binary_values, start=1):
        return [
            {"x": i, "y": i, "description": f"Binding site"}
            for i in range(start, len(binary_values) + start)
            if binary_values[i - start] == 1
        ]

    # Generate data for RNA binding probabilities and binary sites
    features_data = {}

    # Process hybridRNAbind data
    features_data["hybridRNAbind_prob"] = json.dumps(
        create_data_points(rna_binding["hybridRNAbind"]["RNA_propensity"])
    )
    features_data["hybridRNAbind_binary"] = json.dumps(
        create_binary_sites(rna_binding["hybridRNAbind"]["RNA_binary"])
    )

    # Process iDRNA-TF data
    features_data["iDRNA_TF_prob"] = json.dumps(
        create_data_points(rna_binding["iDRNA-TF"]["RNA_propensity"])
    )
    features_data["iDRNA_TF_binary"] = json.dumps(
        create_binary_sites(rna_binding["iDRNA-TF"]["RNA_binary"])
    )

    # Fetch protein features from UniProt
    url = f"https://www.ebi.ac.uk/proteins/api/features/{uniprot}"
    features = requests.get(url).json()

    # Process domain and site features
    domain_list = []
    for feature in features.get("features", []):
        if feature["category"] == "DOMAINS_AND_SITES":
            domain_list.append(
                {
                    "x": int(feature.get("begin", 0)),
                    "y": int(feature.get("end", 0)),
                    "description": feature.get("description", ""),
                }
            )
    features_data["domains"] = json.dumps(domain_list)

    # Process PTM features
    ptm_list = []
    for feature in features.get("features", []):
        if feature["category"] == "PTM":
            ptm_list.append(
                {
                    "x": int(feature.get("begin", 0)),
                    "y": int(feature.get("end", 0)),
                    "description": feature.get("description", ""),
                }
            )
    features_data["ptm"] = json.dumps(ptm_list)

    # Process InterPro features
    features_data["interpro"] = json.dumps(interpro_list)

    # Process MRO information if available
    mro_features = []
    if mro_info:
        for role in mro_info:
            mro_features.append(
                {
                    "x": role.get("start", 0),
                    "y": role.get("end", 0),
                    "description": role.get("description", "Moonlighting role"),
                }
            )
    features_data["mro"] = json.dumps(mro_features)

    # Generate feature viewer configuration
    feature_configs = [
        {
            "data": "hybridRNAbind_prob",
            "name": "HybridRNAbind Prob",
            "className": "rnabinding",
            "color": "#008B8D",
            "type": "line",
            "height": "5",
        },
        {
            "data": "iDRNA_TF_prob",
            "name": "iDRNA-TF Prob",
            "className": "rnabinding",
            "color": "#0066CC",
            "type": "line",
            "height": "5",
        },
        {
            "data": "hybridRNAbind_binary",
            "name": "HybridRNAbind Sites",
            "className": "rnasites",
            "color": "#008000",
            "type": "rect",
        },
        {
            "data": "iDRNA_TF_binary",
            "name": "iDRNA-TF Sites",
            "className": "rnasites",
            "color": "#00AA00",
            "type": "rect",
        },
        {
            "data": "domains",
            "name": "Domains and Sites",
            "className": "domains",
            "color": "#FFA500",
            "type": "rect",
        },
        {
            "data": "ptm",
            "name": "PTM",
            "className": "ptm",
            "color": "#FF0000",
            "type": "rect",
        },
        {
            "data": "interpro",
            "name": "InterPro",
            "className": "interpro",
            "color": "#9932CC",
            "type": "rect",
        },
    ]

    # Add MRO features if available
    if mro_info:
        feature_configs.append(
            {
                "data": "mro",
                "name": "Moonlighting Roles",
                "className": "mro",
                "color": "#FF6347",
                "type": "rect",
            }
        )

    # Generate JavaScript feature configurations
    js_features = []
    for config in feature_configs:
        js_features.append(f"""
          ft.addFeature({{
            data: {
            config["data"].startswith("JSON.parse") and config["data"] or config["data"]
        },
            name: "{config["name"]}",
            className: "{config["className"]}",
            color: "{config["color"]}",
            type: "{config["type"]}"{"," if "height" in config else ""}{
            f' height: "{config["height"]}"' if "height" in config else ""
        }
          }});
        """)

    # Generate the HTML
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

          // Data variables
          {"; ".join([f"var {key} = {value}" for key, value in features_data.items()])}

          // Add features
          {"".join(js_features)}
        </script>
      </body>
    </html>
    """

    return html_code

def add_viewer_label(viewer, label_text, viewer_position):
    # Position the label in the lower left corner
    viewer.addLabel(label_text,
                   {'position': {'x': -15, 'y': -15, 'z': 0},  # Position in lower left
                    'backgroundColor': 'rgba(0,0,0,0.6)',
                    'fontColor': 'white',
                    'fontSize': 14,
                    'fontOpacity': 1.0,
                    'borderThickness': 0.0,
                    'inFront': True,
                    'showBackground': True},
                   viewer=viewer_position)
def main():
    # Title of the app
    st.title("Protein Feature Viewer")
    import pandas as pd

    # Definisci i nomi delle colonne
    column_names = [
        "UniProtKB-AC",
        "MD5",
        "Length",
        "Analysis",
        "Signature Accession",
        "Signature Description",
        "Start",
        "End",
        "E-Value",
        "Status",
        "Date",
        "InterPro Accession",
        "InterPro Description",
        "GO Terms",
        "Pathways",
    ]

    # Leggi il file TSV utilizzando i nomi delle colonne
    df = pd.read_csv(
        "data/interpro_humanMPs.tsv", sep="\t", names=column_names, header=None
    )


    # Two tabs for Primary Sequence and Protein Structure
    hybridRNAbind = parse_hybridRNAbind("data/resultshybridRNAbind.txt")
    iDRNA_TF = parse_iDRNA_TF("data/resultsiDRNA-TF.txt")
    # Create selectbox with protein IDs
    hybridRNAbind_ids = list(hybridRNAbind.keys())

    moonlight_proteins = st.session_state["gene_lists"][
        "humanMPs(MoonProt AND MultiTaskProtDB)"
    ]
    protein_ids = list(set(hybridRNAbind_ids).intersection(set(moonlight_proteins)))

    selected_protein = st.selectbox("Select protein:", protein_ids)
    if selected_protein:

        hybridRNAbind_data = hybridRNAbind[selected_protein]
        iDRNA_TF_data = iDRNA_TF[selected_protein]

        sequence = hybridRNAbind_data["sequence"]
        #make rna_binding diz with keys hybridRNAbind and iDRNA-TF and values hybridRNAbind_data and iDRNA_TF_data
        rna_binding = {
            "hybridRNAbind": hybridRNAbind_data,
            "iDRNA-TF": iDRNA_TF_data,
        }



    tab_seq, tab_structure, tab_molstar = st.tabs(["Primary Sequence", "Protein Structure", "Molstar Viewer"])

    with tab_seq:
        st.subheader("Primary Sequence Feature Viewer")



        if selected_protein:

            #extract rows from df with selected protein as UniProtKB-AC
            interpro_rows=df[df["UniProtKB-AC"]==selected_protein]
            interpro_list=[]
            for _, row in interpro_rows.iterrows():
                interpro_list.append(
                    {
                        "x": int(row["Start"]),
                        "y": int(row["End"]),
                        "description": f"{row['Signature Description']} ({row['Signature Accession']}, {row['Analysis']})",
                    }
                )
            row=st.session_state["df"][st.session_state["df"]["UniProtKB-AC"]==selected_protein]
            # Display protein information
            st.write(f"**Protein ID:** {selected_protein}")
            #st.write(f"**UniProtKB Name:** {row['UniProtKB  name'].values[0]}")
            st.write(f"**Gene Name:** {row['Entrez'].values[0]}")
            st.write(f"**Protein names (long):** {row['Protein names'].values[0]}")

            with st.expander("Show GO terms"):
                st.write(f"**GO Biological Process:** {row['Gene Ontology (biological process)'].values[0]}")
                st.write(f"**GO Cellular Component:** {row['Gene Ontology (cellular component)'].values[0]}")
                st.write(f"**GO Molecular Function:** {row['Gene Ontology (molecular function)'].values[0]}")


            # Create and display the feature viewer
            html_code = create_feature_viewer2(selected_protein,sequence,rna_binding,interpro_list,max_points=1515)
            components.html(html_code, height=500, scrolling=True)
            st.download_button(label="Download HTML", data=html_code, file_name="feature_viewer.html", mime="text/html")



    with tab_structure:
        st.subheader("3D Protein Structure Viewer")

        # Create columns for controls and viewer
        col1, col2 = st.columns([1, 4])

        with col1:
            st.subheader("Visualization Controls")

            # Style selector - will trigger updates automatically when changed
            style = st.selectbox("Select Style", ["stick", "sphere", "cartoon"])

            # Color scheme options
            color_options = ["RNA binding HybridRNAbind", "RNA binding iDRNA-TF", "RNA bindtreamling PST-PRNA", "spectrum"]
            color_scheme = st.selectbox("Color Scheme", color_options)

            # Background color
            bg_color = st.color_picker("Background Color", "#FFFFFF")

            # Surface options
            surface_opacity = st.slider("Surface Opacity", 0.0, 1.0, 0.0)
            surface_color = st.color_picker("Surface Color", "#FFFFFF")

        with col2:
            # Fetch data once - this could be cached to improve performance
            @st.cache_data  # Use cache to avoid redundant API calls
            def get_protein_data(protein_name):
                return fetch_structure(protein_name)

            pdb_data = get_protein_data(selected_protein)

            # Initialize viewer
            viewer = py3Dmol.view(data=pdb_data,width=800, height=600)

            # Handle RNA Propensity coloring
            if color_scheme.startswith("RNA binding"):
                # Map RNA propensity scores to the structure
                if (color_scheme == "RNA binding HybridRNAbind"):

                    rna_scores = hybridRNAbind[selected_protein]["RNA_propensity"]
                elif (color_scheme == "RNA binding iDRNA-TF"):
                    rna_scores = iDRNA_TF[selected_protein]["RNA_propensity"]
                elif (color_scheme == "RNA binding PST-PRNA"):
                    # Load PST-PRNA data
                    pstprna_file = (
                        f"data/pstprna/AF-{selected_protein}-F1-model_v4.pdb_result.txt"
                    )
                    pst_data = parse_pstprna(pstprna_file)
                    rna_scores = pst_data["RNA_propensity"]
                    # Create color mapping function
                def get_color_for_score(score):
                    # Convert score to RGB hex string (white to red)
                    # Score 0 -> #FFFFFF (white)
                    # Score 1 -> #FF0000 (red)
                    r = 255
                    g = b = max(0, 255 - int(score * 255))
                    return f"0x{r:02x}{g:02x}{b:02x}"

                # Apply colors to each residue
                for i, score in enumerate(rna_scores):
                    residue_position = i + 1  # 1-based indexing for PDB
                    color = get_color_for_score(score)

                    # Apply color to this residue
                    viewer.setStyle(
                        {"resi": str(residue_position)},
                        {style: {"color": color}},
                    )
            else:
                # Handle other color schemes
                if color_scheme == "spectrum":
                    style_options = {"color": color_scheme}
                else:
                    style_options = {"colorscheme": color_scheme}
                viewer.setStyle({}, style_options)

            # Set background
            viewer.setBackgroundColor(bg_color)

            # Add surface
            viewer.addSurface(
                py3Dmol.VDW,
                {"opacity": surface_opacity, "color": surface_color},
            )

            # Handle highlighting
            viewer.setHoverable(
                {},
                True,
                """function(atom,viewer,event,container) {
                    if(!atom.label) {
                        atom.label = viewer.addLabel(
                            atom.resn + " " + atom.resi + " (" + atom.atom + ")",
                            {
                                position: atom,
                                backgroundColor: 'mintcream',
                                fontColor:'black',
                                fontSize: 12,
                                padding: 2
                            }
                        );
                    }}""",
                """function(atom,viewer) { 
                    if(atom.label) {
                        viewer.removeLabel(atom.label);
                        delete atom.label;
                    }
                }""",
            )
            viewer.zoomTo()

            # Show the molecule
            showmol(viewer,width=800, height=600)
        pstprna_file = (
            f"data/pstprna/AF-{selected_protein}-F1-model_v4.pdb_result.txt"
        )
        pst_data = parse_pstprna(pstprna_file)

        # Get RNA scores
        hybrid_scores = hybridRNAbind[selected_protein]["RNA_propensity"]
        iDRNA_scores = iDRNA_TF[selected_protein]["RNA_propensity"]
        pst_scores = pst_data["RNA_propensity"]

        def get_hex_color(score):
            red = 255
            gb_val = max(0, 255 - int(score * 255))
            return f"0x{red:02x}{gb_val:02x}{gb_val:02x}"

        # Create 2x2 grid viewer
        grid_viewer = py3Dmol.view(width=1200, height=600, viewergrid=(2, 2))

        # Model 1 - HybridRNAbind
        grid_viewer.addModel(pdb_data, "pdb", viewer=(0, 0))
        for i, s in enumerate(hybrid_scores):
            grid_viewer.setStyle(
                {"resi": str(i + 1)},
                {style: {"color": get_hex_color(s)}},
                viewer=(0, 0),
            )

        # Model 2 - iDRNA-TF
        grid_viewer.addModel(pdb_data, "pdb", viewer=(0, 1))
        for i, s in enumerate(iDRNA_scores):
            grid_viewer.setStyle(
                {"resi": str(i + 1)},
                {style: {"color": get_hex_color(s)}},
                viewer=(0, 1),
            )

        # Model 3 - PST-PRNA
        grid_viewer.addModel(pdb_data, "pdb", viewer=(1, 0))
        for i, s in enumerate(pst_scores):
            grid_viewer.setStyle(
                {"resi": str(i + 1)},
                {style: {"color": get_hex_color(s)}},
                viewer=(1, 0),
            )

        # Model 4 - fallback spectrum
        grid_viewer.addModel(pdb_data, "pdb", viewer=(1, 1))
        grid_viewer.setStyle(
            {},
            {
                style: {
                    "colorscheme": {
                        "prop": "partialCharge",
                        "gradient": "rwb",
                        "min": -0.006,
                        "max": 0.006,
                    }
                }
            },viewer=(1, 1)
        )
        #grid_viewer.setStyle({}, {style: {"color": "spectrum"}}, viewer=(1, 1))

        # Global settings
        for r in range(2):
            for c in range(2):
                grid_viewer.setBackgroundColor(bg_color, viewer=(r, c))
                grid_viewer.addSurface(
                    py3Dmol.VDW,
                    {"opacity": surface_opacity, "color": surface_color},
                    viewer=(r, c),
                )
                grid_viewer.setHoverable(
                            {},
                            True,
                            """function(atom,viewer,event,container) {
                                if(!atom.label) {
                                    atom.label = viewer.addLabel(
                                        atom.resn + " " + atom.resi + " (" + atom.atom + ")",
                                        {
                                            position: atom,
                                            backgroundColor: 'mintcream',
                                            fontColor:'black',
                                            fontSize: 12,
                                            padding: 2
                                        }
                                    );
                                }}""",
                            """function(atom,viewer) { 
                                if(atom.label) {
                                    viewer.removeLabel(atom.label);
                                    delete atom.label;
                                }
                            }""",
                            viewer=(r, c),
                        )
        grid_viewer.zoomTo()
        add_viewer_label(grid_viewer, "HybridRNAbind", (0, 0))
        add_viewer_label(grid_viewer, "iDRNA-TF", (0, 1))
        add_viewer_label(grid_viewer, "PST-PRNA", (1, 0))
        add_viewer_label(grid_viewer, "Partial Charge", (1, 1))
        showmol(grid_viewer, width=1200, height=600)
        # After displaying the grid viewer
        # After displaying the grid viewer
        st.markdown(
            """
            <div style="display: flex; justify-content: flex-start; margin: 10px 0; gap: 20px;">
                <div style="text-align: left;">
                    <div style="background: linear-gradient(to right, blue, white, red); height: 20px; width: 200px;"></div>
                    <div style="display: flex; justify-content: space-between; width: 200px;">
                        <span>-0.006</span>
                        <span>0</span>
                        <span>+0.006</span>
                    </div>
                    <div style="font-size: 12px;">Partial Charge (blue→negative, red→positive)</div>
                </div>
                <div style="text-align: left;">
                    <div style="background: linear-gradient(to right, white, red); height: 20px; width: 200px;"></div>
                    <div style="display: flex; justify-content: space-between; width: 200px;">
                        <span>0</span>
                        <span>0.5</span>
                        <span>1</span>
                    </div>
                    <div style="font-size: 12px;">RNA Propensity (white→low, red→high)</div>
                </div>
            </div>
            """,
            unsafe_allow_html=True,
        )



if __name__ == "__main__":
    initialize_page()
    main()
