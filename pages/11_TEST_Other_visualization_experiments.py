import streamlit as st
import streamlit.components.v1 as components

st.header("iCn3D Viewer via iframe")

# URL for the iCn3D viewer with a specific structure (mmdbid=1tup) and parameters.
url = "https://www.ncbi.nlm.nih.gov/Structure/icn3d/?mmdbid=6IQ6&width=300&height=300&closepopup=1&showcommand=0&shownote=0&mobilemenu=1&showtitle=0"

# Embed the iframe with the specified height.
components.iframe(url, height=500)

st.divider()


st.header("Embedding iCn3D as a JavaScript Widget")

widget_html = """
<!DOCTYPE html>
<html>
  <head>
    <meta charset="UTF-8">
    <!-- Include required CSS -->
    <link rel="stylesheet" href="https://www.ncbi.nlm.nih.gov/Structure/icn3d/lib/jquery-ui.min.css">
    <link rel="stylesheet" href="https://www.ncbi.nlm.nih.gov/Structure/icn3d/icn3d.css">
    <!-- Include required JS libraries -->
    <script src="https://www.ncbi.nlm.nih.gov/Structure/icn3d/lib/jquery.min.js"></script>
    <script src="https://www.ncbi.nlm.nih.gov/Structure/icn3d/lib/jquery-ui.min.js"></script>
    <script src="https://www.ncbi.nlm.nih.gov/Structure/icn3d/lib/threeClass.min.js"></script>
    <script src="https://www.ncbi.nlm.nih.gov/Structure/icn3d/icn3d.min.js"></script>
  </head>
  <body>
    <!-- Div where the iCn3D viewer will be rendered -->
    <div id="icn3dwrap" style="width: 100%; height: 500px;"></div>
    <script type="text/javascript">
      $(document).ready(async function() {
          var cfg = {
              divid: 'icn3dwrap',
              width: '100%',
              height: '100%',
              resize: true,
              rotate: 'right',
              mobilemenu: true,
              showcommand: false,
              showtitle: false
          };
          // Set a structure by its MMDB ID
          cfg['mmdbid'] = '1tup';
          var icn3dui = new icn3d.iCn3DUI(cfg);
          // Render the 3D structure
          await icn3dui.show3DStructure();
      });
    </script>
  </body>
</html>
"""

components.html(widget_html, height=550, scrolling=True)

st.divider()
import streamlit as st
import streamlit.components.v1 as components

st.title("Protein Feature Visualization with Nightingale")

# Paste your full HTML/JS code into the multi-line string.
html_code = """
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>Nightingale Protein Features</title>
  <script type="importmap">
    {
      "imports": {
        "@nightingale-elements/": "https://cdn.jsdelivr.net/npm/@nightingale-elements/"
      }
    }
  </script>
  <style>
    td {
      padding: 5px;
    }
    td:first-child {
      background-color: lightcyan;
      font: 0.8em sans-serif;
      white-space: nowrap;
    }
    td:nth-child(2) {
      background-color: aliceblue;
    }
    tr:nth-child(-n + 3) > td {
      background-color: transparent;
    }
  </style>
  <script type="module">
    import "@nightingale-elements/nightingale-sequence@latest";
    import "@nightingale-elements/nightingale-track@latest";
    import "@nightingale-elements/nightingale-manager@latest";
    import "@nightingale-elements/nightingale-navigation@latest";
    import "@nightingale-elements/nightingale-colored-sequence@latest";
    import "@nightingale-elements/nightingale-linegraph-track@latest";

    const accession = "P05067";

    // Load feature and variation data from the Proteins API
    const featuresData = await (
      await fetch("https://www.ebi.ac.uk/proteins/api/features/" + accession)
    ).json();
    const variationData = await (
      await fetch("https://www.ebi.ac.uk/proteins/api/variation/" + accession)
    ).json();

    customElements.whenDefined("nightingale-sequence").then(() => {
      const seq = document.querySelector("#sequence");
      seq.data = featuresData.sequence;
    });

    customElements.whenDefined("nightingale-colored-sequence").then(() => {
      const coloredSeq = document.querySelector("#colored-sequence");
      coloredSeq.data = featuresData.sequence;
    });

    customElements.whenDefined("nightingale-track").then(() => {
      // Nightingale expects "start" rather than the API's "begin"
      const features = featuresData.features.map((ft) => ({
        ...ft,
        start: ft.start || ft.begin,
      }));

      // Filter the data for each feature type and assign to the relevant track data
      const domain = document.querySelector("#domain");
      domain.data = features.filter(({ type }) => type === "DOMAIN");

      const region = document.querySelector("#region");
      region.data = features.filter(({ type }) => type === "REGION");

      const site = document.querySelector("#site");
      site.data = features.filter(({ type }) => type === "SITE");

      const binding = document.querySelector("#binding");
      binding.data = features.filter(({ type }) => type === "BINDING");

      const chain = document.querySelector("#chain");
      chain.data = features.filter(({ type }) => type === "CHAIN");

      const disulfide = document.querySelector("#disulfide-bond");
      disulfide.data = features.filter(({ type }) => type === "DISULFID");

      const betaStrand = document.querySelector("#beta-strand");
      betaStrand.data = features.filter(({ type }) => type === "STRAND");
    });

    customElements.whenDefined("nightingale-linegraph-track").then(() => {
      // Count the variants by position
      const count = {};
      for (const feature of variationData.features) {
        const position = feature.begin;
        count[position] = (count[position] || 0) + 1;
      }
      const max = Math.max(...Object.values(count));
      const values = Object.entries(count).map(([position, value]) => ({
        position: +position,
        value: +value,
      }));
      // Set the variants line graph track data
      const variants = document.querySelector("#variants");
      variants.data = [
        {
          color: "grey",
          values,
          range: [0, max],
        },
      ];
    });
  </script>
</head>
<body>
  <nightingale-manager>
    <table>
      <tbody>
        <tr>
          <td></td>
          <td>
            <nightingale-navigation
              id="navigation"
              min-width="800"
              height="40"
              length="770"
              display-start="1"
              display-end="770"
              margin-color="white"
            ></nightingale-navigation>
          </td>
        </tr>
        <tr>
          <td></td>
          <td>
            <nightingale-sequence
              id="sequence"
              min-width="800"
              height="40"
              length="770"
              display-start="1"
              display-end="770"
              margin-color="white"
              highlight-event="onmouseover"
            ></nightingale-sequence>
          </td>
        </tr>
        <tr>
          <td></td>
          <td>
            <nightingale-colored-sequence
              id="colored-sequence"
              min-width="800"
              height="15"
              length="770"
              display-start="1"
              display-end="770"
              scale="hydrophobicity-scale"
              margin-color="white"
              highlight-event="onmouseover"
            ></nightingale-colored-sequence>
          </td>
        </tr>
        <tr>
          <td>Domain</td>
          <td>
            <nightingale-track
              id="domain"
              min-width="800"
              height="15"
              length="770"
              display-start="1"
              display-end="770"
              margin-color="aliceblue"
              highlight-event="onmouseover"
            ></nightingale-track>
          </td>
        </tr>
        <tr>
          <td>Region</td>
          <td>
            <nightingale-track
              id="region"
              min-width="800"
              height="15"
              length="770"
              display-start="1"
              display-end="770"
              margin-color="aliceblue"
              highlight-event="onmouseover"
            ></nightingale-track>
          </td>
        </tr>
        <tr>
          <td>Site</td>
          <td>
            <nightingale-track
              id="site"
              min-width="800"
              height="15"
              length="770"
              display-start="1"
              display-end="770"
              margin-color="aliceblue"
              highlight-event="onmouseover"
            ></nightingale-track>
          </td>
        </tr>
        <tr>
          <td>Chain</td>
          <td>
            <nightingale-track
              id="chain"
              layout="non-overlapping"
              min-width="800"
              height="15"
              length="770"
              display-start="1"
              display-end="770"
              margin-color="aliceblue"
              highlight-event="onmouseover"
            ></nightingale-track>
          </td>
        </tr>
        <tr>
          <td>Binding site</td>
          <td>
            <nightingale-track
              id="binding"
              min-width="800"
              height="15"
              length="770"
              display-start="1"
              display-end="770"
              margin-color="aliceblue"
              highlight-event="onmouseover"
            ></nightingale-track>
          </td>
        </tr>
        <tr>
          <td>Disulfide bond</td>
          <td>
            <nightingale-track
              id="disulfide-bond"
              layout="non-overlapping"
              min-width="800"
              height="15"
              length="770"
              display-start="1"
              display-end="770"
              margin-color="aliceblue"
              highlight-event="onmouseover"
            ></nightingale-track>
          </td>
        </tr>
        <tr>
          <td>Beta strand</td>
          <td>
            <nightingale-track
              id="beta-strand"
              min-width="800"
              height="15"
              length="770"
              display-start="1"
              display-end="770"
              margin-color="aliceblue"
              highlight-event="onmouseover"
            ></nightingale-track>
          </td>
        </tr>
        <tr>
          <td>Variants</td>
          <td>
            <nightingale-linegraph-track
              id="variants"
              min-width="800"
              height="15"
              length="770"
              display-start="1"
              display-end="770"
              margin-color="aliceblue"
              highlight-event="onmouseover"
            ></nightingale-linegraph-track>
          </td>
        </tr>
      </tbody>
    </table>
  </nightingale-manager>
</body>
</html>
"""

# Embed the HTML snippet into the Streamlit app.
components.html(html_code, height=900, scrolling=True)

st.divider()

import streamlit as st
import streamlit.components.v1 as components

def nightingale_viewer(accession):
    st.title(
        "Protein Feature Visualization with Nightingale – Click to See Component Names"
    )

    # The HTML snippet below embeds the Nightingale visualization.
    # Each Nightingale component has an inline onclick event that calls showComponentName(name)
    # and displays the component's name in the "component-info" div.

    html_code = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
      <meta charset="UTF-8">
      <title>Nightingale Protein Features with Clickable Names</title>
      <script type="importmap">
        {
          "imports": {
            "@nightingale-elements/": "https://cdn.jsdelivr.net/npm/@nightingale-elements/"
          }
        }
      </script>
      <style>
        /* Table and cell styles */
        td {
          padding: 5px;
          vertical-align: middle;
        }
        td:first-child {
          background-color: lightcyan;
          font: 0.8em sans-serif;
          white-space: nowrap;
        }
        td:nth-child(2) {
          background-color: aliceblue;
        }
        tr:nth-child(-n + 3) > td {
          background-color: transparent;
        }
        /* Tooltip styles from previous example (optional) */
        .tooltip {
          position: relative;
          display: inline-block;
          cursor: help;
        }
        .tooltip .tooltiptext {
          visibility: hidden;
          width: 200px;
          background-color: #555;
          color: #fff;
          text-align: center;
          border-radius: 6px;
          padding: 5px 8px;
          position: absolute;
          z-index: 1;
          bottom: 125%;
          left: 50%;
          margin-left: -100px;
          opacity: 0;
          transition: opacity 0.3s;
          font-size: 0.75em;
        }
        .tooltip:hover .tooltiptext {
          visibility: visible;
          opacity: 1;
        }
        /* Style for the information display area */
        #component-info {
          margin-top: 10px;
          padding: 10px;
          background-color: #f0f8ff;
          border: 1px solid #ccc;
          font-size: 1em;
          font-weight: bold;
        }
      </style>
      <script type="module">
        import "@nightingale-elements/nightingale-sequence@latest";
        import "@nightingale-elements/nightingale-track@latest";
        import "@nightingale-elements/nightingale-manager@latest";
        import "@nightingale-elements/nightingale-navigation@latest";
        import "@nightingale-elements/nightingale-colored-sequence@latest";
        import "@nightingale-elements/nightingale-linegraph-track@latest";
    
        // Helper function to show component names on click.
        window.showComponentName = function(name) {
          document.getElementById('component-info').textContent = "Component: " + name;
        };
    
        const accession = {accession};
    
        // Load feature and variation data from the Proteins API
        const featuresData = await (
          await fetch("https://www.ebi.ac.uk/proteins/api/features/" + accession)
        ).json();
        const variationData = await (
          await fetch("https://www.ebi.ac.uk/proteins/api/variation/" + accession)
        ).json();
    
        customElements.whenDefined("nightingale-sequence").then(() => {
          const seq = document.querySelector("#sequence");
          seq.data = featuresData.sequence;
        });
    
        customElements.whenDefined("nightingale-colored-sequence").then(() => {
          const coloredSeq = document.querySelector("#colored-sequence");
          coloredSeq.data = featuresData.sequence;
        });
    
        customElements.whenDefined("nightingale-track").then(() => {
          // Adjust features: use "start" instead of API's "begin"
          const features = featuresData.features.map((ft) => ({
            ...ft,
            start: ft.start || ft.begin,
          }));
    
          document.querySelector("#domain").data = features.filter(({ type }) => type === "DOMAIN");
          document.querySelector("#region").data = features.filter(({ type }) => type === "REGION");
          document.querySelector("#site").data = features.filter(({ type }) => type === "SITE");
          document.querySelector("#binding").data = features.filter(({ type }) => type === "BINDING");
          document.querySelector("#chain").data = features.filter(({ type }) => type === "CHAIN");
          document.querySelector("#disulfide-bond").data = features.filter(({ type }) => type === "DISULFID");
          document.querySelector("#beta-strand").data = features.filter(({ type }) => type === "STRAND");
        });
    
        customElements.whenDefined("nightingale-linegraph-track").then(() => {
          // Count the variants by position
          const count = {};
          for (const feature of variationData.features) {
            const position = feature.begin;
            count[position] = (count[position] || 0) + 1;
          }
          const max = Math.max(...Object.values(count));
          const values = Object.entries(count).map(([position, value]) => ({
            position: +position,
            value: +value,
          }));
          document.querySelector("#variants").data = [
            {
              color: "grey",
              values,
              range: [0, max],
            },
          ];
        });
      </script>
    </head>
    <body>
      <nightingale-manager>
        <table>
          <tbody>
            <tr>
              <td></td>
              <td>
                <nightingale-navigation
                  id="navigation"
                  min-width="800"
                  height="40"
                  length="770"
                  display-start="1"
                  display-end="770"
                  margin-color="white"
                  onclick="showComponentName('Navigation')"
                ></nightingale-navigation>
              </td>
            </tr>
            <tr>
              <td></td>
              <td>
                <nightingale-sequence
                  id="sequence"
                  min-width="800"
                  height="40"
                  length="770"
                  display-start="1"
                  display-end="770"
                  margin-color="white"
                  highlight-event="onmouseover"
                  onclick="showComponentName('Sequence')"
                ></nightingale-sequence>
              </td>
            </tr>
            <tr>
              <td></td>
              <td>
                <nightingale-colored-sequence
                  id="colored-sequence"
                  min-width="800"
                  height="15"
                  length="770"
                  display-start="1"
                  display-end="770"
                  scale="hydrophobicity-scale"
                  margin-color="white"
                  highlight-event="onmouseover"
                  onclick="showComponentName('Colored Sequence')"
                ></nightingale-colored-sequence>
              </td>
            </tr>
            <tr>
              <td>
                <div class="tooltip">Domain
                  <span class="tooltiptext">Domains are conserved functional segments.</span>
                </div>
              </td>
              <td>
                <nightingale-track
                  id="domain"
                  min-width="800"
                  height="15"
                  length="770"
                  display-start="1"
                  display-end="770"
                  margin-color="aliceblue"
                  highlight-event="onmouseover"
                  onclick="showComponentName('Domain')"
                ></nightingale-track>
              </td>
            </tr>
            <tr>
              <td>
                <div class="tooltip">Region
                  <span class="tooltiptext">Regions indicate broader segments.</span>
                </div>
              </td>
              <td>
                <nightingale-track
                  id="region"
                  min-width="800"
                  height="15"
                  length="770"
                  display-start="1"
                  display-end="770"
                  margin-color="aliceblue"
                  highlight-event="onmouseover"
                  onclick="showComponentName('Region')"
                ></nightingale-track>
              </td>
            </tr>
            <tr>
              <td>
                <div class="tooltip">Site
                  <span class="tooltiptext">Sites mark specific residues.</span>
                </div>
              </td>
              <td>
                <nightingale-track
                  id="site"
                  min-width="800"
                  height="15"
                  length="770"
                  display-start="1"
                  display-end="770"
                  margin-color="aliceblue"
                  highlight-event="onmouseover"
                  onclick="showComponentName('Site')"
                ></nightingale-track>
              </td>
            </tr>
            <tr>
              <td>
                <div class="tooltip">Chain
                  <span class="tooltiptext">Chains indicate continuous segments.</span>
                </div>
              </td>
              <td>
                <nightingale-track
                  id="chain"
                  layout="non-overlapping"
                  min-width="800"
                  height="15"
                  length="770"
                  display-start="1"
                  display-end="770"
                  margin-color="aliceblue"
                  highlight-event="onmouseover"
                  onclick="showComponentName('Chain')"
                ></nightingale-track>
              </td>
            </tr>
            <tr>
              <td>
                <div class="tooltip">Binding site
                  <span class="tooltiptext">Binding sites for ligands.</span>
                </div>
              </td>
              <td>
                <nightingale-track
                  id="binding"
                  min-width="800"
                  height="15"
                  length="770"
                  display-start="1"
                  display-end="770"
                  margin-color="aliceblue"
                  highlight-event="onmouseover"
                  onclick="showComponentName('Binding Site')"
                ></nightingale-track>
              </td>
            </tr>
            <tr>
              <td>
                <div class="tooltip">Disulfide bond
                  <span class="tooltiptext">Stabilizes protein structure.</span>
                </div>
              </td>
              <td>
                <nightingale-track
                  id="disulfide-bond"
                  layout="non-overlapping"
                  min-width="800"
                  height="15"
                  length="770"
                  display-start="1"
                  display-end="770"
                  margin-color="aliceblue"
                  highlight-event="onmouseover"
                  onclick="showComponentName('Disulfide Bond')"
                ></nightingale-track>
              </td>
            </tr>
            <tr>
              <td>
                <div class="tooltip">Beta strand
                  <span class="tooltiptext">Elements of secondary structure.</span>
                </div>
              </td>
              <td>
                <nightingale-track
                  id="beta-strand"
                  min-width="800"
                  height="15"
                  length="770"
                  display-start="1"
                  display-end="770"
                  margin-color="aliceblue"
                  highlight-event="onmouseover"
                  onclick="showComponentName('Beta Strand')"
                ></nightingale-track>
              </td>
            </tr>
            <tr>
              <td>
                <div class="tooltip">Variants
                  <span class="tooltiptext">Sequence alterations or mutations.</span>
                </div>
              </td>
              <td>
                <nightingale-linegraph-track
                  id="variants"
                  min-width="800"
                  height="15"
                  length="770"
                  display-start="1"
                  display-end="770"
                  margin-color="aliceblue"
                  highlight-event="onmouseover"
                  onclick="showComponentName('Variants')"
                ></nightingale-linegraph-track>
              </td>
            </tr>
          </tbody>
        </table>
        <!-- This div displays the name of the clicked component -->
        <div id="component-info">Click a component to see its name.</div>
      </nightingale-manager>
    </body>
    </html>
    """

    # Embed the HTML snippet in the Streamlit app.
    components.html(html_code, height=1000, scrolling=True)
