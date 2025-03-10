def nightingale_viewer(accession):
    st.title(
        "Protein Feature Visualization with Nightingale – Click to See Component Names"
    )

    # The HTML snippet below embeds the Nightingale visualization.
    # Each Nightingale component has an inline onclick event that calls showComponentName(name)
    # and displays the component's name in the "component-info" div.
    html_code = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
          <meta charset="UTF-8">
          <title>Nightingale Protein Features with Clickable Names</title>
          <script type="importmap">
            {{
              "imports": {{
                "@nightingale-elements/": "https://cdn.jsdelivr.net/npm/@nightingale-elements/"
              }}
            }}
          </script>
          <style>
            /* Table and cell styles */
            td {{
              padding: 5px;
              vertical-align: middle;
            }}
            td:first-child {{
              background-color: lightcyan;
              font: 0.8em sans-serif;
              white-space: nowrap;
            }}
            td:nth-child(2) {{
              background-color: aliceblue;
            }}
            tr:nth-child(-n + 3) > td {{
              background-color: transparent;
            }}
            /* Tooltip styles from previous example (optional) */
            .tooltip {{
              position: relative;
              display: inline-block;
              cursor: help;
            }}
            .tooltip .tooltiptext {{
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
            }}
            .tooltip:hover .tooltiptext {{
              visibility: visible;
              opacity: 1;
            }}
            /* Style for the information display area */
            #component-info {{
              margin-top: 10px;
              padding: 10px;
              background-color: #f0f8ff;
              border: 1px solid #ccc;
              font-size: 1em;
              fontlatest";
            import "@nightingale-elements/nightingale-navigation@latest";
            import "@nightingale-elements/nightingale-colored-sequence@latest";
            import "@nightingale-elements/nightingale-linegraph-track@latest";

            // Helper function to show component names on click.
            window.showComponentName = function(name) {{
              document.getElementById('component-info').textContent = "Component: " + name;
            }};

            const accession = "{accession}";

            // Load feature and variation data from the Proteins API
            const featuresData = await (
              await fetch("https://www.ebi.ac.uk/proteins/api/features/" + accession)
            ).json();
            const variationData = await (
              await fetch("https://www.ebi.ac.uk/proteins/api/variation/" + accession)
            ).json();

            customElements.whenDefined("nightingale-sequence").then(() => {{
              const seq = document.querySelector("#sequence");
              seq.data = featuresData.sequence;
            }});

            customElements.whenDefined("nightingale-colored-sequence").then(() => {{
              const coloredSeq = document.querySelector("#colored-sequence");
              coloredSeq.data = featuresData.sequence;
            }});

            customElements.whenDefined("nightingale-track").then(() => {{
              // Adjust features: use "start" instead of API's "begin"
              const features = featuresData.features.map((ft) => ({{
                ...ft,
                start: ft.start || ft.begin,
              }}));

              document.querySelector("#domain").data = features.filter(({type}) => type === "DOMAIN");
              document.querySelector("#region").data = features.filter(({type}) => type === "REGION");
              document.querySelector("#site").data = features.filter(({type}) => type === "SITE");
              document.querySelector("#binding").data = features.filter(({type}) => type === "BINDING");
              document.querySelector("#chain").data = features.filter(({type}) => type === "CHAIN");
              document.querySelector("#disulfide-bond").data = features.filter(({type}) => type === "DISULFID");
              document.querySelector("#beta-strand").data = features.filter(({type}) => type === "STRAND");
            }});

            customElements.whenDefined("nightingale-linegraph-track").then(() => {{
              // Count the variants by position
              const count = {{}};
              for (const feature of variationData.features) {{
                const position = feature.begin;
                count[position] = (count[position] || 0) + 1;
              }}
              const max = Math.max(...Object.values(count));
              const values = Object.entries(count).map(([position, value]) => ({{
                position: +position,
                value: +value,
              }}));
              document.querySelector("#variants").data = [
                {{
                  color: "grey",
                  values,
                  range: [0, max],
                }},
              ];
            }});
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
    #make downlaodable
    #st.download_button(html_code, label="Download HTML", file_name="nightingale_viewer.html", mime="text/html")