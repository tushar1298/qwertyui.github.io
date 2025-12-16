import streamlit as st
import py3Dmol
import pandas as pd
import io
import zipfile
import re  # <-- NEW: for cleaning PubMed values

from supabase import create_client

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski, QED

# ----------------------------------------------------
# Page setup
# ----------------------------------------------------
st.set_page_config(
    page_title="NucLigs Database",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ----------------------------------------------------
# CSS Styling
# ----------------------------------------------------
st.markdown(
    """
    <style>
    /* Global Container Padding */
    .block-container {
        padding-top: 3.5rem;
        padding-bottom: 3rem;
    }

    /* Scrollable Containers */
    .meta-scroll {
        max-height: 70vh;
        overflow-y: auto;
        padding-right: 12px;
    }
    .meta-scroll::-webkit-scrollbar {
        width: 8px;
    }
    .meta-scroll::-webkit-scrollbar-track {
        background: #f1f1f1;
        border-radius: 4px;
    }
    .meta-scroll::-webkit-scrollbar-thumb {
        background: #c1c1c1;
        border-radius: 4px;
    }
    .meta-scroll::-webkit-scrollbar-thumb:hover {
        background: #a8a8a8;
    }

    /* Section Cards in DB View */
    .feature-card {
        background-color: #ffffff;
        border: 1px solid #e0e0e0;
        border-radius: 10px;
        padding: 15px;
        margin-bottom: 15px;
        box-shadow: 0 2px 5px rgba(0,0,0,0.03);
    }
    .feature-card h5 {
        color: #2c3e50;
        font-size: 0.95rem;
        font-weight: 700;
        margin-bottom: 12px;
        border-bottom: 2px solid #f0f2f6;
        padding-bottom: 8px;
    }
    
    /* Reference Card */
    .ref-card {
        background-color: #fcfcfc;
        border-left: 4px solid #3498db;
        padding: 15px;
        margin-bottom: 15px;
        border-radius: 6px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
    }
    .ref-title {
        font-weight: 700;
        color: #2c3e50;
        margin-bottom: 8px;
        font-size: 1.0rem;
        line-height: 1.4;
    }
    .ref-meta {
        font-size: 0.85rem;
        color: #555;
        margin-bottom: 8px;
    }

    /* ID Highlight Card */
    .id-card {
        background-color: #f8f9fa;
        border-left: 5px solid #4CAF50;
        padding: 15px;
        border-radius: 8px;
        margin-bottom: 15px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.05);
    }
    .id-label {
        font-size: 0.75rem;
        color: #666;
        font-weight: 700;
        text-transform: uppercase;
        letter-spacing: 0.05em;
    }
    .id-value {
        font-size: 1.5rem;
        font-weight: 800;
        color: #2c3e50;
        margin: 4px 0 8px 0;
        line-height: 1.2;
    }
    .id-sub {
        font-size: 1.0rem;
        color: #555;
        font-weight: 500;
        font-style: italic;
    }

    /* Data rows */
    .data-row {
        display: flex;
        justify-content: space-between;
        margin-bottom: 8px;
        align-items: flex-start;
    }
    .data-label {
        font-size: 0.85rem;
        color: #666;
        font-weight: 600;
        min-width: 140px;
        flex-shrink: 0;
        margin-right: 15px;
        padding-top: 2px;
    }
    .data-value {
        font-family: 'Source Code Pro', monospace;
        font-size: 0.9rem;
        color: #222;
        font-weight: 600;
        text-align: right;
        word-break: break-word;
        flex-grow: 1;
    }

    /* Tighter layout for rows inside reference cards (PubMed spacing) */
    .ref-card .data-row {
        justify-content: flex-start;
        align-items: center;
        margin-bottom: 4px;
    }
    .ref-card .data-label {
        min-width: 55px;
        margin-right: 4px;
    }
    .ref-card .data-value {
        text-align: left;
    }

    .reference-text {
        font-size: 0.75rem;
        color: #999;
        margin-left: 5px;
        font-weight: 400;
    }

    /* Homepage Cards */
    .home-card {
        padding: 25px;
        border-radius: 12px;
        border: 1px solid #eef2f5;
        background-color: white;
        text-align: left;
        transition: all 0.2s ease;
        box-shadow: 0 4px 6px rgba(0,0,0,0.02);
        height: 100%;
    }
    .home-card:hover {
        transform: translateY(-5px);
        box-shadow: 0 12px 20px rgba(0,0,0,0.08);
        border-color: #4CAF50;
    }
    .home-card h3 {
        color: #2c3e50;
        font-size: 1.2rem;
        margin-bottom: 15px;
        font-weight: 700;
    }
    .home-card p {
        color: #555;
        font-size: 0.95rem;
        line-height: 1.6;
    }
    
    /* Buttons */
    .stButton > button {
        width: 100%;
        border-radius: 8px;
        font-weight: 600;
        padding-top: 0.5rem;
        padding-bottom: 0.5rem;
    }
    
    /* Tabs Styling */
    .stTabs [data-baseweb="tab-list"] {
        gap: 10px;
    }
    .stTabs [data-baseweb="tab"] {
        height: 50px;
        white-space: pre-wrap;
        background-color: #f8f9fa;
        border-radius: 8px;
        color: #495057;
        font-weight: 600;
        padding: 0 20px;
    }
    .stTabs [aria-selected="true"] {
        background-color: #e8f5e9;
        color: #2e7d32;
        border-bottom: 2px solid #2e7d32;
    }
    
    /* Badge Styles */
    .badge-pass { background-color: #d4edda; color: #155724; padding: 4px 10px; border-radius: 12px; font-size: 0.75rem; font-weight: 700; }
    .badge-fail { background-color: #f8d7da; color: #721c24; padding: 4px 10px; border-radius: 12px; font-size: 0.75rem; font-weight: 700; }
    
    /* Sidebar Logo Adjustment */
    [data-testid="stSidebar"] img {
        margin-bottom: 0px; 
        border-radius: 10px;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# ----------------------------------------------------
# Supabase Configuration
# ----------------------------------------------------

NEXT_PUBLIC_SUPABASE_URL="https://heuzgnhlrumyfcfigoon.supabase.co"
NEXT_PUBLIC_SUPABASE_PUBLISHABLE_DEFAULT_KEY="sb_publishable_AM951Hs4gISMnct_hoTOkA_CnjMPj97"


BUCKET_NAME = "NucLigs_PDBs"
METADATA_BUCKET = "codes"
METADATA_FILENAME = "NucLigs_metadata.xlsx"
METADATA_REF_FILENAME = "references.xlsx"

@st.cache_resource
def init_supabase():
    try:
        return create_client(NEXT_PUBLIC_SUPABASE_URL, NEXT_PUBLIC_SUPABASE_PUBLISHABLE_DEFAULT_KEY)
    except Exception:
        return None

supabase = init_supabase()

LOGO_URL = "https://raw.githubusercontent.com/tushar1298/qwertyui/main/NucLigs.png"

# ----------------------------------------------------
# Data & Compute Functions
# ----------------------------------------------------
@st.cache_data(ttl=0)
def load_metadata():
    if not supabase:
        return pd.DataFrame()
    try:
        data_bytes = supabase.storage.from_(METADATA_BUCKET).download(METADATA_FILENAME)
        df = pd.read_excel(io.BytesIO(data_bytes))
        df.columns = [str(c).strip().lower().replace(" ", "_") for c in df.columns]
        return df
    except Exception as e:
        st.sidebar.error(f"Error loading metadata from Supabase: {e}")
        return pd.DataFrame()

@st.cache_data(ttl=0)
def load_references():
    if not supabase:
        return pd.DataFrame()
    try:
        data_bytes = supabase.storage.from_(METADATA_BUCKET).download(METADATA_REF_FILENAME)
        df = pd.read_excel(io.BytesIO(data_bytes), sheet_name=0)
        df.columns = [str(c).strip().lower().replace(" ", "_") for c in df.columns]
        return df
    except Exception:
        return pd.DataFrame()

@st.cache_data(ttl=0)
def get_ids_from_metadata():
    df = load_metadata()
    if not df.empty and 'nl' in df.columns:
        return sorted(df['nl'].dropna().astype(str).unique().tolist())
    return []

def fetch_pdb_from_supabase(filename_or_id: str) -> str | None:
    if not supabase:
        return None
    try:
        filename = filename_or_id if filename_or_id.lower().endswith(".pdb") else f"{filename_or_id}.pdb"
        data_bytes = supabase.storage.from_(BUCKET_NAME).download(filename)
        return data_bytes.decode('utf-8')
    except Exception:
        return None

def calculate_esol(mol, logp, mw, rb, aromatic_rings):
    try:
        num_heavy = mol.GetNumHeavyAtoms()
        num_aromatic = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
        aromatic_prop = num_aromatic / num_heavy if num_heavy > 0 else 0
        esol = 0.16 - (0.63 * logp) - (0.0062 * mw) + (0.066 * rb) - (0.74 * aromatic_prop)
        return esol
    except Exception:
        return None

def compute_physchem(mol) -> dict:
    props = {}
    if mol is None:
        return props
    try:
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        h_acc = Lipinski.NumHAcceptors(mol)
        h_don = Lipinski.NumHDonors(mol)
        rb = Lipinski.NumRotatableBonds(mol)
        formula = rdMolDescriptors.CalcMolFormula(mol)
        charge = Chem.GetFormalCharge(mol)
        chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        qed = QED.qed(mol)
        aromatic_rings = Lipinski.NumAromaticRings(mol)
        esol = calculate_esol(mol, logp, mw, rb, aromatic_rings)
        
        violations = 0
        if mw > 500: violations += 1
        if logp > 5: violations += 1
        if h_don > 5: violations += 1
        if h_acc > 10: violations += 1

        props["Formula"] = formula
        props["Charge"] = str(charge)
        props["Chiral Centers"] = str(chiral_centers)
        props["Mol Wt"] = f"{mw:.2f}"
        props["LogP"] = f"{logp:.2f}"
        props["TPSA"] = f"{rdMolDescriptors.CalcTPSA(mol):.2f}"
        props["QED"] = f"{qed:.3f}"
        props["ESOL (LogS)"] = f"{esol:.2f}" if esol is not None else "N/A"
        props["H-Acc"] = h_acc
        props["H-Don"] = h_don
        props["Rot. Bonds"] = rb
        props["Arom. Rings"] = aromatic_rings
        props["Sat. Rings"] = Lipinski.NumSaturatedRings(mol)
        props["Atoms"] = mol.GetNumAtoms()
        props["F-Csp3"] = f"{rdMolDescriptors.CalcFractionCSP3(mol):.2f}"
        props["Lipinski Violations"] = violations
        props["_RDKitMol"] = mol
    except Exception:
        pass
    return props

def show_3d_pdb(pdb_text: str, style_choice: str = "Stick", bg_color: str = "white"):
    """Render 3D viewer and add 'Save PNG Snapshot' button."""
    view = py3Dmol.view(width=900, height=700)
    view.addModel(pdb_text, "pdb")
    
    if style_choice == "Stick":
        view.setStyle({"stick": {"colorscheme": "greenCarbon"}})
    elif style_choice == "Sphere":
        view.setStyle({"sphere": {"colorscheme": "greenCarbon"}})
    elif style_choice == "Line":
        view.setStyle({"line": {"colorscheme": "greenCarbon"}})
    else:
        view.setStyle({"stick": {"colorscheme": "greenCarbon"}})
    
    view.zoomTo()
    view.setBackgroundColor(bg_color)
    html = view._make_html()

    injected = html.replace(
        "</body>",
        """
        <div style="text-align:center;margin-top:8px;">
          <button onclick="viewer.png()"
                  style="padding:6px 14px;border-radius:6px;border:1px solid #ccc;
                         background:#f8f9fa;cursor:pointer;font-size:13px;">
            📷 Save PNG Snapshot
          </button>
        </div>
        </body>
        """
    )

    st.components.v1.html(injected, height=740)

def render_row(label, value, ref=None, help_text=None):
    ref_html = f'<span class="reference-text">({ref})</span>' if ref else ''
    st.markdown(
        f"""<div class="data-row">
               <span class="data-label">{label}</span>
               <span class="data-value">{value}{ref_html}</span>
           </div>""",
        unsafe_allow_html=True
    )

# ---- NEW: clean PubMed field so you only see one ID ----
def format_pubmed(value: object) -> str:
    """Return first numeric PubMed ID from a messy cell."""
    if value is None:
        return ""
    s = str(value)
    if s.lower() == "nan":
        return ""
    # grab first continuous digit block
    m = re.search(r"\d+", s)
    return m.group(0) if m else s.strip()

# ----------------------------------------------------
# PAGE RENDERERS
# ----------------------------------------------------
def render_homepage():
    st.markdown(
        f"""
        <div style="text-align: center; padding-top: 10px;">
            <img src="{LOGO_URL}" width="180" style="margin-bottom: 15px;">
            <h1 style='color: #2c3e50; margin-bottom: 0;'>NucLigs Database</h1>
            <p style='color: #666; font-size: 1.15rem; font-weight: 300;'>
                The Premier Resource for Nucleoside and Nucleotide analog Structures
            </p>
        </div>
        """,
        unsafe_allow_html=True
    )

    st.markdown("<br>", unsafe_allow_html=True)

    st.markdown("""
    <div style='background-color: #f8f9fa; padding: 30px; border-radius: 12px; border-left: 5px solid #4CAF50; box-shadow: 0 4px 15px rgba(0,0,0,0.05); margin-bottom: 40px;'>
        <h3 style='color: #2c3e50; margin-top: 0;'>About the Database</h3>
        <p style='color: #444; font-size: 1.05rem; line-height: 1.6;'>
            The <b>NucLigs Database</b> is a specialized repository designed to facilitate research in the field of nucleic acid targeting by incorporating 
            <b>nucleoside and nucleotide analogs</b> at one place providing a unified platform for detailed analysis and visualization.
        </p>
        <p style='color: #444; font-size: 1.05rem; line-height: 1.6; margin-bottom: 0;'>
            Whether you are involved in <b>rational drug design</b>, <b>structural biology</b>, or <b>cheminformatics</b>, NucLigs offers robust tools to explore 
            the physico-chemical landscape of nucleic acid interactions, supporting the discovery of next-generation therapeutics.
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    f1, f2, f3 = st.columns(3)
    with f1:
        st.markdown("""
        <div class="home-card">
            <h3>3D Visualization</h3>
            <p>Interactive, high-fidelity rendering of ligand-target complexes using Py3Dmol. Inspect binding modes, molecular surfaces, and structural conformations in real-time directly within your browser.</p>
        </div>
        """, unsafe_allow_html=True)
    with f2:
        st.markdown("""
        <div class="home-card">
            <h3>Chemical Profiling</h3>
            <p>Automated calculation of critical molecular descriptors. Access data on Molecular Weight, LogP, TPSA, and Lipinski's Rule of 5 compliance powered by the RDKit cheminformatics engine.</p>
        </div>
        """, unsafe_allow_html=True)
    with f3:
        st.markdown("""
        <div class="home-card">
            <h3>Data Accessibility</h3>
            <p>Seamlessly retrieve standardized structural data. Export ligands and complexes in industry-standard formats (PDB, SDF, MOL2) to integrate directly with your local modeling workflows.</p>
        </div>
        """, unsafe_allow_html=True)
        
    st.markdown("<br><br>", unsafe_allow_html=True)

    _, btn_col, _ = st.columns([1.5, 1, 1.5])
    with btn_col:
        if st.button("Explore the Collection", type="primary", use_container_width=True):
            st.session_state['page'] = 'database'
            st.rerun()
    
    st.markdown(
        "<div style='text-align: center; margin-top: 50px; color: #aaa; font-size: 0.85rem;'>© 2024 NucLigs Database Project • Version 2.2</div>",
        unsafe_allow_html=True
    )

def render_database():
    metadata_df = load_metadata()
    refs_df = load_references()
    all_nuc_ids = get_ids_from_metadata()
    
    # Map ID -> "Name (ID)"
    id_map = {}
    if not metadata_df.empty and 'nl' in metadata_df.columns:
        tmp = metadata_df.set_index('nl')
        name_col = 'names' if 'names' in tmp.columns else ('name' if 'name' in tmp.columns else None)
        if name_col:
            for nid, row in tmp.iterrows():
                chem_name = str(row[name_col]) if pd.notna(row[name_col]) else "Unknown"
                if len(chem_name) > 50:
                    chem_name = chem_name[:47] + "..."
                id_map[str(nid)] = f"{chem_name} ({nid})"

    # Sidebar
    with st.sidebar:
        if st.button("Back to Home"):
            st.session_state['page'] = 'home'
            st.rerun()
            
        st.markdown("---")
        st.markdown("### Structure Finder")

        search_query = st.text_input(
            "Filter database:",
            placeholder="Search NucL ID or Name...",
            label_visibility="collapsed"
        )
        
        is_searching = False
        if search_query:
            is_searching = True
            q = search_query.lower()
            mask = (
                metadata_df['nl'].astype(str).str.lower().str.contains(q, na=False) |
                metadata_df['names'].astype(str).str.lower().str.contains(q, na=False)
            )
            filtered = metadata_df[mask]['nl'].dropna().unique().tolist()
            nuc_ids = sorted([str(x) for x in filtered])
            
            if not nuc_ids:
                st.warning("No matches found.")
                nuc_ids = []
            else:
                st.success(f"Found {len(nuc_ids)} structures")
        else:
            nuc_ids = all_nuc_ids

        if nuc_ids:
            fmt_func = (lambda x: id_map.get(x, x)) if is_searching else (lambda x: x)
            selected_nuc_id = st.selectbox(
                "Select Structure Result:",
                nuc_ids,
                index=0,
                format_func=fmt_func
            )
        else:
            selected_nuc_id = None
        
        # Bulk actions
        if nuc_ids:
            with st.expander("Bulk Actions", expanded=False):
                st.caption("Download multiple structures & data based on your current search.")
                
                download_mode = st.radio(
                    "Download Scope:",
                    ["Select Specific", "All Search Results"],
                    horizontal=True,
                    label_visibility="collapsed"
                )
                
                if download_mode == "Select Specific":
                    default_sel = nuc_ids[:5] if len(nuc_ids) > 0 else []
                    bulk_selected = st.multiselect(
                        "Select structures:",
                        nuc_ids,
                        default=default_sel,
                        format_func=fmt_func
                    )
                else:
                    bulk_selected = nuc_ids
                    st.info(f"Ready to download all {len(bulk_selected)} matching structures.")
                
                fmt = st.selectbox("Format", [".pdb", ".sdf", ".mol"], index=0)
                include_preds = st.checkbox(
                    "Compute Features (Slower)",
                    value=False,
                    help="Runs RDKit calculation for every selected structure."
                )
                
                if st.button("Generate Download Package", use_container_width=True, disabled=len(bulk_selected) == 0):
                    status = st.status("Processing bulk download...", expanded=True)
                    try:
                        zip_buffer = io.BytesIO()
                        csv_data_list = []
                        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
                            progress_bar = status.progress(0)
                            total = len(bulk_selected)
                            for idx, nid in enumerate(bulk_selected):
                                row = metadata_df[metadata_df['nl'].astype(str) == nid]
                                if row.empty:
                                    continue
                                meta_dict = row.iloc[0].to_dict()
                                pdb_fname = str(meta_dict.get('pdbs', ''))
                                content = fetch_pdb_from_supabase(pdb_fname)
                                if not content:
                                    continue
                                final_content = content
                                final_ext = fmt
                                mol_obj = None
                                if fmt != ".pdb" or include_preds:
                                    smiles = str(meta_dict.get('smiles', ''))
                                    try:
                                        if smiles and smiles.lower() != 'nan':
                                            mol_obj = Chem.MolFromSmiles(smiles)
                                        else:
                                            mol_obj = Chem.MolFromPDBBlock(content, sanitize=True, removeHs=False)
                                    except Exception:
                                        mol_obj = None
                                if fmt in [".sdf", ".mol"] and mol_obj:
                                    final_content = Chem.MolToMolBlock(mol_obj)
                                elif fmt in [".sdf", ".mol"] and not mol_obj:
                                    final_content = content
                                    final_ext = ".pdb"
                                zip_file.writestr(f"{nid}{final_ext}", final_content)
                                if include_preds and mol_obj:
                                    preds = compute_physchem(mol_obj)
                                    preds.pop("_RDKitMol", None)
                                    meta_dict.update(preds)
                                csv_data_list.append(meta_dict)
                                progress_bar.progress((idx + 1) / total)
                        
                        status.write("Compressing files...")
                        zip_buffer.seek(0)
                        csv_buffer = None
                        if csv_data_list:
                            df_bulk = pd.DataFrame(csv_data_list)
                            csv_buffer = df_bulk.to_csv(index=False).encode('utf-8')
                        
                        status.update(label="Ready for Download!", state="complete", expanded=True)
                        st.download_button(
                            label=f"Download {len(bulk_selected)} Structures (.zip)",
                            data=zip_buffer,
                            file_name="nucligs_structures.zip",
                            mime="application/zip",
                            use_container_width=True
                        )
                        if csv_buffer:
                            st.download_button(
                                label="Download Data Table (.csv)",
                                data=csv_buffer,
                                file_name="nucligs_data.csv",
                                mime="text/csv",
                                use_container_width=True
                            )
                    except Exception as e:
                        status.update(label="Error processing download", state="error")
                        st.error(f"An error occurred: {str(e)}")

        st.divider()
        st.markdown("**Viewer Settings**")
        bg_mode = st.radio("Background", ["Light", "Dark"], horizontal=True, label_visibility="collapsed")
        bg_color = "white" if bg_mode == "Light" else "#1e1e1e"
        
        style_mode = st.selectbox("Style", ["Stick", "Sphere", "Line"])
        
        st.divider()
        st.caption(f"**Total Entries:** {len(all_nuc_ids)}")

    # Main content
    if not selected_nuc_id:
        st.info("Please search for or select a structure from the sidebar.")
        return

    row = metadata_df[metadata_df['nl'].astype(str) == selected_nuc_id]
    
    if row.empty:
        st.error(f"Metadata not found for ID: {selected_nuc_id}")
        return
    
    pdb_filename = str(row.iloc[0]['pdbs'])
    smiles_str = str(row.iloc[0].get('smiles', ''))
    data = row.iloc[0].to_dict()
    pdb_text = fetch_pdb_from_supabase(pdb_filename)

    if not pdb_text:
        st.error("PDB structure could not be loaded from Supabase.")
        return

    # Build RDKit mol
    mol = None
    if smiles_str and smiles_str.lower() != 'nan' and smiles_str.strip():
        try:
            mol = Chem.MolFromSmiles(smiles_str)
        except Exception:
            mol = None
    if mol is None:
        try:
            mol = Chem.MolFromPDBBlock(pdb_text, sanitize=True, removeHs=False)
        except Exception:
            mol = None

    physchem = compute_physchem(mol)

    col_left, col_right = st.columns([1.5, 1.0])

    # LEFT: 3D viewer + downloads
    with col_left:
        st.subheader(f"3D Visualization: {selected_nuc_id}")
        show_3d_pdb(pdb_text, style_choice=style_mode, bg_color=bg_color)

        st.markdown("##### Export Data")
        d1, d2, d3, d4 = st.columns([1, 1, 1, 1.2]) 
        with d1:
            st.download_button(
                "Download PDB",
                pdb_text,
                f"{selected_nuc_id}.pdb",
                "chemical/x-pdb",
                use_container_width=True
            )
        
        mol_obj = physchem.get("_RDKitMol")
        if mol_obj:
            sdf_data = Chem.MolToMolBlock(mol_obj)
            with d2:
                st.download_button(
                    "Download SDF",
                    sdf_data,
                    f"{selected_nuc_id}.sdf",
                    "chemical/x-mdl-sdfile",
                    use_container_width=True
                )
            with d3:
                st.download_button(
                    "Download MOL2",
                    sdf_data,
                    f"{selected_nuc_id}.mol",
                    "chemical/x-mdl-molfile",
                    use_container_width=True
                )
        else:
            with d2:
                st.button("SDF Unavail.", disabled=True, use_container_width=True)
            with d3:
                st.button("MOL2 Unavail.", disabled=True, use_container_width=True)
        
        with d4:
            full_data = {**data, **physchem}
            full_data.pop("_RDKitMol", None)
            csv_single = pd.DataFrame([full_data]).to_csv(index=False).encode('utf-8')
            st.download_button(
                "Data Profile (.csv)",
                csv_single,
                f"{selected_nuc_id}_data.csv",
                "text/csv",
                use_container_width=True
            )

    # RIGHT: Tabs
    with col_right:
        tab_analysis, tab_metadata, tab_refs = st.tabs(["Chemical Analysis", "Metadata Record", "References"])
        
        # Tab 1: Chemical Analysis
        with tab_analysis:
            st.markdown('<div class="meta-scroll">', unsafe_allow_html=True)
            if physchem:
                st.markdown('<div class="feature-card"><h5>Identity</h5>', unsafe_allow_html=True)
                render_row("Formula", physchem.get("Formula", "-"))
                render_row("Mol Weight", f"{physchem.get('Mol Wt', '-')} da")
                render_row("Formal Charge", physchem.get("Charge", "0"))
                render_row("Stereocenters", physchem.get("Chiral Centers", "0"))
                st.markdown('</div>', unsafe_allow_html=True)

                violations = physchem.get("Lipinski Violations", 0)
                badge_class = "badge-pass" if violations == 0 else "badge-fail"
                badge_text = "PASS (0 Violations)" if violations == 0 else f"FAIL ({violations} Violations)"
                
                st.markdown(
                    f'<div class="feature-card">'
                    f'<div style="display:flex; justify-content:space-between; align-items:center; '
                    f'margin-bottom:10px; border-bottom: 2px solid #f0f2f6; padding-bottom:8px;">'
                    f'<h5 style="margin:0; border:none; padding:0;">Rule of 5</h5>'
                    f'<span class="{badge_class}">{badge_text}</span>'
                    f'</div>',
                    unsafe_allow_html=True
                )
                render_row("LogP", physchem.get("LogP", "-"), "≤ 5")
                render_row("H-Donors", physchem.get("H-Don", "-"), "≤ 5")
                render_row("H-Acceptors", physchem.get("H-Acc", "-"), "≤ 10")
                render_row("Rot. Bonds", physchem.get("Rot. Bonds", "-"))
                render_row("Mol Weight", f"{physchem.get('Mol Wt', '-')} da", "≤ 500")
                st.markdown('</div>', unsafe_allow_html=True)

                st.markdown('<div class="feature-card"><h5>Druglikeness</h5>', unsafe_allow_html=True)
                render_row("QED Score", physchem.get("QED", "-"))
                render_row("Est. Solubility", physchem.get("ESOL (LogS)", "-"))
                render_row("TPSA", f"{physchem.get('TPSA', '-')} Å²")
                render_row("Fraction Csp3", physchem.get("F-Csp3", "-"))
                st.markdown('</div>', unsafe_allow_html=True)
                
                st.markdown('<div class="feature-card"><h5>Ring Systems</h5>', unsafe_allow_html=True)
                render_row("Aromatic Rings", physchem.get("Arom. Rings", "-"))
                render_row("Saturated Rings", physchem.get("Sat. Rings", "-"))
                st.markdown('</div>', unsafe_allow_html=True)
            else:
                st.warning("Unable to compute chemical properties.")
            st.markdown("</div>", unsafe_allow_html=True)

        # Tab 2: Metadata
        with tab_metadata:
            st.markdown('<div class="meta-scroll">', unsafe_allow_html=True)
            if data:
                nl_id = data.get('nl', 'Unknown')
                chem_name = data.get('names', data.get('name', '')) 
                
                st.markdown(
                    f'<div class="id-card">'
                    f'<div class="id-label">NucLigs Identifier</div>'
                    f'<div class="id-value">{nl_id}</div>'
                    f'<div class="id-sub">{chem_name}</div>'
                    f'</div>',
                    unsafe_allow_html=True
                )

                st.markdown('<div class="feature-card"><h5>General Info</h5>', unsafe_allow_html=True)
                exclude_fields = ['nl', 'names', 'name', 'pdbs', 'match', 'smiles', 'inchi', 'description', 'sequence']
                for key, value in data.items():
                    if key not in exclude_fields and str(value).lower() != 'nan':
                        render_row(key.replace('_', ' ').title(), str(value))
                st.markdown('</div>', unsafe_allow_html=True)

                long_fields = ["smiles", "inchi", "description", "sequence"]
                for key in long_fields:
                    if key in data and str(data[key]).lower() != 'nan':
                        with st.expander(key.upper(), expanded=False):
                            st.code(str(data[key]), language="text")
            else:
                st.info("No metadata found.")
            st.markdown("</div>", unsafe_allow_html=True)

        # Tab 3: References
        with tab_refs:
            st.markdown('<div class="meta-scroll">', unsafe_allow_html=True)

            pdb_id = str(data.get('pdbs', 'Unknown')).strip()

            chembl_id = None
            chembl_keys = ['chembl_id', 'chemblid', 'chembl', 'chembl_id_']
            for k in chembl_keys:
                if k in data and pd.notna(data[k]):
                    val = str(data[k]).strip()
                    if val and val.lower() != 'nan':
                        chembl_id = val
                        break

            def norm_chembl(x: str) -> str:
                s = str(x).strip().upper()
                if not s:
                    return ""
                s = s.replace("CHEMBL_", "CHEMBL").replace("CHEMBL ", "CHEMBL")
                if not s.startswith("CHEMBL"):
                    s = "CHEMBL" + s.lstrip("_- ")
                return s

            chembl_norm = norm_chembl(chembl_id) if chembl_id else ""

            st.markdown(
                f"""
                <div class="id-card">
                    <div class="id-label">Structure (PDB)</div>
                    <div class="id-value">{pdb_id}</div>
                    <div class="id-sub">ChEMBL ID: {chembl_norm or 'N/A'}</div>
                </div>
                """,
                unsafe_allow_html=True
            )

            ref_matches = pd.DataFrame()
            if not refs_df.empty:
                candidates = []
                if 'pdbs' in refs_df.columns and pdb_id:
                    pdb_mask = refs_df['pdbs'].astype(str).str.strip().str.upper() == pdb_id.upper()
                    candidates.append(refs_df[pdb_mask])
                if chembl_norm and 'chembl_id' in refs_df.columns:
                    chembl_mask = refs_df['chembl_id'].astype(str).apply(norm_chembl) == chembl_norm
                    candidates.append(refs_df[chembl_mask])
                if candidates:
                    ref_matches = pd.concat(candidates, ignore_index=True).drop_duplicates()

            if not ref_matches.empty:
                for _, ref_row in ref_matches.iterrows():
                    ref_data = ref_row.to_dict()
                    st.markdown('<div class="ref-card">', unsafe_allow_html=True)

                    title = ref_data.get('title', ref_data.get('reference_title', 'N/A'))
                    journal = ref_data.get('journal', ref_data.get('journal_name', 'N/A'))
                    year = ref_data.get('year', ref_data.get('publication_year', 'N/A'))
                    doi = ref_data.get('doi', ref_data.get('doi_id', None))
                    pubmed = ref_data.get('pubmed_id', ref_data.get('pmid', 'N/A'))

                    clean_title = str(title).strip("(),'\"")
                    if clean_title.lower() == 'nan':
                        clean_title = "Untitled Reference"
                    st.markdown(f"<div class='ref-title'>{clean_title}</div>", unsafe_allow_html=True)

                    j_str = str(journal) if str(journal).lower() != 'nan' else ""
                    y_str = str(year) if str(year).lower() != 'nan' else ""
                    if j_str or y_str:
                        meta_str = f"<b>{j_str}</b> ({y_str})"
                        st.markdown(f"<div class='ref-meta'>{meta_str}</div>", unsafe_allow_html=True)

                    row_pdb = ref_data.get('pdbs', pdb_id)
                    row_chembl = ref_data.get('chembl_id', chembl_norm)
                    st.markdown(
                        f"<div class='ref-meta'>PDB: <b>{row_pdb}</b> &nbsp; | &nbsp; ChEMBL: <b>{row_chembl}</b></div>",
                        unsafe_allow_html=True
                    )

                    cols = st.columns([1, 2])
                    with cols[0]:
                        pm_str = format_pubmed(pubmed)   # <-- cleaned PubMed
                        if pm_str:
                            render_row("PubMed", pm_str)

                    with cols[1]:
                        if doi and str(doi).lower() != 'nan':
                            clean_doi = str(doi).strip("(), ")
                            doi_link = clean_doi if clean_doi.startswith("http") else f"https://doi.org/{clean_doi}"
                            st.markdown(
                                f"<span class='data-label'>DOI:</span> "
                                f"<a href='{doi_link}' target='_blank'>{clean_doi}</a>",
                                unsafe_allow_html=True
                            )

                    st.markdown('</div>', unsafe_allow_html=True)
            else:
                st.info("No references found in references.xlsx for this PDB / ChEMBL ID.")

            st.markdown("</div>", unsafe_allow_html=True)

# ----------------------------------------------------
# MAIN ROUTER
# ----------------------------------------------------
if 'page' not in st.session_state:
    st.session_state['page'] = 'home'

if st.session_state['page'] == 'home':
    render_homepage()
else:
    render_database()
