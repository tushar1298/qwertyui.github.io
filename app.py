import streamlit as st
import requests
import py3Dmol
import pandas as pd
import io

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski
from Bio.PDB import PDBParser

# ----------------------------------------------------
# Page setup
# ----------------------------------------------------
st.set_page_config(
    page_title="Molecular Viewer",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ----------------------------------------------------
# Minimal CSS (Just for scrollbars & headers)
# ----------------------------------------------------
st.markdown(
    """
    <style>
    /* Remove top padding to make header flush */
    .block-container {
        padding-top: 2rem;
        padding-bottom: 2rem;
    }

    /* Style for the metadata scroll area */
    .meta-scroll {
        max-height: 75vh;
        overflow-y: auto;
        padding-right: 10px;
    }

    /* Scrollbar styling for Webkit */
    .meta-scroll::-webkit-scrollbar {
        width: 6px;
    }
    .meta-scroll::-webkit-scrollbar-track {
        background: transparent;
    }
    .meta-scroll::-webkit-scrollbar-thumb {
        background-color: #ccc;
        border-radius: 20px;
    }

    /* Subtle divider */
    .sub-divider {
        margin-top: 1rem;
        margin-bottom: 1rem;
        border-top: 1px solid #f0f2f6;
    }
    
    /* Text styling */
    .label-text {
        font-weight: 600;
        font-size: 0.85rem;
        color: #555;
        margin-bottom: 0px;
    }
    .value-text {
        font-size: 0.95rem;
        color: #111;
        margin-bottom: 12px;
        word-wrap: break-word;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# ----------------------------------------------------
# Constants & URLs
# ----------------------------------------------------
GITHUB_API_URL = "https://api.github.com/repos/tushar1298/qwertyui/contents/PDBs"
GITHUB_RAW_BASE = "https://raw.githubusercontent.com/tushar1298/qwertyui/main/PDBs/"
METADATA_URL = "https://raw.githubusercontent.com/tushar1298/qwertyui/main/NucLigs_data_2811.xlsx"

# ----------------------------------------------------
# Data Fetching Functions
# ----------------------------------------------------
@st.cache_data
def list_pdb_files():
    try:
        r = requests.get(GITHUB_API_URL)
        r.raise_for_status()
        files = r.json()
        pdb_files = [
            f["name"]
            for f in files
            if isinstance(f, dict) and f.get("name", "").lower().endswith(".pdb")
        ]
        return sorted(pdb_files)
    except Exception as e:
        st.error(f"Error fetching file list: {e}")
        return []

def fetch_pdb_from_github(filename: str) -> str | None:
    try:
        url = f"{GITHUB_RAW_BASE}{filename}"
        r = requests.get(url)
        if r.status_code == 200 and r.text.strip():
            return r.text
        return None
    except:
        return None

@st.cache_data
def load_metadata():
    try:
        df = pd.read_excel(METADATA_URL)
        df.columns = [c.strip().lower() for c in df.columns]
        return df
    except Exception:
        return pd.DataFrame()

# ----------------------------------------------------
# Computation Functions
# ----------------------------------------------------
def compute_physchem_from_pdb(pdb_text: str) -> dict:
    props = {}
    try:
        mol = Chem.MolFromPDBBlock(pdb_text, sanitize=True, removeHs=False)
        if mol is None:
            return props

        # Use shorter keys for the metric display
        props["Mol Wt"] = f"{Descriptors.MolWt(mol):.2f}"
        props["LogP"] = f"{Crippen.MolLogP(mol):.2f}"
        props["TPSA"] = f"{rdMolDescriptors.CalcTPSA(mol):.2f}"
        props["H-Acc"] = Lipinski.NumHAcceptors(mol)
        props["H-Don"] = Lipinski.NumHDonors(mol)
        props["Rot. Bonds"] = Lipinski.NumRotatableBonds(mol)
        props["Rings"] = Lipinski.RingCount(mol)
        props["Atoms"] = mol.GetNumAtoms()
    except Exception:
        pass
    return props

def compute_biopython_features(pdb_text: str) -> dict:
    feats = {}
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("lig", io.StringIO(pdb_text))
        # Just simple counts
        feats["Bio Res"] = sum(1 for _ in structure.get_residues())
        feats["Bio Atoms"] = sum(1 for _ in structure.get_atoms())
    except Exception:
        pass
    return feats

def find_metadata(metadata_df, pdb_filename):
    if metadata_df.empty:
        return None
        
    pdb_root = pdb_filename.replace(".pdb", "").lower()
    if "pdbs" not in metadata_df.columns:
        return None
    
    metadata_df["match"] = metadata_df["pdbs"].astype(str).str.lower()
    
    # Try exact match
    match = metadata_df[metadata_df["match"] == pdb_filename.lower()]
    if not match.empty:
        return match

    # Try root match
    match = metadata_df[metadata_df["match"] == pdb_root]
    return match if not match.empty else None

# ----------------------------------------------------
# 3D Viewer Function
# ----------------------------------------------------
def show_3d_pdb(pdb_text: str):
    # Increased height for better immersion
    view = py3Dmol.view(width="100%", height=600)
    view.addModel(pdb_text, "pdb")
    view.setStyle({"stick": {"colorscheme": "greenCarbon"}})
    view.zoomTo()
    view.setBackgroundColor("white") # Matches streamlit light theme usually
    html = view._make_html()
    st.components.v1.html(html, height=600)

# ----------------------------------------------------
# APP LOGIC
# ----------------------------------------------------

# 1. SIDEBAR CONTROLS
with st.sidebar:
    st.markdown("### 🧬 Molecule Explorer")
    
    # Load Data
    all_pdb_files = list_pdb_files()
    metadata_df = load_metadata()

    # Search & Select
    search_query = st.text_input("🔍 Search Structure", placeholder="e.g., 1A2B...")
    
    if search_query:
        pdb_files = [p for p in all_pdb_files if search_query.lower() in p.lower()]
        if not pdb_files:
            st.warning("No matches found.")
            pdb_files = all_pdb_files
    else:
        pdb_files = all_pdb_files

    selected_pdb = st.selectbox("Select PDB File", pdb_files, index=0 if pdb_files else None)
    
    st.markdown("---")
    st.caption(f"**Source:** tushar1298/qwertyui")
    st.caption(f"**Total Files:** {len(all_pdb_files)}")
    if selected_pdb:
         st.caption(f"**Viewing:** {selected_pdb}")

# 2. MAIN AREA
if not selected_pdb:
    st.info("Please select a structure from the sidebar to begin.")
else:
    pdb_text = fetch_pdb_from_github(selected_pdb)

    if pdb_text:
        # Compute properties
        physchem = compute_physchem_from_pdb(pdb_text)
        biopy = compute_biopython_features(pdb_text)
        
        # Main Layout: Viewer (Left/Top) vs Data (Right/Bottom)
        col_viewer, col_data = st.columns([1.8, 1])

        with col_viewer:
            st.subheader(f"Structure Visualization")
            show_3d_pdb(pdb_text)

        with col_data:
            st.subheader("Analysis & Metadata")
            
            # Create a scrolling container for the data column
            # We use a custom HTML div wrapper defined in CSS
            st.markdown('<div class="meta-scroll">', unsafe_allow_html=True)
            
            # --- SECTION 1: COMPUTED METRICS ---
            st.markdown("##### ⚗️ Physico-Chemical Properties")
            
            if physchem:
                # Create a grid of metrics
                c1, c2, c3, c4 = st.columns(4)
                c1.metric("Mol Wt", physchem.get("Mol Wt", "-"))
                c2.metric("LogP", physchem.get("LogP", "-"))
                c3.metric("TPSA", physchem.get("TPSA", "-"))
                c4.metric("Rings", physchem.get("Rings", "-"))

                c5, c6, c7, c8 = st.columns(4)
                c5.metric("H-Acc", physchem.get("H-Acc", "-"))
                c6.metric("H-Don", physchem.get("H-Don", "-"))
                c7.metric("Rot Bonds", physchem.get("Rot. Bonds", "-"))
                c8.metric("Atoms", physchem.get("Atoms", "-"))
            else:
                st.warning("Could not compute RDKit properties.")

            st.markdown('<div class="sub-divider"></div>', unsafe_allow_html=True)

            # --- SECTION 2: REPOSITORY METADATA ---
            st.markdown("##### 📄 Descriptive Metadata")
            
            row = find_metadata(metadata_df, selected_pdb)

            if row is not None and not row.empty:
                data = row.iloc[0].to_dict()
                data.pop("match", None)
                
                # Separate long text fields for better formatting
                long_fields = ["names", "smiles", "inchi", "description"]
                
                # Display standard fields first
                for key, value in data.items():
                    if key.lower() not in long_fields:
                        st.markdown(
                            f"""
                            <div class="label-text">{key.replace('_', ' ').title()}</div>
                            <div class="value-text">{value}</div>
                            """, 
                            unsafe_allow_html=True
                        )
                
                # Display long fields (like SMILES) at the bottom
                for key in long_fields:
                    if key in data:
                        st.markdown(f'<div class="label-text">{key.upper()}</div>', unsafe_allow_html=True)
                        st.code(str(data[key]), language="text")
            else:
                st.info("No external metadata found in Excel sheet.")
            
            st.markdown("</div>", unsafe_allow_html=True) # End scroll
