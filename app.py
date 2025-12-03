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
    page_title="Molecular Structure & Metadata Viewer",
    layout="wide",
)

# ----------------------------------------------------
# CSS styling
# ----------------------------------------------------
st.markdown(
    """
    <style>
    .block-container {
        padding-top: 1.2rem;
        padding-bottom: 1.6rem;
        padding-left: 2.0rem;
        padding-right: 2.0rem;
        max-width: 1300px;
    }

    .card {
        padding: 1.0rem 1.2rem;
        border-radius: 0.8rem;
        border: 1px solid #e5e7eb;
        background-color: #ffffff;
        box-shadow: 0 3px 10px rgba(15, 23, 42, 0.06);
    }

    .title-main {
        font-size: 1.9rem;
        font-weight: 700;
        margin-bottom: 0.1rem;
        letter-spacing: 0.02em;
    }

    .subtitle {
        font-size: 0.95rem;
        color: #6b7280;
        margin-bottom: 0.8rem;
    }

    .meta-section-title {
        font-size: 0.95rem;
        font-weight: 600;
        margin: 0.35rem 0 0.25rem 0;
        color: #4b5563;
    }

    .meta-item {
        padding: 0.25rem 0 0.25rem 0;
        border-bottom: 1px solid #f3f4f6;
    }

    .meta-label {
        font-weight: 600;
        color: #374151;
        font-size: 0.88rem;
    }

    .meta-value {
        color: #111827;
        font-size: 0.88rem;
        line-height: 1.35rem;
        word-break: break-word;
        white-space: pre-wrap;
    }

    .meta-value-long {
        font-size: 0.84rem;
        line-height: 1.3rem;
    }

    /* Scrollable metadata panel */
    .scrollable-meta {
        max-height: 520px;
        overflow-y: auto;
        padding-right: 6px;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# ----------------------------------------------------
# GitHub repo settings (PDBs are inside /PDBs/)
# ----------------------------------------------------
GITHUB_API_URL = "https://api.github.com/repos/tushar1298/qwertyui/contents/PDBs"
GITHUB_RAW_BASE = "https://raw.githubusercontent.com/tushar1298/qwertyui/main/PDBs/"

# Metadata file
METADATA_URL = (
    "https://raw.githubusercontent.com/tushar1298/qwertyui/main/NucLigs_data_2811.xlsx"
)

# ----------------------------------------------------
# Load PDB list from GitHub (/PDBs/)
# ----------------------------------------------------
@st.cache_data
def list_pdb_files():
    r = requests.get(GITHUB_API_URL)
    r.raise_for_status()
    files = r.json()
    pdb_files = [
        f["name"]
        for f in files
        if isinstance(f, dict) and f.get("name", "").lower().endswith(".pdb")
    ]
    return sorted(pdb_files)

# ----------------------------------------------------
# Fetch PDB text from GitHub
# ----------------------------------------------------
def fetch_pdb_from_github(filename: str) -> str | None:
    url = f"{GITHUB_RAW_BASE}{filename}"
    r = requests.get(url)
    if r.status_code == 200 and r.text.strip():
        return r.text
    return None

# ----------------------------------------------------
# 3D Viewer
# ----------------------------------------------------
def show_3d_pdb(pdb_text: str):
    view = py3Dmol.view(width=640, height=480)
    view.addModel(pdb_text, "pdb")
    view.setStyle({"stick": {}})
    view.zoomTo()
    html = view._make_html()
    st.components.v1.html(html, height=500)

# ----------------------------------------------------
# RDKit property prediction
# ----------------------------------------------------
def compute_physchem_from_pdb(pdb_text: str) -> dict:
    props = {}
    try:
        mol = Chem.MolFromPDBBlock(pdb_text, sanitize=True, removeHs=False)
        if mol is None:
            return props

        props["MolWt (RDKit)"] = round(Descriptors.MolWt(mol), 3)
        props["ExactMolWt"] = round(Descriptors.ExactMolWt(mol), 3)
        props["Crippen LogP"] = round(Crippen.MolLogP(mol), 3)
        props["TPSA"] = round(rdMolDescriptors.CalcTPSA(mol), 3)
        props["H-bond Acceptors"] = Lipinski.NumHAcceptors(mol)
        props["H-bond Donors"] = Lipinski.NumHDonors(mol)
        props["Rotatable Bonds"] = Lipinski.NumRotatableBonds(mol)
        props["Ring Count"] = Lipinski.RingCount(mol)
        props["Total Atoms"] = mol.GetNumAtoms()
        props["Heavy Atoms"] = mol.GetNumHeavyAtoms()
    except Exception:
        pass
    return props

# ----------------------------------------------------
# Biopython structural features
# ----------------------------------------------------
def compute_biopython_features(pdb_text: str) -> dict:
    feats = {}
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("lig", io.StringIO(pdb_text))
        feats["Bio: Atoms"] = sum(1 for _ in structure.get_atoms())
        feats["Bio: Residues"] = sum(1 for _ in structure.get_residues())
    except Exception:
        pass
    return feats

# ----------------------------------------------------
# Load Excel metadata
# ----------------------------------------------------
@st.cache_data
def load_metadata():
    df = pd.read_excel(METADATA_URL)
    df.columns = [c.strip().lower() for c in df.columns]
    return df

# ----------------------------------------------------
# Match metadata using "pdbs" column
# ----------------------------------------------------
def find_metadata(metadata_df, pdb_filename):
    pdb_root = pdb_filename.replace(".pdb", "").lower()
    if "pdbs" not in metadata_df.columns:
        return None
    metadata_df["match"] = metadata_df["pdbs"].astype(str).str.lower()

    match = metadata_df[metadata_df["match"] == pdb_filename.lower()]
    if not match.empty:
        return match

    match = metadata_df[metadata_df["match"] == pdb_root]
    return match if not match.empty else None

# ----------------------------------------------------
# UI HEADER
# ----------------------------------------------------
st.markdown('<div class="title-main">Molecular Structure & Metadata Viewer</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="subtitle">'
    "Browse molecular PDB structures stored in GitHub with linked metadata and auto-computed chemical features."
    "</div>",
    unsafe_allow_html=True,
)

# Load data
all_pdb_files = list_pdb_files()
metadata_df = load_metadata()

# ----------------------------------------------------
# Top selection bar with search
# ----------------------------------------------------
with st.container():
    st.markdown('<div class="card">', unsafe_allow_html=True)
    sel_col, info_col = st.columns([2.3, 1.2])

    with sel_col:
        search = st.text_input("Filter structures (type part of ID/name)", "")
        if search:
            pdb_files = [p for p in all_pdb_files if search.lower() in p.lower()]
            if not pdb_files:
                st.warning("No structures match this filter.")
                pdb_files = all_pdb_files
        else:
            pdb_files = all_pdb_files

        selected_pdb = st.selectbox("Select structure", pdb_files)

    with info_col:
        st.markdown(
            "**Repository:** `tushar1298/qwertyui/PDBs`  \n"
            f"**Metadata file:** `{METADATA_URL.split('/')[-1]}`"
        )

    st.markdown("</div>", unsafe_allow_html=True)

st.markdown("")  # small spacer

# ----------------------------------------------------
# MAIN CONTENT
# ----------------------------------------------------
pdb_text = fetch_pdb_from_github(selected_pdb)

if pdb_text:
    physchem = compute_physchem_from_pdb(pdb_text)
    biopy = compute_biopython_features(pdb_text)

    col_left, col_right = st.columns([2.2, 1])

    # ---------------- 3D VIEWER ----------------
    with col_left:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown(f"#### 3D Structure · `{selected_pdb}`")
        show_3d_pdb(pdb_text)
        st.markdown("</div>", unsafe_allow_html=True)

    # ---------------- METADATA PANEL ----------------
    with col_right:
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.markdown("#### Metadata & Properties")

        st.markdown('<div class="scrollable-meta">', unsafe_allow_html=True)

        row = find_metadata(metadata_df, selected_pdb)

        if row is not None and not row.empty:
            data = row.iloc[0].to_dict()
            # remove helper column if present
            data.pop("match", None)

            # show metadata first
            st.markdown('<div class="meta-section-title">Metadata</div>', unsafe_allow_html=True)

            # define which fields are "long text"
            long_fields = {"names", "smiles", "inchi"}

            for key, value in data.items():
                label = key.title()
                cls = "meta-value-long" if key.lower() in long_fields else "meta-value"
                st.markdown(
                    f"""
                    <div class="meta-item">
                        <div class="meta-label">{label}</div>
                        <div class="{cls}">{value}</div>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )
        else:
            st.info("No metadata found for this entry.")

        # Predicted properties
        st.markdown(
            '<br><div class="meta-section-title">Predicted physico-chemical properties</div>',
            unsafe_allow_html=True,
        )
        merged_props = {**physchem, **biopy}

        for key, value in merged_props.items():
            st.markdown(
                f"""
                <div class="meta-item">
                    <div class="meta-label">{key}</div>
                    <div class="meta-value">{value}</div>
                </div>
                """,
                unsafe_allow_html=True,
            )

        st.markdown("</div>", unsafe_allow_html=True)  # close scrollable-meta
        st.markdown("</div>", unsafe_allow_html=True)  # close card
