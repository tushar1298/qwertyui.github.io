import streamlit as st
import requests
import py3Dmol
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski
from Bio.PDB import PDBParser
import io

# ----------------------------------------------------
# Page setup
# ----------------------------------------------------
st.set_page_config(
    page_title="Molecular Structure & Metadata Viewer",
    layout="wide",
)

# Minimal, clean CSS tweaks
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
        margin-bottom: 0.35rem;
        color: #4b5563;
    }

    .meta-item {
        padding: 0.35rem 0 0.25rem 0;
        border-bottom: 1px solid #f3f4f6;
    }

    .meta-label {
        font-weight: 600;
        color: #374151;
        font-size: 0.9rem;
    }

    .meta-value {
        color: #111827;
        font-size: 0.9rem;
        line-height: 1.35rem;
        word-break: break-word;
        white-space: pre-wrap;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# ----------------------------------------------------
# GitHub repo settings
# ----------------------------------------------------
GITHUB_API_URL = "https://api.github.com/repos/tushar1298/qwertyui/contents"
GITHUB_RAW_BASE = "https://raw.githubusercontent.com/tushar1298/qwertyui/main"

# Your metadata file in the repo
METADATA_URL = (
    "https://raw.githubusercontent.com/tushar1298/qwertyui/main/NucLigs_data_2811.xlsx"
)

# ----------------------------------------------------
# Load PDB list from GitHub
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
    try:
        r = requests.get(url)
        if r.status_code == 200 and r.text.strip():
            return r.text
        else:
            st.error(f"Could not fetch PDB file: {filename}")
            return None
    except Exception as e:
        st.error(f"Error fetching PDB: {e}")
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
# Physico-chemical property prediction (RDKit)
# ----------------------------------------------------
def compute_physchem_from_pdb(pdb_text: str) -> dict:
    props = {}
    try:
        mol = Chem.MolFromPDBBlock(pdb_text, sanitize=True, removeHs=False)
        if mol is None:
            return props

        props["MolWt (RDKit)"] = round(Descriptors.MolWt(mol), 3)
        props["ExactMolWt"] = round(Descriptors.ExactMolWt(mol), 3)
        props["LogP (Crippen)"] = round(Crippen.MolLogP(mol), 3)
        props["TPSA"] = round(rdMolDescriptors.CalcTPSA(mol), 3)
        props["H-bond Acceptors"] = Lipinski.NumHAcceptors(mol)
        props["H-bond Donors"] = Lipinski.NumHDonors(mol)
        props["Rotatable Bonds"] = Lipinski.NumRotatableBonds(mol)
        props["Ring Count"] = Lipinski.RingCount(mol)
        props["Total Atoms"] = mol.GetNumAtoms()
        props["Heavy Atoms"] = mol.GetNumHeavyAtoms()
    except Exception:
        # fail silently; just return whatever we have
        pass
    return props

# ----------------------------------------------------
# Extra structural info (Biopython)
# ----------------------------------------------------
def compute_biopython_features(pdb_text: str) -> dict:
    feats = {}
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("lig", io.StringIO(pdb_text))
        num_atoms = sum(1 for _ in structure.get_atoms())
        num_residues = sum(1 for _ in structure.get_residues())
        feats["Bio: Atoms"] = num_atoms
        feats["Bio: Residues"] = num_residues
    except Exception:
        pass
    return feats

# ----------------------------------------------------
# Load metadata from excel
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
    """
    Metadata identifier column: 'pdbs'
    Values can be 'Name.pdb' or 'Name'.
    """
    pdb_filename = pdb_filename.strip()
    pdb_root = pdb_filename.replace(".pdb", "").strip()

    if "pdbs" not in metadata_df.columns:
        st.error("Column 'pdbs' not found in metadata file.")
        return None

    col_values = metadata_df["pdbs"].astype(str).str.lower()
    needle_full = pdb_filename.lower()
    needle_root = pdb_root.lower()

    match = metadata_df[col_values == needle_full]
    if match.empty:
        match = metadata_df[col_values == needle_root]

    return match if not match.empty else None

# ----------------------------------------------------
# UI Header
# ----------------------------------------------------
st.markdown('<div class="title-main">Molecular Structure & Metadata Viewer</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="subtitle">'
    "Browse molecular PDB structures stored in your GitHub repository and view the "
    "associated metadata and predicted physico-chemical properties."
    "</div>",
    unsafe_allow_html=True,
)

# Load data
try:
    pdb_files = list_pdb_files()
except Exception as e:
    st.error(f"Error listing PDB files from GitHub: {e}")
    pdb_files = []

try:
    metadata_df = load_metadata()
except Exception as e:
    st.error(f"Error loading metadata file: {e}")
    metadata_df = pd.DataFrame()

if not pdb_files:
    st.error("No PDB files found in the GitHub repository.")
    st.stop()

# ----------------------------------------------------
# Top control bar
# ----------------------------------------------------
with st.container():
    st.markdown('<div class="card">', unsafe_allow_html=True)
    col_sel, col_info = st.columns([2.3, 1.2])

    with col_sel:
        selected_pdb = st.selectbox("Select structure", pdb_files, index=0)

    with col_info:
        st.markdown(
            "**Repository:** `tushar1298/qwertyui`  \n"
            f"**Metadata file:** `{METADATA_URL.split('/')[-1]}`"
        )

    st.markdown("</div>", unsafe_allow_html=True)

st.markdown("")  # spacer

# ----------------------------------------------------
# Main content
# ----------------------------------------------------
if selected_pdb:
    pdb_text = fetch_pdb_from_github(selected_pdb)

    if pdb_text:
        # compute properties once
        physchem_props = compute_physchem_from_pdb(pdb_text)
        biopy_props = compute_biopython_features(pdb_text)

        left, right = st.columns([2.2, 1])

        # 3D viewer card
        with left:
            st.markdown('<div class="card">', unsafe_allow_html=True)
            st.markdown(f"#### 3D Structure · `{selected_pdb}`")
            show_3d_pdb(pdb_text)
            st.markdown("</div>", unsafe_allow_html=True)

        # Metadata card (organized vertical layout + predicted props)
        with right:
            st.markdown('<div class="card">', unsafe_allow_html=True)
            st.markdown("#### Metadata & Properties")

            if metadata_df.empty:
                st.info("No metadata file loaded.")
            else:
                row = find_metadata(metadata_df, selected_pdb)

                if row is None or row.empty:
                    st.info("No metadata found for this entry.")
                else:
                    meta = row.iloc[0].to_dict()

                    # ---------- ORIGINAL METADATA ----------
                    preferred_order = [
                        "nl",
                        "names",
                        "smiles",
                        "formula",
                        "inchi",
                        "inchikey",
                        "molecule name",
                        "molecule_name",
                    ]

                    norm_map = {k: k.lower().replace("_", " ") for k in meta.keys()}

                    ordered_keys = []
                    for pref in preferred_order:
                        for original, norm in norm_map.items():
                            if norm == pref and original not in ordered_keys:
                                ordered_keys.append(original)

                    for k in sorted(meta.keys()):
                        if k not in ordered_keys:
                            ordered_keys.append(k)

                    # Core fields
                    st.markdown('<div class="meta-section-title">Core fields</div>', unsafe_allow_html=True)
                    core_norms = [p.lower().replace("_", " ") for p in preferred_order]
                    for k in ordered_keys:
                        if norm_map[k] in core_norms:
                            pretty_label = norm_map[k].title()
                            value = meta[k]
                            st.markdown(
                                f"""
                                <div class="meta-item">
                                    <div class="meta-label">{pretty_label}</div>
                                    <div class="meta-value">{value}</div>
                                </div>
                                """,
                                unsafe_allow_html=True,
                            )

                    # Additional metadata
                    remaining = [
                        k for k in ordered_keys
                        if norm_map[k] not in core_norms
                    ]
                    if remaining:
                        st.markdown("<br>", unsafe_allow_html=True)
                        st.markdown('<div class="meta-section-title">Additional fields</div>', unsafe_allow_html=True)
                        for k in remaining:
                            pretty_label = norm_map[k].title()
                            value = meta[k]
                            st.markdown(
                                f"""
                                <div class="meta-item">
                                    <div class="meta-label">{pretty_label}</div>
                                    <div class="meta-value">{value}</div>
                                </div>
                                """,
                                unsafe_allow_html=True,
                            )

            # ---------- PREDICTED PHYSICO-CHEMICAL PROPERTIES ----------
            if physchem_props or biopy_props:
                st.markdown("<br>", unsafe_allow_html=True)
                st.markdown(
                    '<div class="meta-section-title">Predicted physico-chemical properties</div>',
                    unsafe_allow_html=True,
                )

                merged_props = {**physchem_props, **biopy_props}
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

            st.markdown("</div>", unsafe_allow_html=True)
