import streamlit as st
import requests
import py3Dmol
import pandas as pd

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
    /* Reduce top padding */
    .block-container {
        padding-top: 1.2rem;
        padding-bottom: 2rem;
        padding-left: 2.2rem;
        padding-right: 2.2rem;
    }

    /* Card style for panels */
    .card {
        padding: 1.2rem 1.4rem;
        border-radius: 0.75rem;
        border: 1px solid #e5e7eb;
        background-color: #ffffff;
        box-shadow: 0 4px 14px rgba(15, 23, 42, 0.06);
    }

    /* Section titles */
    h1 {
        font-size: 2.0rem !important;
        font-weight: 700 !important;
        margin-bottom: 0.2rem !important;
    }

    .subtitle {
        font-size: 0.95rem;
        color: #6b7280;
        margin-bottom: 1.2rem;
    }

    /* Selectbox label spacing */
    label[data-baseweb="typography"] {
        font-weight: 500;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# ----------------------------------------------------
# GitHub repo settings
# ----------------------------------------------------
GITHUB_API_URL = "https://api.github.com/repos/tushar1298/qwertyui/contents"
GITHUB_RAW_BASE = "https://raw.githubusercontent.com/tushar1298/qwertyui/main/"

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
    view = py3Dmol.view(width=640, height=520)
    view.addModel(pdb_text, "pdb")
    view.setStyle({"stick": {}})
    view.zoomTo()
    html = view._make_html()
    st.components.v1.html(html, height=540)


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
# UI
# ----------------------------------------------------
st.markdown("### Molecular Structure & Metadata Viewer")
st.markdown(
    '<div class="subtitle">'
    "Browse molecular PDB structures stored in your GitHub repository and view the "
    "associated metadata for each entry."
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

# Top control bar
with st.container():
    st.markdown('<div class="card">', unsafe_allow_html=True)
    col_sel, col_info = st.columns([2, 1])

    with col_sel:
        selected_pdb = st.selectbox("Select structure", pdb_files, index=0)

    with col_info:
        st.markdown(
            "**Repository:** `tushar1298/qwertyui`  \n"
            f"**Metadata file:** `{METADATA_URL.split('/')[-1]}`"
        )

    st.markdown("</div>", unsafe_allow_html=True)

st.markdown("")  # small spacer

if selected_pdb:
    pdb_text = fetch_pdb_from_github(selected_pdb)

    if pdb_text:
        left, right = st.columns([2.2, 1])

        # 3D viewer card
        with left:
            st.markdown('<div class="card">', unsafe_allow_html=True)
            st.markdown(f"#### 3D Structure · `{selected_pdb}`")
            show_3d_pdb(pdb_text)
            st.markdown("</div>", unsafe_allow_html=True)

        # Metadata card (VERTICAL VIEW)
        with right:
            st.markdown('<div class="card">', unsafe_allow_html=True)
            st.markdown("#### Metadata (Details)")

            if metadata_df.empty:
                st.info("No metadata file loaded.")
            else:
                row = find_metadata(metadata_df, selected_pdb)

                if row is None or row.empty:
                    st.info("No metadata found for this entry.")
                else:
                    meta_dict = row.iloc[0].to_dict()

                    # Show each field vertically as key : value
                    for key, value in meta_dict.items():
                        pretty_key = str(key).replace("_", " ").title()
                        st.markdown(
                            f"""
                            <div style="margin-bottom:0.5rem;">
                                <span style="font-weight:600; color:#374151;">{pretty_key}:</span><br>
                                <span style="color:#111827; word-wrap:break-word; font-size:0.92rem;">
                                    {value}
                                </span>
                            </div>
                            """,
                            unsafe_allow_html=True,
                        )

            st.markdown("</div>", unsafe_allow_html=True)
