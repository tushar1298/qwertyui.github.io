import streamlit as st
import requests
import py3Dmol
import pandas as pd

# ----------------------------------------------------
# GitHub repo settings
# ----------------------------------------------------
GITHUB_API_URL = "https://api.github.com/repos/tushar1298/qwertyui/contents"
GITHUB_RAW_BASE = "https://raw.githubusercontent.com/tushar1298/qwertyui/main/"

# 👉 YOUR metadata file name here:
METADATA_URL = "https://raw.githubusercontent.com/tushar1298/qwertyui/main/NucLigs_data_2811.xlsx"


# ----------------------------------------------------
# Load PDB list from GitHub
# ----------------------------------------------------
@st.cache_data
def list_pdb_files():
    r = requests.get(GITHUB_API_URL)
    r.raise_for_status()
    files = r.json()
    pdb_files = [f["name"] for f in files 
                 if isinstance(f, dict) and f.get("name", "").endswith(".pdb")]
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
            st.error(f"❌ Could not fetch PDB file: {filename}")
            return None
    except Exception as e:
        st.error(f"Error fetching PDB: {e}")
        return None


# ----------------------------------------------------
# 3D Viewer
# ----------------------------------------------------
def show_3d_pdb(pdb_text: str):
    view = py3Dmol.view(width=550, height=550)
    view.addModel(pdb_text, "pdb")
    view.setStyle({"stick": {}})
    view.zoomTo()
    html = view._make_html()
    st.components.v1.html(html, height=560)


# ----------------------------------------------------
# Load metadata from excel
# ----------------------------------------------------
@st.cache_data
def load_metadata():
    df = pd.read_excel(METADATA_URL)
    df.columns = [c.strip().lower() for c in df.columns]   # normalize column names
    return df


# ----------------------------------------------------
# Match metadata using "pdbs" column
# ----------------------------------------------------
def find_metadata(metadata_df, pdb_filename):
    """
    metadata column: 'pdbs'
    content example: 'A1.pdb' or just 'A1'
    """
    pdb_filename = pdb_filename.strip()
    pdb_root = pdb_filename.replace(".pdb", "").strip()

    col = "pdbs"
    if col not in metadata_df.columns:
        st.error("❌ Column 'pdbs' not found in metadata file!")
        return None

    # Convert everything to lower string
    col_values = metadata_df[col].astype(str).str.lower()
    needle1 = pdb_filename.lower()
    needle2 = pdb_root.lower()

    match = metadata_df[col_values == needle1]  # exact match with filename
    if match.empty:
        match = metadata_df[col_values == needle2]  # match without .pdb

    return match if not match.empty else None


# ----------------------------------------------------
# STREAMLIT UI
# ----------------------------------------------------
st.title("🧬 GitHub PDB 3D Viewer + Metadata Table")
st.write("Select a PDB file from the repo to view its 3D structure and metadata.")

# Load everything
pdb_files = list_pdb_files()
metadata_df = load_metadata()

if not pdb_files:
    st.error("❌ No PDB files found in your GitHub repo.")
    st.stop()

# Select PDB
selected_pdb = st.selectbox("Select a PDB file", pdb_files)

if selected_pdb:
    pdb_text = fetch_pdb_from_github(selected_pdb)

    if pdb_text:
        col1, col2 = st.columns([2, 1])

        # Left side → 3D viewer
        with col1:
            st.subheader(f"3D Structure: {selected_pdb}")
            show_3d_pdb(pdb_text)

        # Right side → Metadata
        with col2:
            st.subheader("📊 Metadata")
            metadata_row = find_metadata(metadata_df, selected_pdb)

            if metadata_row is None:
                st.info("No metadata found for this PDB.")
            else:
                st.dataframe(metadata_row.reset_index(drop=True))
