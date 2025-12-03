import streamlit as st
import requests
import py3Dmol

# GitHub repo settings
GITHUB_API_URL = "https://api.github.com/repos/tushar1298/qwertyui/contents"
GITHUB_RAW_BASE = "https://raw.githubusercontent.com/tushar1298/qwertyui/main/"

@st.cache_data
def list_pdb_files():
    r = requests.get(GITHUB_API_URL)
    files = r.json()
    pdb_files = [f["name"] for f in files if f["name"].endswith(".pdb")]
    return pdb_files

def fetch_pdb_from_github(filename: str) -> str | None:
    url = f"{GITHUB_RAW_BASE}{filename}"
    try:
        r = requests.get(url)
        if r.status_code == 200:
            return r.text
        else:
            st.error("❌ Could not fetch PDB file.")
            return None
    except Exception as e:
        st.error(f"Error: {e}")
        return None

def show_3d_pdb(pdb_text: str):
    view = py3Dmol.view(width=500, height=500)
    view.addModel(pdb_text, "pdb")
    view.setStyle({"stick": {}})
    view.zoomTo()
    html = view._make_html()
    st.components.v1.html(html, height=520)

st.title("🧬 GitHub PDB Viewer")

pdb_files = list_pdb_files()

if not pdb_files:
    st.warning("No PDB files found in GitHub repo.")
else:
    choice = st.selectbox("Select a PDB file", pdb_files)
    if choice:
        pdb_text = fetch_pdb_from_github(choice)
        if pdb_text:
            show_3d_pdb(pdb_text)
