import streamlit as st
import py3Dmol
import pandas as pd
import io
import zipfile
import re

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
st.markdown("""
<style>
.block-container {padding-top:3.5rem;padding-bottom:3rem;}
.meta-scroll{max-height:70vh;overflow-y:auto;padding-right:12px;}
.feature-card{background:#fff;border:1px solid #e0e0e0;border-radius:10px;padding:15px;margin-bottom:15px}
.ref-card{background:#fcfcfc;border-left:4px solid #3498db;padding:15px;border-radius:6px;margin-bottom:12px}
.data-row{display:flex;justify-content:space-between;margin-bottom:6px}
.data-label{min-width:140px;font-weight:600;color:#666}
.data-value{font-family:monospace;font-weight:600}
.ref-card .data-row{justify-content:flex-start}
.ref-card .data-label{min-width:60px;margin-right:5px}
.id-card{background:#f8f9fa;border-left:5px solid #4CAF50;padding:15px;border-radius:8px;margin-bottom:15px}
</style>
""", unsafe_allow_html=True)

# ----------------------------------------------------
# Supabase Configuration
# ----------------------------------------------------
SUPABASE_URL = "https://heuzgnhlrumyfcfigoon.supabase.co"
SUPABASE_KEY = "sb_secret_UuFsAopmAmHrdvHf6-mGBg_X0QNgMF5"

BUCKET_PDB = "NucLigs_PDBs"
BUCKET_META = "codes"
META_FILE = "NucLigs_metadata.xlsx"
REF_FILE = "references.xlsx"

@st.cache_resource
def init_supabase():
    return create_client(SUPABASE_URL, SUPABASE_KEY)

supabase = init_supabase()

# ----------------------------------------------------
# Link helpers (NEW)
# ----------------------------------------------------
def link_chembl(cid):
    if not cid or str(cid).lower() == "nan":
        return "N/A"
    cid = str(cid).strip().upper().replace("CHEMBL_", "CHEMBL")
    return f"<a href='https://www.ebi.ac.uk/chembl/compound_report_card/{cid}/' target='_blank'>{cid}</a>"

def link_drugbank(did):
    if not did or str(did).lower() == "nan":
        return "N/A"
    did = str(did).strip()
    return f"<a href='https://go.drugbank.com/drugs/{did}' target='_blank'>{did}</a>"

def link_pubchem(pid):
    if not pid or str(pid).lower() == "nan":
        return "N/A"
    pid = str(pid).strip()
    return f"<a href='https://pubchem.ncbi.nlm.nih.gov/compound/{pid}' target='_blank'>{pid}</a>"

# ----------------------------------------------------
# Data loading
# ----------------------------------------------------
@st.cache_data
def load_metadata():
    data = supabase.storage.from_(BUCKET_META).download(META_FILE)
    df = pd.read_excel(io.BytesIO(data))
    df.columns = [c.lower().replace(" ", "_") for c in df.columns]
    return df

@st.cache_data
def load_references():
    data = supabase.storage.from_(BUCKET_META).download(REF_FILE)
    df = pd.read_excel(io.BytesIO(data))
    df.columns = [c.lower().replace(" ", "_") for c in df.columns]
    return df

def fetch_pdb(name):
    try:
        if not name.lower().endswith(".pdb"):
            name += ".pdb"
        return supabase.storage.from_(BUCKET_PDB).download(name).decode()
    except:
        return None

# ----------------------------------------------------
# Utilities
# ----------------------------------------------------
def render_row(label, value):
    st.markdown(
        f"<div class='data-row'><div class='data-label'>{label}</div>"
        f"<div class='data-value'>{value}</div></div>",
        unsafe_allow_html=True
    )

def show_3d_pdb(pdb_text):
    view = py3Dmol.view(width=900, height=650)
    view.addModel(pdb_text, "pdb")
    view.setStyle({"stick": {"colorscheme": "greenCarbon"}})
    view.zoomTo()
    html = view._make_html()
    html = html.replace(
        "</body>",
        "<div style='text-align:center;margin-top:8px;'>"
        "<button onclick='viewer.png()'>📷 Save PNG Snapshot</button>"
        "</div></body>"
    )
    st.components.v1.html(html, height=720)

# ----------------------------------------------------
# MAIN APP
# ----------------------------------------------------
df = load_metadata()
refs = load_references()

st.sidebar.title("NucLigs Database")
nid = st.sidebar.selectbox("Select NucL ID", sorted(df["nl"].astype(str)))

row = df[df["nl"].astype(str) == nid].iloc[0]
pdb_text = fetch_pdb(row["pdbs"])

if not pdb_text:
    st.error("PDB not found")
    st.stop()

col1, col2 = st.columns([1.6, 1])

# LEFT
with col1:
    st.subheader(f"3D Structure: {nid}")
    show_3d_pdb(pdb_text)

# RIGHT
with col2:
    tab1, tab2, tab3 = st.tabs(["Chemical", "Metadata", "References"])

    # ---------------- METADATA TAB ----------------
    with tab2:
        st.markdown("<div class='meta-scroll'>", unsafe_allow_html=True)

        st.markdown(
            f"<div class='id-card'><b>{nid}</b><br>{row.get('names','')}</div>",
            unsafe_allow_html=True
        )

        exclude = ["nl", "names", "pdbs", "smiles", "inchi"]

        for k, v in row.items():
            if k in exclude or str(v).lower() == "nan":
                continue

            label = k.replace("_", " ").title()

            if k in ["chembl_id", "chembl"]:
                render_row(label, link_chembl(v))
            elif k in ["drugbank_id", "drugbank"]:
                render_row(label, link_drugbank(v))
            elif k in ["pubchem_id", "pubchem", "cid"]:
                render_row(label, link_pubchem(v))
            else:
                render_row(label, v)

        st.markdown("</div>", unsafe_allow_html=True)

    # ---------------- REFERENCES TAB ----------------
    with tab3:
        st.markdown("<div class='meta-scroll'>", unsafe_allow_html=True)

        chembl = row.get("chembl_id", "")
        pdb = row.get("pdbs", "")

        matches = refs[
            (refs.get("chembl_id","").astype(str) == str(chembl)) |
            (refs.get("pdbs","").astype(str) == str(pdb))
        ]

        if matches.empty:
            st.info("No references found")
        else:
            for _, r in matches.iterrows():
                st.markdown("<div class='ref-card'>", unsafe_allow_html=True)
                st.markdown(f"<b>{r.get('title','')}</b>", unsafe_allow_html=True)
                render_row("Journal", r.get("journal",""))
                render_row("Year", r.get("year",""))
                render_row("PubMed", link_pubchem(r.get("pubmed_id","")))
                st.markdown("</div>", unsafe_allow_html=True)

        st.markdown("</div>", unsafe_allow_html=True)
