# ===================== IMPORTS =====================
import streamlit as st
import py3Dmol
import pandas as pd
import io
import zipfile
import re

from supabase import create_client
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski, QED

# ===================== PAGE SETUP =====================
st.set_page_config(
    page_title="NucLigs Database",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ===================== CSS =====================
st.markdown("""<style>
.block-container { padding-top:3.5rem; padding-bottom:3rem; }
.meta-scroll { max-height:70vh; overflow-y:auto; padding-right:12px; }
.feature-card { background:#fff; border:1px solid #e0e0e0; border-radius:10px; padding:15px; margin-bottom:15px; }
.ref-card { background:#fcfcfc; border-left:4px solid #3498db; padding:15px; margin-bottom:15px; border-radius:6px; }
.id-card { background:#f8f9fa; border-left:5px solid #4CAF50; padding:15px; border-radius:8px; margin-bottom:15px; }
.data-row { display:flex; justify-content:space-between; margin-bottom:6px; }
.data-label { font-size:0.85rem; font-weight:600; min-width:140px; }
.data-value { font-family:monospace; font-size:0.9rem; font-weight:600; }
</style>""", unsafe_allow_html=True)

# ===================== SUPABASE =====================
SUPABASE_URL = "https://heuzgnhlrumyfcfigoon.supabase.co"
SUPABASE_KEY = "YOUR_SUPABASE_KEY"

BUCKET_NAME = "NucLigs_PDBs"
METADATA_BUCKET = "codes"
METADATA_FILENAME = "NucLigs_metadata.xlsx"
METADATA_REF_FILENAME = "references.xlsx"

@st.cache_resource
def init_supabase():
    return create_client(SUPABASE_URL, SUPABASE_KEY)

supabase = init_supabase()
LOGO_URL = "https://raw.githubusercontent.com/tushar1298/qwertyui/main/NucLigs.png"

# ===================== LOADERS =====================
@st.cache_data
def load_metadata():
    data = supabase.storage.from_(METADATA_BUCKET).download(METADATA_FILENAME)
    df = pd.read_excel(io.BytesIO(data))
    df.columns = df.columns.str.lower()
    return df

@st.cache_data
def load_references():
    data = supabase.storage.from_(METADATA_BUCKET).download(METADATA_REF_FILENAME)
    df = pd.read_excel(io.BytesIO(data))
    df.columns = df.columns.str.lower()
    return df

# ===================== HELPERS =====================
def format_pubmed(val):
    if pd.isna(val): return ""
    m = re.search(r"\d+", str(val))
    return m.group(0) if m else ""

def linkify(label, value):
    if not value or str(value).lower() == "nan":
        return value

    v = str(value).strip()
    if v.startswith("CHEMBL"):
        return f"<a href='https://www.ebi.ac.uk/chembl/compound_report_card/{v}' target='_blank'>{v}</a>"
    if v.startswith("DB"):
        return f"<a href='https://go.drugbank.com/drugs/{v}' target='_blank'>{v}</a>"
    if v.isdigit():
        return f"<a href='https://pubchem.ncbi.nlm.nih.gov/compound/{v}' target='_blank'>{v}</a>"
    return v

def render_row(label, value):
    st.markdown(
        f"<div class='data-row'><span class='data-label'>{label}</span>"
        f"<span class='data-value'>{value}</span></div>",
        unsafe_allow_html=True
    )

# ===================== HOME =====================
def render_homepage():
    st.markdown(
        f"<div style='text-align:center'>"
        f"<img src='{LOGO_URL}' width='160'><h1>NucLigs Database</h1>"
        f"</div>", unsafe_allow_html=True
    )
    if st.button("Explore the Collection", use_container_width=True):
        st.session_state.page = "database"

# ===================== DATABASE =====================
def render_database():
    metadata = load_metadata()
    refs = load_references()

    selected = st.sidebar.selectbox("Select Structure", metadata["nl"].astype(str))
    row = metadata[metadata["nl"].astype(str) == selected].iloc[0]

    st.subheader(f"Metadata: {selected}")

    # ---------- Metadata ----------
    with st.expander("General Information", expanded=True):
        for k, v in row.items():
            if str(v).lower() != "nan":
                render_row(k.replace("_", " ").title(), linkify(k, v))

    # ---------- References ----------
    st.subheader("References")
    matches = refs[
        (refs.get("nl","").astype(str) == str(selected)) |
        (refs.get("pdbs","").astype(str) == str(row.get("pdbs",""))) |
        (refs.get("names","").astype(str).str.lower() == str(row.get("names","")).lower())
    ]

    if matches.empty:
        st.info("No references found.")
    else:
        for _, r in matches.iterrows():
            st.markdown("<div class='ref-card'>", unsafe_allow_html=True)
            st.markdown(f"<b>{r.get('title','Untitled')}</b>", unsafe_allow_html=True)

            if pd.notna(r.get("doi")):
                doi = r["doi"]
                st.markdown(f"[DOI](https://doi.org/{doi})")

            pm = format_pubmed(r.get("pubmed_id"))
            if pm:
                st.markdown(f"[PubMed](https://pubmed.ncbi.nlm.nih.gov/{pm}/)")

            st.markdown("</div>", unsafe_allow_html=True)

# ===================== ROUTER =====================
if "page" not in st.session_state:
    st.session_state.page = "home"

if st.session_state.page == "home":
    render_homepage()
else:
    render_database()
