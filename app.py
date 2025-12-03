import streamlit as st
import pandas as pd
import psycopg2
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import py3Dmol

# ---------------- DB CONNECTION ----------------

@st.cache_resource
def get_connection():
    db = st.secrets["db"]
    conn = psycopg2.connect(
        host=db["host"],
        database=db["dbname"],
        user=db["user"],
        password=db["password"],
        port=db["port"],
    )
    return conn

@st.cache_data
def load_molecules():
    conn = get_connection()
    df = pd.read_sql("SELECT * FROM NucLigs_data;", conn)
    return df

# ---------------- 2D & 3D FUNCTIONS ----------------

def show_2d(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        st.error("Invalid SMILES")
        return
    img = Draw.MolToImage(mol, size=(300, 300))
    st.image(img, caption="2D Structure")

def smiles_to_3d(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    return Chem.MolToMolBlock(mol)

def show_3d(smiles):
    molblock = smiles_to_3d(smiles)
    if molblock is None:
        st.error("Cannot generate 3D structure")
        return
    view = py3Dmol.view(width=400, height=400)
    view.addModel(molblock, "mol")
    view.setStyle({"stick": {}})
    view.zoomTo()
    st.components.v1.html(view._make_html(), height=450)

# ---------------- STREAMLIT UI ----------------

st.title("🧬 Molecule Viewer (2D + 3D)")

df = load_molecules()

st.sidebar.header("Search Molecule")
search = st.sidebar.text_input("Search by NL or name")

if search:
    results = df[df["NL"].str.contains(search, case=False) | 
                 df["names"].str.contains(search, case=False)]
else:
    results = df

st.sidebar.write(f"Total found: {len(results)}")

selected = st.sidebar.selectbox(
    "Choose a molecule",
    results.index,
    format_func=lambda x: f"{results.at[x,'NL']}  -  {results.at[x,'names']}",
)

row = results.loc[selected]

st.subheader(f"**{row['names']}** ({row['NL']})")

# show properties
st.write(row)

# 2D view
st.subheader("2D Structure")
show_2d(row["smiles"])

# 3D view
st.subheader("3D Structure")
show_3d(row["smiles"])
