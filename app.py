import streamlit as st
import requests
import py3Dmol
import pandas as pd
import io

# Supabase imports
# NOTE: You must run 'pip install supabase' to use this
from supabase import create_client, Client

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski, QED
from Bio.PDB import PDBParser

# ----------------------------------------------------
# Page setup
# ----------------------------------------------------
st.set_page_config(
    page_title="NucLigs Database",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ----------------------------------------------------
# CSS Styling (Enhanced)
# ----------------------------------------------------
st.markdown(
    """
    <style>
    /* Global Container Padding */
    .block-container {
        padding-top: 2rem;
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

    /* Section Cards */
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

    /* Data Points */
    .data-row {
        display: flex;
        justify-content: space-between;
        margin-bottom: 8px;
        align-items: center;
    }
    .data-label {
        font-size: 0.85rem;
        color: #666;
        font-weight: 500;
    }
    .data-value {
        font-family: 'Source Code Pro', monospace;
        font-size: 0.9rem;
        color: #222;
        font-weight: 600;
    }
    .reference-text {
        font-size: 0.75rem;
        color: #999;
        margin-left: 5px;
        font-weight: 400;
    }

    /* Status Badges */
    .badge-pass {
        background-color: #d4edda;
        color: #155724;
        padding: 4px 10px;
        border-radius: 12px;
        font-size: 0.75rem;
        font-weight: 700;
    }
    .badge-fail {
        background-color: #f8d7da;
        color: #721c24;
        padding: 4px 10px;
        border-radius: 12px;
        font-size: 0.75rem;
        font-weight: 700;
    }

    /* Sidebar Logo Adjustment */
    [data-testid="stSidebar"] img {
        margin-bottom: 20px;
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
# Derived from your DB host: db.heuzgnhlrumyfcfigoon.supabase.co
SUPABASE_URL = "https://heuzgnhlrumyfcfigoon.supabase.co"

# ⚠️ PLACEHOLDER: Enter your Supabase 'anon' public key here
# You can find this in Supabase Dashboard -> Project Settings -> API
SUPABASE_KEY = "sb_secret_UuFsAopmAmHrdvHf6-mGBg_X0QNgMF5"

BUCKET_NAME = "NucLigs_PDBs"

# Initialize Client
# We use st.cache_resource to initialize the connection once
@st.cache_resource
def init_supabase():
    try:
        return create_client(SUPABASE_URL, SUPABASE_KEY)
    except Exception as e:
        return None

supabase = init_supabase()

# ----------------------------------------------------
# External URLs (Metadata/Logo)
# ----------------------------------------------------
METADATA_URL = "https://raw.githubusercontent.com/tushar1298/qwertyui/main/NucLigs_metadata.xlsx"
LOGO_URL = "https://raw.githubusercontent.com/tushar1298/qwertyui/main/NucLigs.png"

# ----------------------------------------------------
# Data Fetching Functions
# ----------------------------------------------------
@st.cache_data(ttl=600)
def list_pdb_files():
    """List PDB files from Supabase Storage Bucket"""
    if not supabase:
        st.error("Supabase client not initialized. Check API Key.")
        return []
        
    try:
        # Fetch list of files from the bucket
        # Note: Depending on folder structure, you might need to adjust path
        res = supabase.storage.from_(BUCKET_NAME).list()
        
        # 'res' is typically a list of dicts/objects
        files = []
        if res:
            for f in res:
                # Supabase storage list returns objects with 'name'
                name = f.get('name', '')
                if name.lower().endswith(".pdb"):
                    files.append(name)
        return sorted(files)
    except Exception as e:
        st.error(f"Error fetching file list from Supabase: {e}")
        return []

def fetch_pdb_from_supabase(filename: str) -> str | None:
    """Download specific PDB file content from Supabase"""
    if not supabase:
        return None
        
    try:
        # Download returns bytes
        data_bytes = supabase.storage.from_(BUCKET_NAME).download(filename)
        # Decode bytes to string for Py3Dmol and RDKit
        return data_bytes.decode('utf-8')
    except Exception as e:
        st.error(f"Error downloading {filename}: {e}")
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
def calculate_esol(mol, logp, mw, rb, aromatic_rings):
    """
    Estimate solubility (ESOL)
    LogS = 0.16 - 0.63(cLogP) - 0.0062(MW) + 0.066(RB) - 0.74(AromaticProportion)
    """
    try:
        num_heavy = mol.GetNumHeavyAtoms()
        num_aromatic = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
        aromatic_prop = num_aromatic / num_heavy if num_heavy > 0 else 0
        
        # ESOL Formula
        esol = 0.16 - (0.63 * logp) - (0.0062 * mw) + (0.066 * rb) - (0.74 * aromatic_prop)
        return esol
    except:
        return None

def compute_physchem_from_pdb(pdb_text: str) -> dict:
    props = {}
    try:
        mol = Chem.MolFromPDBBlock(pdb_text, sanitize=True, removeHs=False)
        if mol is None:
            return props

        # 1. Basic Descriptors
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        h_acc = Lipinski.NumHAcceptors(mol)
        h_don = Lipinski.NumHDonors(mol)
        rb = Lipinski.NumRotatableBonds(mol)
        
        # 2. Formula & Charge
        formula = rdMolDescriptors.CalcMolFormula(mol)
        charge = Chem.GetFormalCharge(mol)
        
        # 3. Stereochemistry
        chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))

        # 4. Advanced Descriptors
        qed = QED.qed(mol)
        aromatic_rings = Lipinski.NumAromaticRings(mol)
        
        # 5. ESOL Calculation
        esol = calculate_esol(mol, logp, mw, rb, aromatic_rings)

        # 6. Lipinski Violations
        # MW <= 500, LogP <= 5, H-Don <= 5, H-Acc <= 10
        violations = 0
        if mw > 500: violations += 1
        if logp > 5: violations += 1
        if h_don > 5: violations += 1
        if h_acc > 10: violations += 1

        # Store nicely formatted strings
        props["Formula"] = formula
        props["Charge"] = str(charge)
        props["Chiral Centers"] = str(chiral_centers)
        props["Mol Wt"] = f"{mw:.2f}"
        props["LogP"] = f"{logp:.2f}"
        props["TPSA"] = f"{rdMolDescriptors.CalcTPSA(mol):.2f}"
        props["QED"] = f"{qed:.3f}"
        props["ESOL (LogS)"] = f"{esol:.2f}" if esol else "N/A"
        
        props["H-Acc"] = h_acc
        props["H-Don"] = h_don
        props["Rot. Bonds"] = rb
        props["Arom. Rings"] = aromatic_rings
        props["Sat. Rings"] = Lipinski.NumSaturatedRings(mol)
        props["Atoms"] = mol.GetNumAtoms()
        props["F-Csp3"] = f"{rdMolDescriptors.CalcFractionCSP3(mol):.2f}"
        props["Lipinski Violations"] = violations
        
        # Store Raw mol for export
        props["_RDKitMol"] = mol
        
    except Exception:
        pass
    return props

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
def show_3d_pdb(pdb_text: str, bg_color: str = "white"):
    view = py3Dmol.view(width="100%", height=700)
    view.addModel(pdb_text, "pdb")
    view.setStyle({"stick": {"colorscheme": "greenCarbon"}})
    view.zoomTo()
    view.setBackgroundColor(bg_color)
    html = view._make_html()
    st.components.v1.html(html, height=700)

# ----------------------------------------------------
# Helper to render a data row
# ----------------------------------------------------
def render_row(label, value, ref=None, help_text=None):
    tooltip = f'title="{help_text}"' if help_text else ''
    ref_html = f'<span class="reference-text">({ref})</span>' if ref else ''
    st.markdown(
        f"""
        <div class="data-row" {tooltip}>
            <span class="data-label">{label}</span>
            <span class="data-value">{value}{ref_html}</span>
        </div>
        """, 
        unsafe_allow_html=True
    )

# ----------------------------------------------------
# APP LOGIC
# ----------------------------------------------------

# 1. SIDEBAR CONTROLS
with st.sidebar:
    st.image(LOGO_URL, use_container_width=True)
    st.markdown("### NucLigs Database")
    
    # Load Data
    all_pdb_files = list_pdb_files()
    metadata_df = load_metadata()

    # --- ENHANCED SEARCH SECTION ---
    st.markdown("#### 🕵️ Overall Search")
    search_query = st.text_input("Filter database:", placeholder="Enter structure ID...", label_visibility="collapsed")
    
    # Filter Logic
    if search_query:
        pdb_files = [p for p in all_pdb_files if search_query.lower() in p.lower()]
        if not pdb_files:
            st.warning("No matches found.")
            pdb_files = []
        else:
            st.success(f"Found {len(pdb_files)} structures")
    else:
        pdb_files = all_pdb_files

    # List Display (Sidebar)
    if pdb_files:
        selected_pdb = st.selectbox("Select Structure Result:", pdb_files, index=0)
    else:
        selected_pdb = None
    
    st.divider()
    
    # Viewer Settings in Sidebar
    st.markdown("**Viewer Settings**")
    bg_mode = st.radio("Background", ["Light", "Dark"], horizontal=True, label_visibility="collapsed")
    bg_color = "white" if bg_mode == "Light" else "#1e1e1e"

    st.divider()
    
    if supabase is None:
        st.error("⚠️ Supabase connection failed. Check Key.")
    
    st.caption(f"**Total Entries:** {len(all_pdb_files)}")
    st.caption("v1.4.0 • Powered by Supabase & RDKit")

# 2. MAIN AREA
if not selected_pdb:
    st.info("👈 Please search for or select a structure from the sidebar.")
else:
    # Fetch from SUPABASE instead of GitHub
    pdb_text = fetch_pdb_from_supabase(selected_pdb)

    if pdb_text:
        # Compute properties
        physchem = compute_physchem_from_pdb(pdb_text)
        
        # Main Layout: 3 Columns
        col_viewer, col_preds, col_meta = st.columns([2.2, 0.9, 0.9])

        # --- COLUMN 1: 3D VIEWER & DOWNLOADS ---
        with col_viewer:
            st.subheader(f"3D Visualization: {selected_pdb}")
            show_3d_pdb(pdb_text, bg_color)
            
            # --- DOWNLOAD OPTIONS ---
            st.markdown("##### 📥 Export Data")
            d1, d2, d3 = st.columns(3)
            
            # 1. PDB Download
            with d1:
                st.download_button(
                    label="Download .PDB",
                    data=pdb_text,
                    file_name=selected_pdb,
                    mime="chemical/x-pdb",
                    use_container_width=True
                )
            
            # 2. SDF Download (Convert on fly)
            mol_obj = physchem.get("_RDKitMol")
            if mol_obj:
                sdf_data = Chem.MolToMolBlock(mol_obj) # V2000 mol block standard
                with d2:
                    st.download_button(
                        label="Download .SDF",
                        data=sdf_data,
                        file_name=selected_pdb.replace('.pdb', '.sdf'),
                        mime="chemical/x-mdl-sdfile",
                        use_container_width=True
                    )
                with d3:
                    st.download_button(
                        label="Download .MOL",
                        data=sdf_data,
                        file_name=selected_pdb.replace('.pdb', '.mol'),
                        mime="chemical/x-mdl-molfile",
                        use_container_width=True
                    )
            else:
                with d2:
                    st.button("SDF Unavail.", disabled=True, use_container_width=True)
                with d3:
                    st.button("MOL Unavail.", disabled=True, use_container_width=True)


        # --- COLUMN 2: SCIENTIFIC PREDICTIONS ---
        with col_preds:
            st.subheader("Chemical Analysis")
            st.markdown('<div class="meta-scroll">', unsafe_allow_html=True)
            
            if physchem:
                # 1. Identity
                st.markdown('<div class="feature-card">', unsafe_allow_html=True)
                st.markdown("<h5>🧪 Identity</h5>", unsafe_allow_html=True)
                render_row("Formula", physchem.get("Formula", "-"))
                render_row("Mol Weight", f"{physchem.get('Mol Wt', '-')} da")
                render_row("Formal Charge", physchem.get("Charge", "0"))
                render_row("Stereocenters", physchem.get("Chiral Centers", "0"))
                st.markdown('</div>', unsafe_allow_html=True)

                # 2. Lipinski Rules with REFERENCE DATA
                violations = physchem.get("Lipinski Violations", 0)
                badge_class = "badge-pass" if violations == 0 else "badge-fail"
                badge_text = "PASS (0 Violations)" if violations == 0 else f"FAIL ({violations} Violations)"
                
                st.markdown('<div class="feature-card">', unsafe_allow_html=True)
                st.markdown(f"""
                    <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:10px; border-bottom: 2px solid #f0f2f6; padding-bottom:8px;">
                        <h5 style="margin:0; border:none; padding:0;">⚖️ Rule of 5</h5>
                        <span class="{badge_class}">{badge_text}</span>
                    </div>
                """, unsafe_allow_html=True)
                
                # Added References here
                render_row("LogP", physchem.get("LogP", "-"), ref="≤ 5")
                render_row("H-Donors", physchem.get("H-Don", "-"), ref="≤ 5")
                render_row("H-Acceptors", physchem.get("H-Acc", "-"), ref="≤ 10")
                render_row("Rot. Bonds", physchem.get("Rot. Bonds", "-")) # RB often cited as <= 10 but not strictly Ro5
                render_row("Mol Weight", f"{physchem.get('Mol Wt', '-')} da", ref="≤ 500")
                
                st.markdown('</div>', unsafe_allow_html=True)

                # 3. Druglikeness
                st.markdown('<div class="feature-card">', unsafe_allow_html=True)
                st.markdown("<h5>💊 Druglikeness</h5>", unsafe_allow_html=True)
                render_row("QED Score", physchem.get("QED", "-"))
                render_row("Est. Solubility", physchem.get("ESOL (LogS)", "-"))
                render_row("TPSA", f"{physchem.get('TPSA', '-')} Å²")
                render_row("Fraction Csp3", physchem.get("F-Csp3", "-"))
                st.markdown('</div>', unsafe_allow_html=True)
                
                # 4. Ring System
                st.markdown('<div class="feature-card">', unsafe_allow_html=True)
                st.markdown("<h5>💍 Ring Systems</h5>", unsafe_allow_html=True)
                render_row("Aromatic Rings", physchem.get("Arom. Rings", "-"))
                render_row("Saturated Rings", physchem.get("Sat. Rings", "-"))
                st.markdown('</div>', unsafe_allow_html=True)
            
            else:
                st.warning("Unable to compute chemical properties for this structure.")
            st.markdown("</div>", unsafe_allow_html=True)


        # --- COLUMN 3: METADATA ---
        with col_meta:
            st.subheader("Metadata Record")
            st.markdown('<div class="meta-scroll">', unsafe_allow_html=True)
            
            row = find_metadata(metadata_df, selected_pdb)

            if row is not None and not row.empty:
                data = row.iloc[0].to_dict()
                data.pop("match", None)
                
                long_fields = ["names", "smiles", "inchi", "description", "sequence"]
                
                # Standard fields Card
                st.markdown('<div class="feature-card">', unsafe_allow_html=True)
                st.markdown("<h5>📋 General Info</h5>", unsafe_allow_html=True)
                for key, value in data.items():
                    if key.lower() not in long_fields and str(value).lower() != 'nan':
                        render_row(key.replace('_', ' ').title(), str(value))
                st.markdown('</div>', unsafe_allow_html=True)

                # Long fields in Expanders to save space
                for key in long_fields:
                    if key in data and str(data[key]).lower() != 'nan':
                        with st.expander(key.upper(), expanded=False):
                            st.code(str(data[key]), language="text")
            else:
                st.info("No repository metadata found for this ID.")
            
            st.markdown("</div>", unsafe_allow_html=True)
