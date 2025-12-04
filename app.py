import streamlit as st
import requests
import py3Dmol
import pandas as pd
import io
import zipfile

# Supabase imports
# NOTE: You must run 'pip install supabase' to use this
from supabase import create_client, Client

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski, QED
from Bio.PDB import PDBParser

# ----------------------------------------------------
# Page setup (Must be first)
# ----------------------------------------------------
st.set_page_config(
    page_title="NucLigs Database",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ----------------------------------------------------
# CSS Styling
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

    /* Section Cards in DB View */
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

    /* ID Highlight Card */
    .id-card {
        background-color: #f8f9fa;
        border-left: 5px solid #4CAF50;
        padding: 15px;
        border-radius: 8px;
        margin-bottom: 15px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.05);
    }
    .id-label {
        font-size: 0.75rem;
        color: #666;
        font-weight: 700;
        text-transform: uppercase;
        letter-spacing: 0.05em;
    }
    .id-value {
        font-size: 1.5rem;
        font-weight: 800;
        color: #2c3e50;
        margin: 4px 0 8px 0;
        line-height: 1.2;
    }
    .id-sub {
        font-size: 1.0rem;
        color: #555;
        font-weight: 500;
        font-style: italic;
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

    /* Homepage Cards */
    .home-card {
        padding: 25px;
        border-radius: 12px;
        border: 1px solid #eef2f5;
        background-color: white;
        text-align: left;
        transition: all 0.2s ease;
        box-shadow: 0 4px 6px rgba(0,0,0,0.02);
        height: 100%;
    }
    .home-card:hover {
        transform: translateY(-5px);
        box-shadow: 0 12px 20px rgba(0,0,0,0.08);
        border-color: #4CAF50;
    }
    .home-card h3 {
        color: #2c3e50;
        font-size: 1.2rem;
        margin-bottom: 15px;
        font-weight: 700;
    }
    .home-card p {
        color: #555;
        font-size: 0.95rem;
        line-height: 1.6;
    }
    
    /* Buttons */
    .stButton > button {
        width: 100%;
        border-radius: 8px;
        font-weight: 600;
        padding-top: 0.5rem;
        padding-bottom: 0.5rem;
    }
    
    /* Badge Styles */
    .badge-pass { background-color: #d4edda; color: #155724; padding: 4px 10px; border-radius: 12px; font-size: 0.75rem; font-weight: 700; }
    .badge-fail { background-color: #f8d7da; color: #721c24; padding: 4px 10px; border-radius: 12px; font-size: 0.75rem; font-weight: 700; }
    
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
SUPABASE_URL = "https://heuzgnhlrumyfcfigoon.supabase.co"
SUPABASE_KEY = "sb_secret_UuFsAopmAmHrdvHf6-mGBg_X0QNgMF5"
BUCKET_NAME = "NucLigs_PDBs"

@st.cache_resource
def init_supabase():
    try:
        return create_client(SUPABASE_URL, SUPABASE_KEY)
    except Exception as e:
        return None

supabase = init_supabase()

# ----------------------------------------------------
# External URLs
# ----------------------------------------------------
METADATA_URL = "https://raw.githubusercontent.com/tushar1298/qwertyui/main/NucLigs_metadata.xlsx"
LOGO_URL = "https://raw.githubusercontent.com/tushar1298/qwertyui/main/NucLigs.png"

# ----------------------------------------------------
# Data & Compute Functions
# ----------------------------------------------------
@st.cache_data
def load_metadata():
    try:
        df = pd.read_excel(METADATA_URL)
        df.columns = [str(c).strip().lower() for c in df.columns]
        return df
    except Exception:
        return pd.DataFrame()

@st.cache_data
def get_ids_from_metadata():
    df = load_metadata()
    if not df.empty and 'nl' in df.columns:
        return sorted(df['nl'].dropna().astype(str).unique().tolist())
    return []

def fetch_pdb_from_supabase(filename_or_id: str) -> str | None:
    if not supabase: return None
    try:
        filename = filename_or_id if filename_or_id.lower().endswith(".pdb") else f"{filename_or_id}.pdb"
        data_bytes = supabase.storage.from_(BUCKET_NAME).download(filename)
        return data_bytes.decode('utf-8')
    except Exception as e:
        return None

def calculate_esol(mol, logp, mw, rb, aromatic_rings):
    try:
        num_heavy = mol.GetNumHeavyAtoms()
        num_aromatic = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
        aromatic_prop = num_aromatic / num_heavy if num_heavy > 0 else 0
        esol = 0.16 - (0.63 * logp) - (0.0062 * mw) + (0.066 * rb) - (0.74 * aromatic_prop)
        return esol
    except:
        return None

def compute_physchem(mol) -> dict:
    props = {}
    if mol is None: return props
    try:
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        h_acc = Lipinski.NumHAcceptors(mol)
        h_don = Lipinski.NumHDonors(mol)
        rb = Lipinski.NumRotatableBonds(mol)
        formula = rdMolDescriptors.CalcMolFormula(mol)
        charge = Chem.GetFormalCharge(mol)
        chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        qed = QED.qed(mol)
        aromatic_rings = Lipinski.NumAromaticRings(mol)
        esol = calculate_esol(mol, logp, mw, rb, aromatic_rings)
        
        violations = 0
        if mw > 500: violations += 1
        if logp > 5: violations += 1
        if h_don > 5: violations += 1
        if h_acc > 10: violations += 1

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
        props["_RDKitMol"] = mol
    except Exception: pass
    return props

def show_3d_pdb(pdb_text: str, bg_color: str = "white"):
    view = py3Dmol.view(width="100%", height=700)
    view.addModel(pdb_text, "pdb")
    view.setStyle({"stick": {"colorscheme": "greenCarbon"}})
    view.zoomTo()
    view.setBackgroundColor(bg_color)
    html = view._make_html()
    st.components.v1.html(html, height=700)

def render_row(label, value, ref=None, help_text=None):
    ref_html = f'<span class="reference-text">({ref})</span>' if ref else ''
    st.markdown(f"""<div class="data-row"><span class="data-label">{label}</span><span class="data-value">{value}{ref_html}</span></div>""", unsafe_allow_html=True)

# ----------------------------------------------------
# PAGE RENDERERS
# ----------------------------------------------------

def render_homepage():
    # Centered Header
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        # Reduced logo size using width parameter
        st.image(LOGO_URL, width=250)
        st.markdown("<h1 style='text-align: center; color: #2c3e50; margin-bottom: 0;'>NucLigs Database</h1>", unsafe_allow_html=True)
        st.markdown("<p style='text-align: center; color: #666; font-size: 1.15rem; font-weight: 300;'>The Premier Resource for Nucleic Acid Ligand Structures</p>", unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)

    # About Section
    st.markdown("""
    <div style='background-color: #f8f9fa; padding: 30px; border-radius: 12px; border-left: 5px solid #4CAF50; box-shadow: 0 4px 15px rgba(0,0,0,0.05); margin-bottom: 40px;'>
        <h3 style='color: #2c3e50; margin-top: 0;'>About the Database</h3>
        <p style='color: #444; font-size: 1.05rem; line-height: 1.6;'>
            The <b>NucLigs Database</b> is a specialized repository designed to facilitate research in the field of nucleic acid targeting. 
            It aggregates structural data of small molecule ligands bound to DNA and RNA targets, providing a unified platform for detailed analysis and visualization.
        </p>
        <p style='color: #444; font-size: 1.05rem; line-height: 1.6; margin-bottom: 0;'>
            Whether you are involved in <b>rational drug design</b>, <b>structural biology</b>, or <b>cheminformatics</b>, NucLigs offers robust tools to explore the physico-chemical landscape of nucleic acid interactions, supporting the discovery of next-generation therapeutics.
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # Feature Cards
    f1, f2, f3 = st.columns(3)
    
    with f1:
        st.markdown("""
        <div class="home-card">
            <h3>🧊 3D Visualization</h3>
            <p>Interactive, high-fidelity rendering of ligand-target complexes using Py3Dmol. Inspect binding modes, molecular surfaces, and structural conformations in real-time directly within your browser.</p>
        </div>
        """, unsafe_allow_html=True)
        
    with f2:
        st.markdown("""
        <div class="home-card">
            <h3>⚗️ Chemical Profiling</h3>
            <p>Automated calculation of critical molecular descriptors. Access data on Molecular Weight, LogP, TPSA, and Lipinski's Rule of 5 compliance powered by the RDKit cheminformatics engine.</p>
        </div>
        """, unsafe_allow_html=True)
        
    with f3:
        st.markdown("""
        <div class="home-card">
            <h3>📂 Data Accessibility</h3>
            <p>Seamlessly retrieve standardized structural data. Export ligands and complexes in industry-standard formats (PDB, SDF, MOL) to integrate directly with your local modeling workflows.</p>
        </div>
        """, unsafe_allow_html=True)
        
    st.markdown("<br><br>", unsafe_allow_html=True)

    # Primary Action
    _, btn_col, _ = st.columns([1.5, 1, 1.5])
    with btn_col:
        if st.button("🚀 Explore the Collection", type="primary", use_container_width=True):
            st.session_state['page'] = 'database'
            st.rerun()
    
    st.markdown("<div style='text-align: center; margin-top: 50px; color: #aaa; font-size: 0.85rem;'>© 2024 NucLigs Database Project • Version 2.1</div>", unsafe_allow_html=True)

def render_database():
    metadata_df = load_metadata()
    all_nuc_ids = get_ids_from_metadata()

    # Sidebar
    with st.sidebar:
        if st.button("⬅️ Back to Home"):
            st.session_state['page'] = 'home'
            st.rerun()
            
        st.markdown("---")
        st.markdown("### 🔍 Finder")

        search_query = st.text_input("Filter database:", placeholder="Search NucL ID or Name...", label_visibility="collapsed")
        
        # Filter Logic
        if search_query:
            # Filter matches from both 'nl' and 'names' columns
            q = search_query.lower()
            mask = (
                metadata_df['nl'].astype(str).str.lower().str.contains(q, na=False) | 
                metadata_df['names'].astype(str).str.lower().str.contains(q, na=False)
            )
            filtered_matches = metadata_df[mask]['nl'].dropna().unique().tolist()
            nuc_ids = sorted([str(x) for x in filtered_matches])
            
            if not nuc_ids:
                st.warning("No matches found.")
                nuc_ids = []
            else:
                st.success(f"Found {len(nuc_ids)} structures")
        else:
            nuc_ids = all_nuc_ids

        # Selection
        if nuc_ids:
            selected_nuc_id = st.selectbox("Select Structure Result:", nuc_ids, index=0)
        else:
            selected_nuc_id = None
        
        # --- BULK ACTIONS (New Feature) ---
        if nuc_ids:
            with st.expander("📦 Bulk Actions", expanded=False):
                st.caption("Download multiple structures & data based on your current search.")
                
                # Selection for bulk
                default_sel = nuc_ids[:10] if len(nuc_ids) < 20 else []
                bulk_selected = st.multiselect("Select structures to download:", nuc_ids, default=default_sel)
                
                fmt = st.selectbox("Format", [".pdb", ".sdf", ".mol"], index=0)
                include_preds = st.checkbox("Compute Features (Slower)", value=False, help="Runs RDKit calculation for every selected structure.")
                
                if st.button("Generate Download Package", use_container_width=True, disabled=len(bulk_selected)==0):
                    status = st.status("Processing bulk download...", expanded=True)
                    
                    try:
                        # 1. Prepare ZIP and CSV
                        zip_buffer = io.BytesIO()
                        csv_data_list = []
                        
                        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
                            progress_bar = status.progress(0)
                            total = len(bulk_selected)
                            
                            for idx, nid in enumerate(bulk_selected):
                                # 1. Fetch Data
                                row = metadata_df[metadata_df['nl'].astype(str) == nid]
                                if row.empty: continue
                                
                                meta_dict = row.iloc[0].to_dict()
                                pdb_fname = str(meta_dict.get('pdbs', ''))
                                content = fetch_pdb_from_supabase(pdb_fname)
                                
                                if not content: continue
                                
                                # 2. Process Structure for Format
                                final_content = content
                                final_ext = fmt
                                mol_obj = None
                                
                                # Convert if needed (SDF/MOL) or if Features requested
                                if fmt != ".pdb" or include_preds:
                                    smiles = str(meta_dict.get('smiles', ''))
                                    try:
                                        if smiles and smiles.lower() != 'nan':
                                            mol_obj = Chem.MolFromSmiles(smiles)
                                        else:
                                            mol_obj = Chem.MolFromPDBBlock(content, sanitize=True, removeHs=False)
                                    except:
                                        mol_obj = None

                                if fmt in [".sdf", ".mol"] and mol_obj:
                                    final_content = Chem.MolToMolBlock(mol_obj)
                                elif fmt in [".sdf", ".mol"] and not mol_obj:
                                    # Fallback if conversion fails
                                    final_content = content # Keep PDB
                                    final_ext = ".pdb"

                                # 3. Add to Zip
                                zip_file.writestr(f"{nid}{final_ext}", final_content)
                                
                                # 4. Compute Features if requested
                                if include_preds and mol_obj:
                                    preds = compute_physchem(mol_obj)
                                    # Remove non-serializable RDKit obj
                                    preds.pop("_RDKitMol", None)
                                    meta_dict.update(preds)
                                
                                csv_data_list.append(meta_dict)
                                progress_bar.progress((idx + 1) / total)
                        
                        status.write("Compressing files...")
                        zip_buffer.seek(0)
                        
                        # Generate CSV buffer
                        csv_buffer = None
                        if csv_data_list:
                            df_bulk = pd.DataFrame(csv_data_list)
                            csv_buffer = df_bulk.to_csv(index=False).encode('utf-8')
                        
                        status.update(label="Ready for Download!", state="complete", expanded=True)
                        
                        # Show Download Buttons
                        st.download_button(
                            label=f"📥 Download {len(bulk_selected)} Structures (.zip)",
                            data=zip_buffer,
                            file_name="nucligs_structures.zip",
                            mime="application/zip",
                            use_container_width=True
                        )
                        
                        if csv_buffer:
                            st.download_button(
                                label="📥 Download Data Table (.csv)",
                                data=csv_buffer,
                                file_name="nucligs_data.csv",
                                mime="text/csv",
                                use_container_width=True
                            )
                            
                    except Exception as e:
                        status.update(label="Error processing download", state="error")
                        st.error(f"An error occurred: {str(e)}")

        st.divider()
        st.markdown("**Viewer Settings**")
        bg_mode = st.radio("Background", ["Light", "Dark"], horizontal=True, label_visibility="collapsed")
        bg_color = "white" if bg_mode == "Light" else "#1e1e1e"
        st.divider()
        st.caption(f"**Total Entries:** {len(all_nuc_ids)}")

    # Main Content
    if not selected_nuc_id:
        st.info("👈 Please search for or select a structure from the sidebar.")
        return

    # Resolve Data
    row = metadata_df[metadata_df['nl'].astype(str) == selected_nuc_id]
    
    pdb_text = None
    data = {}
    smiles_str = None
    
    if not row.empty:
        pdb_filename = str(row.iloc[0]['pdbs'])
        smiles_str = str(row.iloc[0].get('smiles', ''))
        data = row.iloc[0].to_dict()
        pdb_text = fetch_pdb_from_supabase(pdb_filename)
    else:
        st.error(f"Metadata not found for ID: {selected_nuc_id}")
        return

    if pdb_text:
        # Create Molecule
        mol = None
        if smiles_str and smiles_str.lower() != 'nan' and smiles_str.strip():
            try:
                mol = Chem.MolFromSmiles(smiles_str)
            except:
                mol = None
        
        if mol is None:
            try:
                mol = Chem.MolFromPDBBlock(pdb_text, sanitize=True, removeHs=False)
            except:
                mol = None

        physchem = compute_physchem(mol)
        
        # UI Columns
        col_viewer, col_preds, col_meta = st.columns([2.2, 0.9, 0.9])

        # 1. Viewer
        with col_viewer:
            st.subheader(f"3D Visualization: {selected_nuc_id}")
            show_3d_pdb(pdb_text, bg_color)
            
            st.markdown("##### 📥 Export Data")
            d1, d2, d3, d4 = st.columns([1, 1, 1, 1.2]) # Adjusted columns
            with d1:
                st.download_button("Download .PDB", pdb_text, f"{selected_nuc_id}.pdb", "chemical/x-pdb", use_container_width=True)
            
            mol_obj = physchem.get("_RDKitMol")
            if mol_obj:
                sdf_data = Chem.MolToMolBlock(mol_obj)
                with d2: st.download_button("Download .SDF", sdf_data, f"{selected_nuc_id}.sdf", "chemical/x-mdl-sdfile", use_container_width=True)
                with d3: st.download_button("Download .MOL", sdf_data, f"{selected_nuc_id}.mol", "chemical/x-mdl-molfile", use_container_width=True)
            else:
                with d2: st.button("SDF Unavail.", disabled=True, use_container_width=True)
                with d3: st.button("MOL Unavail.", disabled=True, use_container_width=True)
            
            # Single CSV Download
            with d4:
                # Merge data for download
                full_data = {**data, **physchem}
                full_data.pop("_RDKitMol", None)
                csv_single = pd.DataFrame([full_data]).to_csv(index=False).encode('utf-8')
                st.download_button("📄 Data Profile (.csv)", csv_single, f"{selected_nuc_id}_data.csv", "text/csv", use_container_width=True)

        # 2. Analysis
        with col_preds:
            st.subheader("Chemical Analysis")
            st.markdown('<div class="meta-scroll">', unsafe_allow_html=True)
            if physchem:
                st.markdown('<div class="feature-card"><h5>🧪 Identity</h5>', unsafe_allow_html=True)
                render_row("Formula", physchem.get("Formula", "-"))
                render_row("Mol Weight", f"{physchem.get('Mol Wt', '-')} da")
                render_row("Formal Charge", physchem.get("Charge", "0"))
                render_row("Stereocenters", physchem.get("Chiral Centers", "0"))
                st.markdown('</div>', unsafe_allow_html=True)

                violations = physchem.get("Lipinski Violations", 0)
                badge_class = "badge-pass" if violations == 0 else "badge-fail"
                badge_text = "PASS (0 Violations)" if violations == 0 else f"FAIL ({violations} Violations)"
                
                st.markdown(f'<div class="feature-card"><div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:10px; border-bottom: 2px solid #f0f2f6; padding-bottom:8px;"><h5 style="margin:0; border:none; padding:0;">⚖️ Rule of 5</h5><span class="{badge_class}">{badge_text}</span></div>', unsafe_allow_html=True)
                render_row("LogP", physchem.get("LogP", "-"), "≤ 5")
                render_row("H-Donors", physchem.get("H-Don", "-"), "≤ 5")
                render_row("H-Acceptors", physchem.get("H-Acc", "-"), "≤ 10")
                render_row("Rot. Bonds", physchem.get("Rot. Bonds", "-"))
                render_row("Mol Weight", f"{physchem.get('Mol Wt', '-')} da", "≤ 500")
                st.markdown('</div>', unsafe_allow_html=True)

                st.markdown('<div class="feature-card"><h5>💊 Druglikeness</h5>', unsafe_allow_html=True)
                render_row("QED Score", physchem.get("QED", "-"))
                render_row("Est. Solubility", physchem.get("ESOL (LogS)", "-"))
                render_row("TPSA", f"{physchem.get('TPSA', '-')} Å²")
                render_row("Fraction Csp3", physchem.get("F-Csp3", "-"))
                st.markdown('</div>', unsafe_allow_html=True)
                
                st.markdown('<div class="feature-card"><h5>💍 Ring Systems</h5>', unsafe_allow_html=True)
                render_row("Aromatic Rings", physchem.get("Arom. Rings", "-"))
                render_row("Saturated Rings", physchem.get("Sat. Rings", "-"))
                st.markdown('</div>', unsafe_allow_html=True)
            else:
                st.warning("Unable to compute chemical properties.")
            st.markdown("</div>", unsafe_allow_html=True)

        # 3. Metadata
        with col_meta:
            st.subheader("Metadata Record")
            st.markdown('<div class="meta-scroll">', unsafe_allow_html=True)
            if data:
                nl_id = data.get('nl', 'Unknown')
                chem_name = data.get('names', data.get('name', '')) 
                
                st.markdown(f'<div class="id-card"><div class="id-label">NucLigs Identifier</div><div class="id-value">{nl_id}</div><div class="id-sub">{chem_name}</div></div>', unsafe_allow_html=True)

                st.markdown('<div class="feature-card"><h5>📋 General Info</h5>', unsafe_allow_html=True)
                exclude_fields = ['nl', 'names', 'name', 'pdbs', 'match', 'smiles', 'inchi', 'description', 'sequence']
                for key, value in data.items():
                    if key not in exclude_fields and str(value).lower() != 'nan':
                        render_row(key.replace('_', ' ').title(), str(value))
                st.markdown('</div>', unsafe_allow_html=True)

                long_fields = ["smiles", "inchi", "description", "sequence"]
                for key in long_fields:
                    if key in data and str(data[key]).lower() != 'nan':
                        with st.expander(key.upper(), expanded=False):
                            st.code(str(data[key]), language="text")
            else:
                st.info("No metadata found.")
            st.markdown("</div>", unsafe_allow_html=True)

# ----------------------------------------------------
# MAIN ROUTER
# ----------------------------------------------------

# Initialize Session State
if 'page' not in st.session_state:
    st.session_state['page'] = 'home'

# Route
if st.session_state['page'] == 'home':
    render_homepage()
else:
    render_database()
