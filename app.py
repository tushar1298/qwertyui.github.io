import io
import requests
import pandas as pd
import streamlit as st
import py3Dmol

# ------------- PAGE CONFIG -----------------
st.set_page_config(
    page_title="NucLigs PDB Extractor",
    page_icon="🧬",
    layout="wide",
)

st.title("🧬 NucLigs PDB Extractor (GitHub + Streamlit)")

st.markdown(
    """
Upload your NucLigs Excel/CSV file (with a **`pdbs`** header) and this app will:

1. Read the table and let you search entries.
2. For a selected row, build the PDB file URL from your **GitHub raw** base path.
3. Download the PDB, let you **view it in 3D** and **download** it.
"""
)

# ------------- CONFIG: GITHUB RAW BASE -----------------
st.sidebar.header("GitHub PDB Source")

default_base = "https://raw.githubusercontent.com/<USER>/<REPO>/<BRANCH>/pdbs"
github_base = st.sidebar.text_input(
    "GitHub raw base URL (directory that contains your .pdb files)",
    value=default_base,
    help=(
        "Example:\n"
        "https://raw.githubusercontent.com/username/NucLigs-PDBs/main/pdbs\n\n"
        "The app will append '<value_in_pdbs_column>.pdb' to this URL."
    ),
)

if github_base.endswith("/"):
    github_base = github_base[:-1]

# ------------- FILE UPLOAD -----------------
st.sidebar.header("Upload database file")
uploaded = st.sidebar.file_uploader(
    "Excel/CSV with 'pdbs' column", type=["xlsx", "xls", "csv"]
)

@st.cache_data
def load_table(file) -> pd.DataFrame:
    """Load Excel/CSV into DataFrame."""
    if file.name.lower().endswith((".xlsx", ".xls")):
        return pd.read_excel(file)
    else:
        return pd.read_csv(file)

def find_column(columns, targets):
    """Return first existing column whose lower-stripped name is in targets."""
    lower_map = {c.lower().strip(): c for c in columns}
    for t in targets:
        if t in lower_map:
            return lower_map[t]
    return None

df = None
if uploaded:
    try:
        df = load_table(uploaded)
    except Exception as e:
        st.error(f"Error reading file: {e}")

if df is None:
    st.info("⬅ Upload your Excel/CSV in the sidebar to begin.")
    st.stop()

# ------------- COLUMN DETECTION (like your HTML) -------------
cols = df.columns

pdb_col = find_column(cols, {"pdbs"})
if pdb_col is None:
    st.error("Could not find a 'pdbs' column in the uploaded file. "
             "Please check the header name.")
    st.stop()

nl_col = find_column(cols, {"nl"})
name_col = find_column(cols, {"names", "name"})

st.success(
    f"Detected columns → pdbs: **{pdb_col}** "
    f"{'· NL: **' + nl_col + '**' if nl_col else ''} "
    f"{'· names: **' + name_col + '**' if name_col else ''}"
)

# ------------- SEARCH + SELECTION UI -----------------
left, right = st.columns([1, 2], gap="large")

with left:
    st.subheader("Search & select")

    search_text = st.text_input("Search (NL / name / pdbs)", "")
    mask = pd.Series(True, index=df.index)

    if search_text:
        s = search_text.lower()
        parts = []
        if nl_col:
            parts.append(df[nl_col].astype(str).str.lower().str.contains(s, na=False))
        if name_col:
            parts.append(df[name_col].astype(str).str.lower().str.contains(s, na=False))
        parts.append(df[pdb_col].astype(str).str.lower().str.contains(s, na=False))
        mask = parts[0]
        for m in parts[1:]:
            mask |= m

    filtered = df[mask].copy()
    st.caption(f"{len(filtered)} entries matched.")

    # Build labels for selectbox
    def make_label(idx):
        row = filtered.loc[idx]
        nl = str(row[nl_col]) if nl_col else ""
        nm = str(row[name_col]) if name_col else ""
        pdb_id = str(row[pdb_col])
        label_bits = [b for b in [nl, nm, pdb_id] if b]
        return " | ".join(label_bits)

    selected_idx = st.selectbox(
        "Choose a row",
        options=filtered.index.tolist(),
        format_func=make_label,
    )

    row = filtered.loc[selected_idx]

with right:
    st.subheader("Selected entry")

    # Show metadata
    st.markdown("**Row details**")
    st.dataframe(row.to_frame("value"))

    # Determine PDB filename
    raw_pdb_id = str(row[pdb_col]).strip()
    if raw_pdb_id.lower().endswith(".pdb"):
        pdb_filename = raw_pdb_id
    else:
        pdb_filename = raw_pdb_id + ".pdb"

    pdb_url = f"{github_base}/{pdb_filename}"

    st.markdown(f"**PDB ID:** `{raw_pdb_id}`  →  **URL:** `{pdb_url}`")

    # ------------- FETCH PDB FROM GITHUB -------------
    fetch = st.button("Fetch PDB from GitHub")

    pdb_text = None
    if fetch:
        if "<USER>" in github_base:
            st.error("Please edit the GitHub base URL to your real repository first.")
        else:
            with st.spinner("Downloading PDB from GitHub..."):
                try:
                    resp = requests.get(pdb_url, timeout=15)
                    if resp.status_code == 200:
                        pdb_text = resp.text
                        st.success("PDB fetched successfully ✅")
                    else:
                        st.error(f"GitHub returned status {resp.status_code}. "
                                 "Check filename or base URL.")
                except Exception as e:
                    st.error(f"Error fetching PDB: {e}")

    # Keep PDB string between reruns using session_state
    if pdb_text is not None:
        st.session_state["last_pdb"] = pdb_text
        st.session_state["last_pdb_name"] = pdb_filename

    if "last_pdb" in st.session_state:
        pdb_text = st.session_state["last_pdb"]
        pdb_filename = st.session_state.get("last_pdb_name", pdb_filename)

        st.markdown("### Download PDB")
        st.download_button(
            "Download PDB file",
            data=pdb_text,
            file_name=pdb_filename,
            mime="chemical/x-pdb",
        )

        st.markdown("### 3D View")
        try:
            view = py3Dmol.view(width=600, height=450)
            view.addModel(pdb_text, "pdb")
            view.setStyle({"stick": {}})
            view.zoomTo()
            html = view._make_html()
            st.components.v1.html(html, height=480)
        except Exception as e:
            st.warning(f"Could not render 3D view: {e}")
