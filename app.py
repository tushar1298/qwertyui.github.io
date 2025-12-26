import streamlit as st
from supabase import create_client

# --------------------------------------------------
# CONFIG (use Streamlit Secrets)
# --------------------------------------------------
SUPABASE_URL = "https://heuzgnhlrumyfcfigoon.supabase.co"
SUPABASE_KEY = "sb_secret_UuFsAopmAmHrdvHf6-mGBg_X0QNgMF5"

BUCKET_NAME = "codes"
REMOTE_APP_FILENAME = "app.py"   # this is the REAL app in Supabase

# --------------------------------------------------
# LOAD REMOTE APP CODE
# --------------------------------------------------
@st.cache_resource
def load_remote_source():
    supabase = create_client(SUPABASE_URL, SUPABASE_KEY)
    data = supabase.storage.from_(BUCKET_NAME).download(REMOTE_APP_FILENAME)
    return data.decode("utf-8")

# --------------------------------------------------
# MAIN LOADER
# --------------------------------------------------
def main():
    try:
        source_code = load_remote_source()

        # 🔐 Execute in SAFE, ISOLATED namespace
        exec(
            source_code,
            {
                "__name__": "__main__",  # prevents recursion
                "st": st,                # Streamlit context
            },
        )

    except Exception as e:
        st.error("❌ Failed to load application from Supabase")
        st.exception(e)

# --------------------------------------------------
main()
