import streamlit as st
from supabase import create_client

SUPABASE_URL = "https://heuzgnhlrumyfcfigoon.supabase.co"
SUPABASE_KEY = st.secrets["SUPABASE_KEY"]

BUCKET_NAME = "codes"
SOURCE_FILENAME = "app1.py"

@st.cache_resource
def load_remote_app():
    supabase = create_client(SUPABASE_URL, SUPABASE_KEY)
    data = supabase.storage.from_(BUCKET_NAME).download(SOURCE_FILENAME)
    return data.decode("utf-8")

def main():
    try:
        source_code = load_remote_app()

        # 🔴 VERY IMPORTANT: execute in isolated namespace
        exec(source_code, {"st": st})

    except Exception as e:
        st.error("Failed to load remote application")
        st.exception(e)

main()
