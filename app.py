import streamlit as st
from supabase import create_client

# --------------------------------------------------------------------------------
# GITHUB LOADER SCRIPT
# --------------------------------------------------------------------------------
# This script serves as a shell to fetch and run the actual application 
# stored securely in your Supabase Storage. Upload this file as 'app.py' 
# to your GitHub repository.
# --------------------------------------------------------------------------------

# Configuration
# Tip: For extra security in production, store these in Streamlit Secrets (st.secrets)

SUPABASE_URL = "https://heuzgnhlrumyfcfigoon.supabase.co"
SUPABASE_KEY = "sb_secret_UuFsAopmAmHrdvHf6-mGBg_X0QNgMF5"



BUCKET_NAME = "codes"
SOURCE_FILENAME = "app.py"

def main():
    try:
        # 1. Initialize Supabase Client
        supabase = create_client(SUPABASE_URL, SUPABASE_KEY)

        # 2. Fetch the source code from the 'codes' bucket
        print(f"Fetching {SOURCE_FILENAME} from Supabase...")
        response = supabase.storage.from_(BUCKET_NAME).download(SOURCE_FILENAME)
        
        # 3. Decode the byte response to a UTF-8 string
        source_code = response.decode('utf-8')

        # 4. Execute the fetched code
        #    We pass 'globals()' to ensure the executed code has access to 
        #    the necessary imports and Streamlit context.
        exec(source_code, globals())

    except Exception as e:
        # If fetching fails, show a generic error (or specific one for debugging)
        st.error("Failed to load the application from the remote server.")
        st.expander("Error Details").code(str(e))

if __name__ == "__main__":
    main()
