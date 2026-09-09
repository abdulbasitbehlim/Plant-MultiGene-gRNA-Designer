import streamlit as st

st.set_page_config(page_title="Live CRISPR Apps", page_icon="🔗", layout="centered")

st.title("🔗 Live CRISPR Apps")
st.write("Use the buttons below to open the deployed tools or inspect their source code on GitHub.")

st.subheader("🌿 Plant MultiGene gRNA Designer")
st.write("Design and validate shared SpCas9 guides for multiple homologous plant genes.")
st.link_button("Open Plant MultiGene App", "https://plant-multigene-grna-designer.onrender.com/", use_container_width=True)
st.link_button("Plant MultiGene GitHub", "https://github.com/abdulbasitbehlim/Plant-MultiGene-gRNA-Designer", use_container_width=True)

st.divider()

st.subheader("🧬 OpenCRISPR-1 gRNA Designer")
st.write("Design and validate conservative NGG-compatible guide RNAs for OpenCRISPR-1.")
st.link_button("Open OpenCRISPR-1 App", "https://opencrispr1-grna-designer.onrender.com/", use_container_width=True)
st.link_button("OpenCRISPR-1 GitHub", "https://github.com/abdulbasitbehlim/OpenCRISPR1-gRNA-Designer", use_container_width=True)

st.info("Render free services may sleep after inactivity, so the first load can take longer while the service wakes up.")
