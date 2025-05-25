import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014
import astropy.units as u
import numpy as np

st.set_page_config(page_title="AstroLens", layout="wide")

st.title("🔭 AstroLens")
st.subheader("Visual Toolkit for Stellar Cluster Analysis")

# --- File Upload ---
st.sidebar.header("Upload Gaia CSV File")
file = st.sidebar.file_uploader("Choose a Gaia DR3 CSV", type=["csv"])

if file:
    df = pd.read_csv(file)

    # --- Column Selection (Assumed Gaia Format) ---
    col1, col2 = st.columns(2)
    with col1:
        x_col = st.selectbox("Color (x-axis)", df.columns, index=df.columns.get_loc('bp_rp') if 'bp_rp' in df.columns else 0)
        y_col = st.selectbox("Magnitude (y-axis)", df.columns, index=df.columns.get_loc('phot_g_mean_mag') if 'phot_g_mean_mag' in df.columns else 1)
    with col2:
        pmra_col = st.selectbox("μ_α*", df.columns, index=df.columns.get_loc('pmra') if 'pmra' in df.columns else 0)
        pmdec_col = st.selectbox("μ_δ", df.columns, index=df.columns.get_loc('pmdec') if 'pmdec' in df.columns else 1)
        parallax_col = st.selectbox("Parallax", df.columns, index=df.columns.get_loc('parallax') if 'parallax' in df.columns else 2)

    # --- Filtering Options ---
    st.sidebar.subheader("Filters")
    pmra_range = st.sidebar.slider("Proper Motion RA (μ_α*)", float(df[pmra_col].min()), float(df[pmra_col].max()), (float(df[pmra_col].min()), float(df[pmra_col].max())))
    pmdec_range = st.sidebar.slider("Proper Motion DEC (μ_δ)", float(df[pmdec_col].min()), float(df[pmdec_col].max()), (float(df[pmdec_col].min()), float(df[pmdec_col].max())))
    parallax_range = st.sidebar.slider("Parallax", float(df[parallax_col].min()), float(df[parallax_col].max()), (float(df[parallax_col].min()), float(df[parallax_col].max())))

    filtered = df[
        (df[pmra_col].between(*pmra_range)) &
        (df[pmdec_col].between(*pmdec_range)) &
        (df[parallax_col].between(*parallax_range))
    ]

    # --- CMD Plot ---
    st.markdown("### 📈 Color-Magnitude Diagram (CMD)")
    fig1, ax1 = plt.subplots()
    ax1.scatter(filtered[x_col], filtered[y_col], s=5, alpha=0.6, color='navy')
    ax1.set_xlabel(f"{x_col}")
    ax1.set_ylabel(f"{y_col}")
    ax1.invert_yaxis()
    st.pyplot(fig1)

    # --- Proper Motion Plot ---
    st.markdown("### 💫 Proper Motion Plot")
    fig2, ax2 = plt.subplots()
    ax2.scatter(filtered[pmra_col], filtered[pmdec_col], s=5, alpha=0.6, color='darkgreen')
    ax2.set_xlabel(f"{pmra_col}")
    ax2.set_ylabel(f"{pmdec_col}")
    st.pyplot(fig2)

    # --- Orbit Simulation ---
    st.markdown("### 🌀 Orbit Simulation")

    with st.expander("Input Initial Cluster Parameters"):
        ra = st.number_input("RA (deg)", value=250.423)
        dec = st.number_input("Dec (deg)", value=-26.537)
        dist = st.number_input("Distance (kpc)", value=35.0)
        pmra_val = st.number_input("Proper Motion RA (mas/yr)", value=-2.5)
        pmdec_val = st.number_input("Proper Motion DEC (mas/yr)", value=-1.8)
        vrad = st.number_input("Radial Velocity (km/s)", value=30.0)
        run_sim = st.button("Simulate Orbit")

    if run_sim:
        # Convert to galpy-compatible units
        ra_u = ra * u.deg
        dec_u = dec * u.deg
        dist_u = dist * u.kpc
        pmra_u = pmra_val * u.mas/u.yr
        pmdec_u = pmdec_val * u.mas/u.yr
        vrad_u = vrad * u.km/u.s

        # Create orbit object
        o = Orbit(vxvv=[ra_u, dec_u, dist_u, pmra_u, pmdec_u, vrad_u], radec=True)

        # Integrate orbit forward in time
        ts = np.linspace(0, 2, 1000) * u.Gyr
        o.integrate(ts, MWPotential2014)

        st.write("Final position:", o.x(ts[-1]), o.y(ts[-1]))
        # Plot
        fig3 = plt.figure(figsize=(6, 6))
        o.plot(d1='x', d2='y', overplot=True)
        # Set axis limits based on your orbit's extent
        x_vals = o.x(ts)
        y_vals = o.y(ts)

        x_min, x_max = np.min(x_vals), np.max(x_vals)
        y_min, y_max = np.min(y_vals), np.max(y_vals)

        # Add some padding
        padding = 0.1 * max(x_max - x_min, y_max - y_min)
        plt.xlim(x_min - padding, x_max + padding)
        plt.ylim(y_min - padding, y_max + padding)

        plt.xlabel("X (kpc)")
        plt.ylabel("Y (kpc)")
        plt.title("Orbit in XY Plane")
        plt.gca().set_aspect('equal', adjustable='box')

        plt.locator_params(axis='x', nbins=10)  # or fewer
        plt.locator_params(axis='y', nbins=10)


        st.pyplot(fig3)

    # --- Download Filtered Data ---
    csv = filtered.to_csv(index=False).encode('utf-8')
    st.download_button("Download Filtered Data", csv, "filtered_gaia.csv", "text/csv")

else:
    st.warning("👈 Please upload a Gaia CSV file to begin.")