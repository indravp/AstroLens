import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014, NFWPotential, LogarithmicHaloPotential, PlummerPotential
import astropy.units as u
import numpy as np

st.set_page_config(page_title="AstroLens", layout="wide")

st.title("🔭 AstroLens")
st.subheader("Visual Toolkit for Stellar Cluster Analysis")

st.sidebar.header("Upload Gaia CSV File")
file = st.sidebar.file_uploader("Choose a Gaia DR3 CSV", type=["csv"])

# --- Orbit Simulation ---
st.markdown("### 🌀 Orbit Simulation")

with st.expander("Input Initial Cluster Parameters"):
    ra = st.number_input("RA (deg)", value=250.423)
    dec = st.number_input("Dec (deg)", value=-26.537)
    dist = st.number_input("Distance (kpc)", value=35.0)
    pmra_val = st.number_input("Proper Motion RA (mas/yr)", value=-2.5)
    pmdec_val = st.number_input("Proper Motion DEC (mas/yr)", value=-1.8)
    vrad = st.number_input("Radial Velocity (km/s)", value=30.0)
    direction = st.radio("Orbit Direction", options=["Forward", "Backward"])
    integration_time_myr = st.slider("Integration time |t| (Million Years)", 10, 10000, 1000, step=10)
    projection = st.selectbox("Choose Orbit Projection:", 
                              options=["x-y", "R-z", "x-z"])
    # Galactic Potential Selector
    potential_choice = st.selectbox(
        "Choose Galactic Potential",
        ("MWPotential2014", "NFWPotential", "LogarithmicHaloPotential", "PlummerPotential"))

    # Initialize the selected potential
    if potential_choice == "MWPotential2014":
        selected_potential = MWPotential2014
    elif potential_choice == "NFWPotential":
        selected_potential = NFWPotential(a=16, amp=0.35)
    elif potential_choice == "LogarithmicHaloPotential":
        selected_potential = LogarithmicHaloPotential(q=1.0, core=1.0)
    elif potential_choice == "PlummerPotential":
        selected_potential = PlummerPotential(amp=1.0, b=0.5)

    # Dictionary of descriptions
    potential_descriptions = {
        "MWPotential2014": "A composite potential representing the Milky Way, including bulge, disk, and dark matter halo components. Used widely in galpy-based orbit studies.",
        "NFWPotential": "Navarro-Frenk-White profile — a common dark matter halo model derived from cosmological simulations. Cuspy inner density profile.",
        "LogarithmicHaloPotential": "Simple model of a flattened dark matter halo with a logarithmic potential. Often used for exploring basic halo dynamics.",
        "PlummerPotential": "Spherically symmetric, used to model globular clusters or soft-core galaxies. Smooth, finite central density."
    }

    # Show description
    st.markdown(f"**Description:** {potential_descriptions[potential_choice]}")

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
    
    # User input for integration time
    integration_time_gyr = integration_time_myr / 1000

    if direction == "Forward":
        ts = np.linspace(0, integration_time_gyr, 1000) * u.Gyr
    else:
        ts = np.linspace(0, -integration_time_gyr, 1000) * u.Gyr
    
    o.integrate(ts, selected_potential)
    
    # Plot
    if projection == "x-y":
        coordinates = ['x', 'y']
        x_vals = o.x(ts)
        y_vals = o.y(ts)
        xlabel, ylabel = "X (kpc)", "Y (kpc)"
    elif projection == "R-z":
        coordinates = ['R', 'z']
        x_vals = o.R(ts)
        y_vals = o.z(ts)
        xlabel, ylabel = "R (kpc)", "Z (kpc)"
    else:  # x-z
        coordinates = ['x', 'z']
        x_vals = o.x(ts)
        y_vals = o.z(ts)
        xlabel, ylabel = "X (kpc)", "Z (kpc)"
    
    t = ts.to(u.Gyr).value

    orbit_data = pd.DataFrame({
        'Time (Gyr)': t,
        'X (kpc)': o.x(ts),
        'Y (kpc)': o.y(ts),
        'Z (kpc)': o.z(ts),
        'R (kpc)': o.R(ts),
        'VX (km/s)': o.vx(ts),
        'VY (km/s)': o.vy(ts),
        'VZ (km/s)': o.vz(ts)
    })

    orbit_csv = orbit_data.to_csv(index=False)


    fig, ax = plt.subplots(figsize=(6, 5))
    sc = ax.scatter(x_vals, y_vals, c=t, cmap='viridis', s=2)
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('Time (Gyr)', fontsize=12)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"Orbit Projection: {projection} over {integration_time_myr} Myr")
    st.pyplot(fig)

    st.download_button("Download Orbit Data", orbit_csv, "orbit_data.csv", "text/csv")

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

    # --- Download Filtered Data ---
    csv = filtered.to_csv(index=False).encode('utf-8')
    st.download_button("Download Filtered Data", csv, "filtered_gaia.csv", "text/csv")

else:
    st.warning("👈 Please open the sidebar(if closed) to upload a Gaia CSV file to begin.")
