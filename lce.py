import requests
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from io import StringIO

def get_ztf_lightcurve(ra_deg, dec_deg, radius_arcsec=2, band="r", fmt="csv"):
    """
    Query ZTF IRSA API for lightcurves around given RA, Dec.
    """
    base_url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"
    radius_deg = radius_arcsec / 3600.0  # arcsec → deg
    params = {
        "POS": f"CIRCLE {ra_deg} {dec_deg} {radius_deg}",
        "BANDNAME": band,
        "FORMAT": fmt.lower()
    }
    response = requests.get(base_url, params=params)
    response.raise_for_status()
    return pd.read_csv(StringIO(response.text), comment="#")

# -------------------------------
# Load your catalogue
# -------------------------------
catalogue = pd.read_csv("magnetic_cataclysmic_variables.txt")

all_lightcurves = []

# ✅ Loop over ALL objects (remove .head(20))
for _, row in catalogue.iterrows():
    name = row["Name"]
    ra_str, dec_str = row["RA(J2000)"], row["DEC(J2000)"]

    try:
        # Convert RA/Dec to decimal degrees
        coord = SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg))
        ra_deg, dec_deg = coord.ra.deg, coord.dec.deg

        print(f"Querying ZTF for {name} at RA={ra_deg:.4f}, Dec={dec_deg:.4f}")
        df = get_ztf_lightcurve(ra_deg, dec_deg, radius_arcsec=2, band="r")

        if not df.empty:
            df["object"] = name
            all_lightcurves.append(df)
            print(f"  ✅ Got {len(df)} points for {name}")
        else:
            print(f"  ⚠ No data returned for {name}")

    except Exception as e:
        print(f"  ❌ Failed for {name}: {e}")

# -------------------------------
# Save everything into one CSV
# -------------------------------
if all_lightcurves:
    merged = pd.concat(all_lightcurves, ignore_index=True)
    merged.to_csv("lightcurves_all_r.csv", index=False)
    print(f"\n✅ Saved r-band lightcurves for {len(catalogue)} stars ({len(merged)} total points)")
else:
    print("\n❌ No lightcurves fetched.")
