#  ==========================================================
# Import required Python libraries
# ==========================================================
import requests
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from io import StringIO

# ==========================================================
# Function to query ZTF (Zwicky Transient Facility) lightcurves
# ==========================================================
def get_ztf_lightcurve(ra_deg, dec_deg, radius_arcsec=1.2, band="r", fmt="csv"):
    """
    Query the ZTF IRSA API for lightcurves around a given RA, Dec position.

    Parameters:
        ra_deg (float): Right Ascension in degrees
        dec_deg (float): Declination in degrees
        radius_arcsec (float): Search radius around target (arcseconds)
        band (str): Photometric band to fetch ('g', 'r', 'i')
        fmt (str): Response format (default: CSV)
    """
    base_url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"
    radius_deg = radius_arcsec / 3600.0  # convert arcsec → degrees

    params = {
        "POS": f"CIRCLE {ra_deg} {dec_deg} {radius_deg}",
        "BANDNAME": band,
        "FORMAT": fmt.lower(),
    }

    # Fetch data
    response = requests.get(base_url, params=params, timeout=60) # Send request to IRSA ZTF API (wait ≤60s) and raise error if the response is not successful
    response.raise_for_status()

    # Return as pandas DataFrame
    return pd.read_csv(StringIO(response.text), comment="#")

# ==========================================================
# Step 1: Load catalogue of objects (.csv)
# ==========================================================
catalogue = pd.read_csv("negatives.csv")  # Columns: objectId, ra, dec
print(f"📘 Loaded catalogue with {len(catalogue)} objects")

all_lightcurves = []  # container for all results

# ==========================================================
# Step 2: Loop through catalogue and fetch lightcurves
# ==========================================================
for i, row in catalogue.head(1500).iterrows():
    object_id = str(row["objectId"])
    ra_deg = float(row["ra"])
    dec_deg = float(row["dec"])

    print(f"\n🔍 ({i+1}) Querying ZTF for {object_id} at RA={ra_deg:.4f}, Dec={dec_deg:.4f}")

    try:
        # Fetch r-band data
        df = get_ztf_lightcurve(ra_deg, dec_deg, radius_arcsec=1.2, band="r")

        if not df.empty:
            df["object"] = object_id
            all_lightcurves.append(df)
            print(f"  ✅ Got {len(df)} points for {object_id}")
        else:
            print(f"  ⚠ No data returned for {object_id}")

    except requests.exceptions.Timeout:
        print(f"  ⏱️ Timeout for {object_id} — skipping.")
        continue
    except Exception as e:
        print(f"  ❌ Failed for {object_id}: {e}")
        continue

# ==========================================================
# Step 3: Save raw merged results
# ==========================================================
if all_lightcurves:
    merged = pd.concat(all_lightcurves, ignore_index=True)
    merged.to_csv("negatives_lightcurves_raw_r.parquet", index=False)
    print(f"\n✅ Saved raw lightcurves with {len(merged)} total photometric points.")

    # ======================================================
    # Step 4: Group by object for compact, CSV-friendly format
    # ======================================================
    grouped = merged.groupby("object", sort=False, as_index=False).agg(list)

    # Convert list columns → single-line strings separated by ;
    for col in grouped.columns:
        if col != "object":
            grouped[col] = grouped[col].apply(lambda x: ";".join(map(str, x)))

    grouped.to_csv("negatives_lightcurves_grouped_r_clean.parquet", index=False)
    print("✅ Grouped output saved to 'negatives_lightcurves_grouped_r_clean.parquet'")

else:
    print("\n❌ No lightcurves were fetched from ZTF.")
