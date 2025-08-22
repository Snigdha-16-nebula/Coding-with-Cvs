import requests              # used here to talk to the ZTF API
import pandas as pd          # Pandas for handling tabular data (reading/writing CSVs, merging, etc.)
from astropy.coordinates import SkyCoord   # Astropy helper to convert RA/Dec strings to decimal degrees
import astropy.units as u    # Unit system from Astropy (e.g., hours → degrees)
from io import StringIO      # Lets us treat raw text as if it were a file (for reading CSV from API response)

def get_ztf_lightcurve(ra_deg, dec_deg, radius_arcsec=2, band="g", fmt="csv"):
    """
    Function: Query ZTF’s IRSA API to download lightcurves near a sky position.
    Inputs:
        ra_deg, dec_deg  → sky coordinates in decimal degrees
        radius_arcsec    → how big a circular search region (default 2 arcsec)
        band             → filter band (e.g., 'g', 'r', 'i')
        fmt              → data format to request (default CSV)
    Output:
        Pandas DataFrame with the lightcurve data
    """
    base_url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"   # API endpoint

    radius_deg = radius_arcsec / 3600.0  # Convert arcseconds → degrees for API usage

    # Dictionary of query parameters for the API
    params = {
        "POS": f"CIRCLE {ra_deg} {dec_deg} {radius_deg}",  # Search a circle at RA, Dec, radius
        "BANDNAME": band,                                 # Which band to fetch
        "FORMAT": fmt.lower()                             # Data format in lowercase (CSV)
    }

    # Send GET request to the server with the parameters
    response = requests.get(base_url, params=params)
    response.raise_for_status()  # If something went wrong (404, timeout), throw an error

    # API response is plain text → treat it like a CSV file and load into pandas
    return pd.read_csv(StringIO(response.text), comment="#")

# -------------------------------
# Step 1: Load your catalogue file
# -------------------------------
catalogue = pd.read_csv("magnetic_cataclysmic_variables.txt")  # Read list of objects (name, RA, Dec)

all_lightcurves = []  # Empty list to collect dataframes of lightcurves

# -------------------------------
# Step 2: Only take first 20 objects from catalogue
# -------------------------------
for _, row in catalogue.head(20).iterrows():   # Loop over first 20 rows of the catalogue
    name = row["Name"]                         # Object name
    ra_str, dec_str = row["RA(J2000)"], row["DEC(J2000)"]  # RA and Dec in sexagesimal strings

    try:
        # Convert RA/Dec strings (like "12:34:56") into decimal degrees
        coord = SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg))
        ra_deg, dec_deg = coord.ra.deg, coord.dec.deg

        # Show which object we’re querying (nice progress print)
        print(f"Querying {name} at RA={ra_deg:.4f}, Dec={dec_deg:.4f}")

        # Call the function we defined earlier to get ZTF lightcurve
        df = get_ztf_lightcurve(ra_deg, dec_deg, radius_arcsec=2, band="g")

        # If data came back, attach the object name and save it
        if not df.empty:
            df["object"] = name              # Tag data with object’s name
            all_lightcurves.append(df)       # Store it in the list
            print(f"  ✅ Got {len(df)} points")   # Show success and how many points
        else:
            print(f"  ⚠ No data for {name}")      # Warn if nothing came back

    except Exception as e:
        # If something fails (bad coordinates, network, etc.), print an error
        print(f"  ❌ Failed for {name}: {e}")

# -------------------------------
# Step 3: Save everything as one CSV
# -------------------------------
if all_lightcurves:
    # Merge all the small dataframes into one big dataframe
    merged = pd.concat(all_lightcurves, ignore_index=True)

    # Save combined dataframe into a single CSV file
    merged.to_csv("lightcurves_first20.csv", index=False)

    # Success message with total number of datapoints collected
    print(f"\n✅ Saved lightcurves_first20.csv ({len(merged)} total points)")
else:
    # If no lightcurves were found for any object
    print("\n❌ No lightcurves fetched.")
