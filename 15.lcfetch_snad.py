import requests
import pandas as pd
import time

# Load dataset
df = pd.read_parquet("subset_sample_2.parquet")

print("Total rows:", len(df))

# Group by objectId

groups = df.groupby("objectId")

all_rows = []

# Process all objects

for i, (object_id, group) in enumerate(groups):

    print(f"Processing {i+1}: {object_id}")

    # Extract metadata
    ra = group["ra"].iloc[0]
    dec = group["dec"].iloc[0]
    finkclass = group["finkclass"].iloc[0]

    try:

       
        # Query SNAD

        url = "https://db.ztf.snad.space/api/v3/data/latest/circle/full/json"

        params = {
            "ra": ra,
            "dec": dec,
            "radius_arcsec": 1.2
        }

        r = requests.get(url, params=params)
        r.raise_for_status()

        data = r.json()

        if len(data) == 0:
            continue

    
        # Extract light curve
       

        for name in data:

            lc = data[name].get("lc", [])
            filt = data[name].get("meta", {}).get("filter", None)

            for point in lc:

                all_rows.append({
                    "objectId": object_id,
                    "finkclass": finkclass,
                    "oid": name,
                    "mjd": point.get("mjd"),
                    "mag": point.get("mag"),
                    "magerr": point.get("magerr"),
                    "filter": filt
                })

    except Exception as e:
        print("Error fetching", object_id, e)

    # avoid hitting API too fast
    time.sleep(0.5)


# Create dataframe
final_df = pd.DataFrame(all_rows)

print("Total photometry rows collected:", len(final_df))

# group
grouped_df = final_df.groupby("objectId").agg({
    "finkclass": "first",
    "oid": "first",
    "filter": list,
    "mjd": list,
    "mag": list,
    "magerr": list
}).reset_index()

print("Objects after grouping:", len(grouped_df))

# Save output
grouped_df.to_parquet("snad_lightcurves_grouped_2.parquet")

print("Saved file: snad_lightcurves_grouped_2.parquet")

# Display result



print(grouped_df.head())
