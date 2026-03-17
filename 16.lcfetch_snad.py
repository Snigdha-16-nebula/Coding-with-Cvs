#import useful python libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import requests
import time
#lightcurve extraction
# Load dataset (only first 5 rows for testing)
df_test = pd.read_parquet("my_dataset.parquet")
print("Total rows:", len(df_test))
# Group by objectId
#df_test = df_test0.groupby("objectId", as_index=False).agg(list)
#print ('total number of rows',len(df_test),'\n\n')
#groups = df_test.groupby("objectId")
#null array
all_rows = []
#r_count = 0   # counter for r-band detections
## Process objects
#i for count
for i in range(len(df_test)):
    object_id = df_test["objectId"].iloc[i]
    print(f"Processing {i+1}: {object_id}")
#taking mean ra,dec
    ra = np.mean(df_test["ra"].iloc[i])
    dec = np.mean(df_test["dec"].iloc[i])
#finkclass as first entry
    finkclass = df_test["finkclass"].iloc[i][0]
    try:
        # Query SNAD
        
        url = "https://db.ztf.snad.space/api/v3/data/latest/circle/full/json" #web

        params = {
            "ra": ra,
            "dec": dec,
            "radius_arcsec": 1.2 #small circle
        }

        r = requests.get(url, params=params)
        r.raise_for_status()

        data = r.json()

        if len(data) == 0:
            print("no data found..")    ##if no data returned
            continue
             # Extract light curve
       
        for name in data:   ##loop over objects

            lc = data[name].get("lc", [])    #"lc" key,a list of data points for the current name, defaulting to an empty list [] if it doesn't exis
            filt = data[name].get("meta", {}).get("filter", None)   #navigates nested "meta" dictionary to retrieve the "filter" value, defaulting to None if meta or filter is missing.

            for point in lc:     #loop over biggie lightcurve entries

                all_rows.append({     #adding inputs over objectids
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
final_df .to_parquet("snad_lightcurves_grouped.parquet")
grouped_df.to_parquet("snad_lightcurves_grouped.parquet")

print("Saved file: snad_lightcurves_grouped.parquet")

# Display result


print(grouped_df.head())

