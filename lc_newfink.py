# useful python libraries
import pandas as pd
import requests
import time

# dataset
final_df = pd.read_parquet("subset_sample.parquet")

# null array
all_lightcurves = []

for object_id in final_df["objectId"].unique():

    print(f"Fetching {object_id}")

    url = "https://api.fink-portal.org/api/v1/objects"

    payload = {
        "objectId": object_id,
        "columns": "i:jd,i:magpsf,i:sigmapsf,i:fid"
    }

    try:
        r = requests.post(url, json=payload, timeout=30)
        r.raise_for_status()

        data = r.json()

        if len(data) > 0:
            df_lc = pd.DataFrame(data)

            # Keep only r-band (fid = 2)
            df_lc = df_lc[df_lc["i:fid"] == 2]

            if len(df_lc) > 0:
                df_lc["objectId"] = object_id
                all_lightcurves.append(df_lc)
                print("r-band points:", len(df_lc))
            else:
                print("No r-band data")

        time.sleep(1)

    except Exception as e:
        print("Failed:", e)

# -------------------------------------------------
# Group by objectId
# -------------------------------------------------
if len(all_lightcurves) > 0:

    merged_lc = pd.concat(all_lightcurves, ignore_index=True)

    grouped = merged_lc.groupby("objectId", as_index=False).agg({
        "i:jd": list,
        "i:magpsf": list,
        "i:sigmapsf": list
    })

    # -------------------------------------------------
    # Keep fink_class from original dataset
    # -------------------------------------------------
    labels = final_df[["objectId", "fink_class"]].drop_duplicates()

    grouped = grouped.merge(
        labels,
        on="objectId",
        how="left"
    )

    # Save
    grouped.to_parquet("fink_grouped_rband.parquet", index=False)

    print("Grouped r-band lightcurves saved with fink_class.")

else:
    print("No r-band lightcurves retrieved.")
