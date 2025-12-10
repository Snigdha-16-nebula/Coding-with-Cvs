#funstion to fetch data randomly from files
import pandas as pd

def sample_unique_from_negative_subsets(file_list, samples_per_group, object_col="object", random_state=42):
    """
    Sample unique objects from multiple groups.

    file_list         : list of parquet file paths
    samples_per_group : how many to draw from EACH group
    object_col        : column name containing object IDs
    random_state      : for reproducibility
    """

    # Set to keep track of already selected unique objects
    seen_objects = set()

    # List to store sampled dataframes
    sampled_list = []

    for f in file_list:
        df = pd.read_parquet(f)

        # Remove objects already used
        df_unique = df[~df[object_col].isin(seen_objects)]

        # How many we can sample from this file
        n = min(samples_per_group, len(df_unique))

        # Skip if no new unique objects available
        #if n == 0:
         #   continue

        # Sample n rows
        sampled = df_unique.sample(n, random_state=random_state)

        # Add to list
        sampled_list.append(sampled)

        # Mark these objects as seen
        seen_objects.update(sampled[object_col].tolist())

    # Combine all sampled dataframes
    final_df = pd.concat(sampled_list, ignore_index=True)
    return final_df

