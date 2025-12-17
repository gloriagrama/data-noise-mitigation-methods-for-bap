import pandas as pd

df = pd.read_csv("three_groups_unreactive_with_vdjdb_matches.csv")

# explode rows with multiple PMIDs
df["reference.id"] = df["reference.id"].astype(str)
df = df.assign(reference_id=df["reference.id"].str.split(";")).explode("reference_id")

# count how many of your TCRs come from each study
study_counts = (
    df.groupby(["group", "reference_id"])
      .size()
      .reset_index(name="n_TCRs")
      .sort_values("n_TCRs", ascending=False)
)

print(study_counts)
