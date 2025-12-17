import pandas as pd

df = pd.read_csv("three_groups_unreactive_with_vdjdb_matches.csv")

# drop rows without VDJdb references
df = df.dropna(subset=["reference.id"])

# split multiple PMIDs
df = df.assign(
    reference_id=df["reference.id"].str.split(";")
).explode("reference_id")

# count TCRs per study per group
study_counts = (
    df.groupby(["group", "reference_id"])
      .size()
      .reset_index(name="n_TCRs")
      .sort_values("n_TCRs", ascending=False)
)

# save output
study_counts.to_csv(
    "study_counts_by_group.csv",
    index=False
)

print(study_counts)
