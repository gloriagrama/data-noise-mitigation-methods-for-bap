import pandas as pd

# ------------------------------------------------------------
# Load your curated unreactive TCRs (15 rows)
# ------------------------------------------------------------
selected = pd.read_csv("three_groups_unreactive_15_TCRs.csv")

# ------------------------------------------------------------
# Load full VDJdb database
# ------------------------------------------------------------
vdjdb = pd.read_csv(
    "vdjdb.txt",
    sep="\t",
    low_memory=False
)

# ------------------------------------------------------------
# Create matching keys
# Match on CDR3β + epitope (maximizes recall)
# ------------------------------------------------------------
selected["match_key"] = (
    selected["cdr3_beta_aa"].astype(str) + "|" +
    selected["epitope_aa"].astype(str)
)

vdjdb["match_key"] = (
    vdjdb["cdr3"].astype(str) + "|" +
    vdjdb["antigen.epitope"].astype(str)
)

# ------------------------------------------------------------
# Merge selected TCRs with VDJdb
# ------------------------------------------------------------
merged = selected.merge(
    vdjdb,
    on="match_key",
    how="left",
    suffixes=("_screen", "_vdjdb")
)

# ------------------------------------------------------------
# Flag alpha-chain agreement (when available)
# ------------------------------------------------------------
# Save full reconciliation table
# ------------------------------------------------------------
merged.to_csv(
    "three_groups_unreactive_with_vdjdb_matches.csv",
    index=False
)

# ------------------------------------------------------------
# Print a concise summary to screen
# ------------------------------------------------------------
summary_cols = [
    "group",
    "name",
    "epitope_aa",
    "cdr3_alpha_aa",
    "cdr3_beta_aa",
    "vdjdb.score",
    "mhc.a",
    "reference.id",
    "method"
]

print("\n=== VDJdb reconciliation summary ===\n")
print(merged[summary_cols])
