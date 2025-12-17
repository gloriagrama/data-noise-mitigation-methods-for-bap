import pandas as pd
import numpy as np

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------
df = pd.read_csv("01_05_2025_TCRvdb.csv")

# ------------------------------------------------------------
# Restrict to epitopes used in Fig 3C
# ------------------------------------------------------------
target_epitopes = ["YLQPRTFLL", "GLCTLVAML"]
df = df[df["epitope_aa"].isin(target_epitopes)].copy()

# ------------------------------------------------------------
# Define functional validation (per paper)
# A TCR is considered validated ONLY if:
#   - log2FoldChange > 0
#   - padj < 1e-5
# ------------------------------------------------------------
df["validated"] = (
    (df["log2FoldChange"] > 0) &
    (df["padj"] < 1e-5)
)

# ------------------------------------------------------------
# Keep ONLY unreactive TCRs (grey points in Fig 3C)
# ------------------------------------------------------------
unreactive = df[~df["validated"]].copy()

# ------------------------------------------------------------
# GROUP 1: LOW — clearly non-reactive
#   Strongly negative log2FoldChange
# ------------------------------------------------------------
group_low = (
    unreactive
    .sort_values("log2FoldChange", ascending=True)
    .head(5)
    .copy()
)
group_low["group"] = "low_nonreactive"

# ------------------------------------------------------------
# GROUP 2: MID — unreactive, weak / noisy
#   log2FoldChange closest to 0
# ------------------------------------------------------------
group_mid = (
    unreactive
    .assign(abs_fc=lambda x: x["log2FoldChange"].abs())
    .sort_values("abs_fc", ascending=True)
    .head(5)
    .copy()
)
group_mid["group"] = "mid_unreactive"

# ------------------------------------------------------------
# GROUP 3: HIGH-BUT-FAILED — near misses
#   Positive log2FoldChange, but padj too high
# ------------------------------------------------------------
group_high_failed = (
    unreactive[unreactive["log2FoldChange"] > 0]
    .sort_values(
        ["log2FoldChange", "padj"],
        ascending=[False, True]
    )
    .head(5)
    .copy()
)
group_high_failed["group"] = "high_failed"

# ------------------------------------------------------------
# Combine all groups
# ------------------------------------------------------------
selected = pd.concat(
    [group_low, group_mid, group_high_failed],
    ignore_index=True
)

# ------------------------------------------------------------
# Save output
# ------------------------------------------------------------
selected.to_csv(
    "three_groups_unreactive_15_TCRs.csv",
    index=False
)

# ------------------------------------------------------------
# Display summary
# ------------------------------------------------------------
print(
    selected[
        [
            "group",
            "name",
            "epitope_aa",
            "log2FoldChange",
            "padj",
            "vdjdb_score",
            "cdr3_alpha_aa",
            "cdr3_beta_aa"
        ]
    ]
)
