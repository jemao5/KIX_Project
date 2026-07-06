import pandas as pd

# --- 1. read the ORIGINAL pre-dedup file, aggregate Count per unique sequence ---
df_orig = pd.read_excel("raw_display_data.xlsx", engine="calamine", sheet_name="Full Library_Filtered")   # same read that worked in your filter script
print("original columns:", list(df_orig.columns))   # confirm column names first!

# clean the sequence the SAME way your filter script did: strip whitespace, uppercase
# (your filter used: "".join(str(raw_seq).split()).upper())
df_orig["seq_clean"] = df_orig["SEQUENCE"].apply(lambda s: "".join(str(s).split()).upper())

# drop blanks/invalid rows the way the filter did (optional but keeps it consistent)
df_orig = df_orig[df_orig["seq_clean"].str.fullmatch(r"[ARNDCEQGHILKMFPSTWYV]{10}")]

# groupby-sum enrichment per unique sequence
seq_count = df_orig.groupby("seq_clean")["Count"].sum().reset_index()   # adjust "COUNT" to actual col name
seq_count.columns = ["sequence", "count"]

# --- 2. read your name->sequence map (deduped) ---
df_names = pd.read_csv("tsv_outputs/full_library_out.tsv", sep="\t", header=None, names=["name", "sequence"])
# sequences in the tsv were already sanitized identically by the filter script, so they should match seq_clean

# --- 3. join enrichment onto names via sequence ---
df_out = df_names.merge(seq_count, on="sequence", how="left")
df_out["count"] = df_out["count"].fillna(0).astype(int)

# --- sanity checks ---
n_matched = (df_out["count"] > 0).sum()
print(f"{len(df_out)} names; {n_matched} matched a count; {len(df_out)-n_matched} got 0 (no match)")
print(df_out["count"].describe())

df_out[["name", "count"]].to_csv("full_library_all_metrics/enrichment_count.csv", sep="\t", index=False)
