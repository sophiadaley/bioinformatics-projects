import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load expression table
df = pd.read_csv("data/GSE20437_expression.tsv", sep="\t")

# Show basic info
print("Shape:", df.shape)
print("\nColumns:")
print(df.columns.tolist()[:10])  # first 10 columns
print("\nFirst 5 rows:")
print(df.head())

# Define sample groups

control_samples = [
    "GSM512539","GSM512540","GSM512541","GSM512542","GSM512543",
    "GSM512544","GSM512545","GSM512546","GSM512547","GSM512548",
    "GSM512549","GSM512550","GSM512551","GSM512552","GSM512553",
    "GSM512554","GSM512555","GSM512556",
    "GSM512575","GSM512576","GSM512577","GSM512578","GSM512579","GSM512580"
]

cancer_samples = [
    "GSM512557","GSM512558","GSM512559","GSM512560","GSM512561",
    "GSM512562","GSM512563","GSM512564","GSM512565",
    "GSM512566","GSM512567","GSM512568","GSM512569","GSM512570",
    "GSM512571","GSM512572","GSM512573","GSM512574"
]

print("Control samples:", len(control_samples))
print("Cancer samples:", len(cancer_samples))

# Check that all samples exist in the dataset

all_samples = set(df.columns)

missing_controls = [s for s in control_samples if s not in all_samples]
missing_cancer = [s for s in cancer_samples if s not in all_samples]

print("Missing control samples:", missing_controls)
print("Missing cancer samples:", missing_cancer)

# Calculate mean expression for each group

df["control_mean"] = df[control_samples].mean(axis=1)
df["cancer_mean"] = df[cancer_samples].mean(axis=1)

# Calculate difference (cancer - control)
df["diff"] = df["cancer_mean"] - df["control_mean"]

print(df[["ID_REF", "control_mean", "cancer_mean", "diff"]].head())

# Sort by biggest differences

top_up = df.sort_values("diff", ascending=False).head(10)
top_down = df.sort_values("diff", ascending=True).head(10)

print("\nTop upregulated genes (higher in cancer):")
print(top_up[["ID_REF", "diff"]])

print("\nTop downregulated genes (lower in cancer):")
print(top_down[["ID_REF", "diff"]])

# Plot distribution of differences
plt.figure()
plt.hist(df["diff"], bins=50)

plt.xlabel("Expression Difference (Cancer - Control)")
plt.ylabel("Number of Genes")
plt.title("Distribution of Gene Expression Differences")

plt.savefig("notebooks/diff_distribution.png")  # ← ADD THIS HERE
plt.show()

# Add small value to avoid division issues
df["logFC"] = np.log2((df["cancer_mean"] + 1) / (df["control_mean"] + 1))
df["abs_diff"] = abs(df["diff"])

plt.figure()

# Scatter plot
plt.scatter(df["logFC"], df["abs_diff"], alpha=0.5)

plt.xlabel("Log2 Fold Change")
plt.ylabel("Absolute Expression Difference")
plt.title("Volcano Plot (Simplified)")

# Add threshold lines (optional but nice)
plt.axvline(x=1, linestyle="--")
plt.axvline(x=-1, linestyle="--")

plt.savefig("figures/volcano_plot.png")
plt.show()