import pandas as pd
from scipy.stats import shapiro, mannwhitneyu, ttest_ind

# Load the updated Ks file with tag information
file_path = "Ks_with_tags.txt"
ks_df = pd.read_csv(file_path, sep="\t")

# Separate the Ks values into tag and nontag groups
tag_ks = ks_df[ks_df["Tag_Status"] == "tag"]["Ks"]
nontag_ks = ks_df[ks_df["Tag_Status"] == "nontag"]["Ks"]

# Step 1: Test for normality
tag_normality = shapiro(tag_ks)
nontag_normality = shapiro(nontag_ks)

print("Shapiro-Wilk Test for Normality:")
print(f"Tag group: W={tag_normality.statistic}, p-value={tag_normality.pvalue}")
print(f"Nontag group: W={nontag_normality.statistic}, p-value={nontag_normality.pvalue}")

# Step 2: Choose the test
if tag_normality.pvalue < 0.05 or nontag_normality.pvalue < 0.05:
    # Use Mann-Whitney U test if data is not normal
    stat, p_value = mannwhitneyu(tag_ks, nontag_ks, alternative="two-sided")
    print("\nMann-Whitney U Test:")
    print(f"U-statistic={stat}, p-value={p_value}")
else:
    # Use Welch's t-test if data is normal
    t_stat_welch, p_value_welch = ttest_ind(tag_ks, nontag_ks, equal_var=False)
    print("\nT-Test:")
    print(f"Welch's t-test result: t-statistic={t_stat_welch}, p-value={p_value_welch}")

import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Prepare data for plotting
ks_df["Tag_Status"] = ks_df["Tag_Status"].replace({"tag": "Tag", "nontag": "Non-Tag"})

# Plot
plt.figure(figsize=(10, 6))
sns.boxplot(data=ks_df, x="Tag_Status", y="Ks")
plt.title("Comparison of Ks Values Between Tag and Non-Tag Groups")
plt.xlabel("Group")
plt.ylabel("Ks Value")
plt.text(0.5, max(ks_df["Ks"]), f'Mann-Whitney U p-value = {p_value:.2e}', 
         horizontalalignment='center', verticalalignment='top', fontsize=12, color='red')
plt.savefig("ks_boxplot.png")

plt.figure(figsize=(10, 6))
sns.kdeplot(tag_ks, label="Tag", fill=True, bw_adjust=0.5)
sns.kdeplot(nontag_ks, label="Non-Tag", fill=True, bw_adjust=0.5)
plt.title("Density Plot of Ks Values for Tag and Non-Tag Groups")
plt.xlabel("Ks Value")
plt.ylabel("Density")
plt.legend()
plt.savefig("ks_density.png")
