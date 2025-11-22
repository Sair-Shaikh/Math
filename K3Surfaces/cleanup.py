import pandas as pd
import os 
# Script 1: 
quad_path = "Dataset/point_counts_quad.csv"
cubic_path = "Dataset/point_counts_cubic.csv"

df = pd.read_csv(cubic_path, index_col=0)
df = df[~df.index.duplicated(keep='last')]
df = df[df.index != ""]

df2 = pd.read_csv(quad_path, index_col=0)
df2 = df2[~df2.index.duplicated(keep='last')]
df2 = df2[df2.index != ""]



for N in range(14,21): 
    fname1 = f"Dataset/zeta_functions/point_counts_cubic_bigq_{N}.csv"
    fname2 = f"Dataset/zeta_functions/point_counts_quad_bigq_{N}.csv"

    if os.path.exists(fname1):
        dfN = pd.read_csv(fname1, index_col=0)
        df = pd.concat([df, dfN], ignore_index=True, axis="columns")
    if os.path.exists(fname2):
        dfN = pd.read_csv(fname2, index_col=0)
        df2 = pd.concat([df2, dfN], ignore_index=True, axis="columns")
df.to_csv("Dataset/zeta_functions/point_counts_cubic_alt2.csv")
df2.to_csv("Dataset/zeta_functions/point_counts_quad_alt2.csv")






# df2 = pd.read_csv(quad_path)
# df2  = df2.drop_duplicates(subset=['key'], keep='last')

# df_merged = pd.concat([df2, df], ignore_index=True)

# df.to_csv("Dataset/zeta_functions/point_counts_quad_alt.csv", index=False)
# df2.to_csv("Dataset/zeta_functions/point_counts_cubic_alt.csv", index=False)
# df_merged.to_csv("Dataset/zeta_functions/point_counts_all.csv", index=False)



# Script 3: 
# Quick script to rewrite correction array spanning multiple lines to one line for orchestrator script
# import re

# fname = "Dataset/cpp_coeffs/cube_coeffs_table_12.txt"
# out_name = "Dataset/cpp_coeffs/cube_coeffs_table_122.txt"

# with open(fname) as f:
#     text = f.read()

# # Collapse any Corrections: [ ... ] block (including newlines) to a single line
# new_text = re.sub(
#     r"Corrections:\s*\[(.*?)\]",
#     lambda m: "Corrections: [" + " ".join(m.group(1).split()) + "]",
#     text,
#     flags=re.DOTALL
# )

# with open(out_name, "w") as f:
#     f.write(new_text)