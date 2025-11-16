import pandas as pd

file_path = "Dataset/point_counts_quad.csv"

df = pd.read_csv(file_path, index_col="key")
df = df.drop_duplicates()