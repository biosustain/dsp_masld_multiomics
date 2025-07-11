import pandas as pd

# Load the TSV file
file_path = "clustered2_mice_gut_catalog_stats.tsv"
data = pd.read_csv(file_path, sep="\t")

# Assuming columns 4, 5, and 6 are 0-indexed as column indices 3, 4, 5
col4, col5, col6 = data.columns[3], data.columns[4], data.columns[5]

# Count 0s in each column
count_col4_0 = (data[col4] == 0).sum()
count_col5_0 = (data[col5] == 0).sum()
count_col6_0 = (data[col6] == 0).sum()

# Count rows with 0 in all three columns
count_all_0 = ((data[col4] == 0) & (data[col5] == 0) & (data[col6] == 0)).sum()

# Count rows with 0 in columns 4 and 5 together
count_4_and_5_0 = ((data[col4] == 0) & (data[col5] == 0)).sum()

# Count rows with 0 in column 4 but not in column 5
count_4_0_and_5_not_0 = ((data[col4] == 0) & (data[col5] != 0)).sum()

# Count rows with 0 in column 5 but not in column 4
count_5_0_and_4_not_0 = ((data[col5] == 0) & (data[col4] != 0)).sum()

# Print results
print(f"Column 4 has {count_col4_0} zeros.")
print(f"Column 5 has {count_col5_0} zeros.")
print(f"Column 6 has {count_col6_0} zeros.")
print(f"Rows with 0 in all three columns: {count_all_0}")
print(f"Rows with 0 in columns 4 and 5 together: {count_4_and_5_0}")
print(f"Rows with 0 in column 4 and not in column 5: {count_4_0_and_5_not_0}")
print(f"Rows with 0 in column 5 and not in column 4: {count_5_0_and_4_not_0}")

