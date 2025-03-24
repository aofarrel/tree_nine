import sys
import polars as pl

# pylint: disable=unspecified-encoding

old_malformed = sys.argv[1]
new_malformed = sys.argv[2]

old_decent = f"{old_malformed}_modified.json"
with open(old_malformed, 'r') as in_file:
	lines = in_file.readlines()
with open(old_decent, 'w') as out_file:
	out_file.write("[\n")
	for i, line in enumerate(lines):
		if i < len(lines) - 1:
		#if i < len(lines) - 2:
			out_file.write(line.strip() + ",\n")
		else:
			out_file.write(line.strip() + "\n")
	out_file.write("]\n")
old = pl.read_json(old_decent, schema_overrides={"cluster_id": pl.Utf8, "jurisdictions": pl.Utf8}).sort("cluster_id")

new_decent = f"{new_malformed}_modified.json"
with open(new_malformed, 'r') as in_file:
	lines = in_file.readlines()
with open(new_decent, 'w') as out_file:
	out_file.write("[\n")
	for i, line in enumerate(lines):
		if i < len(lines) - 1:
			out_file.write(line.strip() + ",\n")
		else:
			out_file.write(line.strip() + "\n")
	out_file.write("]\n")
new = pl.read_json(new_decent, schema_overrides={"cluster_id": pl.Utf8}).sort("cluster_id")
inner_join = old.join(new, on="cluster_id", suffix="_new", how="inner").sort("cluster_id")
different_dates = inner_join.filter(pl.col("first_found") != pl.col("first_found_new"))
if not different_dates.is_empty():
	print("\n⚠️ Date mismatch found in 'date_first_found':")
	for row in different_dates.iter_rows(named=True):
		cluster_id = row["cluster_id"]
		date1, date2 = row["first_found"], row["first_found_new"]
		print(f"  {cluster_id}: {date1} → {date2}")


for col in old.columns:
	if col == "cluster_id":
		continue

	col_old, col_new = col, f"{col}_new"
	#with pl.Config(fmt_str_lengths=500, fmt_table_cell_list_len=500):
	#	print(old.select(["cluster_id", col]))
	#	print(new.select(["cluster_id", col]))
	#	print(inner_join.select(["cluster_id", col_old, col_new]))
	#print("----------------")

	different_rows = inner_join.filter(pl.col(col_old) != pl.col(col_new)) # TODO: this seems to ignore nulls
	if different_rows.is_empty() and inner_join.schema[col_old] != pl.List:
		continue
	print(f"Column '{col}' has differences:")
	for row in different_rows.iter_rows(named=True):
		cluster_id = row["cluster_id"]
		val1, val2 = row[col_old], row[col_new]

		if isinstance(val1, list) and isinstance(val2, list):
			# Find added/removed elements for list columns
			set1, set2 = set(val1), set(val2)
			added = set2 - set1
			removed = set1 - set2

			print(f"  {cluster_id}: Added {len(list(added))} samples: {list(added)} \n Removed {list(removed)}")
		else:
			print(f"  {cluster_id}: {val1} → {val2}")
