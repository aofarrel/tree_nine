import sys
import logging
from datetime import datetime 
import polars as pl
today = datetime.utcnow().date()
pl.Config.set_tbl_rows(-1)
pl.Config.set_tbl_cols(-1)
pl.Config.set_tbl_width_chars(150)
pl.Config.set_fmt_str_lengths(5000)
pl.Config.set_fmt_table_cell_list_len(5000)

# pylint: disable=use-list-literal

df = pl.read_ndjson(sys.argv[1], ignore_errors=True)

change_report = list()

for row in df.iter_rows(named=True):
    logging.debug("Checking %s", row['cluster_id'])
    try:
        what_is = set(row["sample_id"])
    except TypeError:
        what_is = set()
    try:
        what_was = set(row["sample_id_previously"])
    except TypeError:
        what_was = set()
    gained = list(what_is - what_was)
    lost = list(what_was - what_is)
    logging.debug("what is: %s", what_is)
    logging.debug("what was: %s", what_was)
    logging.debug("gained: %s", gained)
    logging.debug("lost: %s", lost)
    change_report.append({"cluster": f"{row['cluster_id']}@{row['cluster_distance']}", 
    	"gained": gained, "lost": lost, "kept": list(what_is.intersection(what_was)),
    	"microreact_url": row['microreact_url']})
change_report_df = pl.DataFrame(change_report).with_columns([
    pl.when(pl.col('gained').list.len() == 0).then(None).otherwise(pl.col('gained')).alias("gained"),
    pl.when(pl.col('lost').list.len() == 0).then(None).otherwise(pl.col('lost')).alias("lost"),
    pl.when(pl.col('kept').list.len() == 0).then(None).otherwise(pl.col('kept')).alias("kept"),
    pl.when(pl.col('microreact_url').is_null()).then(None).otherwise(pl.col('microreact_url')).alias("microreact_url"),
])
change_report_df = change_report_df.with_columns([
	pl.when(pl.col('gained').is_not_null()).then(pl.col('gained').list.len()).otherwise(pl.lit(0)).alias("n_gained"),
	pl.when(pl.col('lost').is_not_null()).then(pl.col('lost').list.len()).otherwise(pl.lit(0)).alias("n_lost"),
	pl.when(pl.col('kept').is_not_null()).then(pl.col('kept').list.len()).otherwise(pl.lit(0)).alias("n_kept"),
])
change_report_df = change_report_df.with_columns(n_now=pl.col('n_gained')-pl.col('n_lost')+pl.col('n_kept'))

#print("Finished. Here's how clusters have changed:")
#print(change_report_df)
#change_report_df.write_ndjson(f'change_report{today.isoformat()}.json')

print("Existing clusters that lost samples (note: it's possible to gain and lose)")
lost_samples = change_report_df.filter(pl.col("lost").is_not_null()).select(['cluster', 'n_gained', 'n_lost', 'n_kept', 'microreact_url', 'lost'])
print(lost_samples)

print("Existing clusters that gained samples (note: it's possible to gain and lose)")
gained_samples = change_report_df.filter((pl.col("gained").is_not_null().and_(pl.col("kept").is_not_null()))).select(['cluster', 'n_gained', 'n_lost', 'n_kept', 'n_now', 'microreact_url', 'gained'])
print(gained_samples)

print("Brand new clusters")
new = change_report_df.filter((pl.col("gained").is_not_null().and_(pl.col("kept").is_null()))).select(['cluster', 'n_gained', 'microreact_url', 'gained'])
print(new)

print("Decimated clusters")
decimated = change_report_df.filter((pl.col("lost").is_null()).and_(pl.col("gained").is_null()).and_(pl.col("kept").is_null())).select(['cluster', 'n_gained', 'n_lost', 'n_kept', 'n_now', 'microreact_url'])
print(decimated)

print("Unchanged clusters")
unchanged = change_report_df.filter((pl.col("lost").is_null()).and_(pl.col("gained").is_null()).and_(pl.col("kept").is_not_null())).select(['cluster', 'n_now', 'microreact_url'])
print(unchanged)
