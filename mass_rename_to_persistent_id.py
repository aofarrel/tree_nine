"""
Mass rename nwks, distance matrices, and pbs from their latest (workdir) cluster IDs to their persistent IDs. We
need to do this because find_clusters.py generates clusters (and many outputs) ignorant of their persistent IDs. 
Previously we didn't bother since Microreact is the source-of-truth, but it makes sense to have better backups
that wouldn't need to be "translated" to actually be useful.
"""
VERSION="0.0.2"

# pylint: disable=trailing-whitespace,line-too-long,consider-using-tuple,useless-suppression
# pylint: disable=missing-function-docstring,broad-exception-caught,wrong-import-position
import json
import sys
from pathlib import Path
import argparse
from concurrent.futures import ThreadPoolExecutor

MAX_WORKERS = 8 


def rename_file(workdir_id, cluster_id, extension):
    old_filename = Path(f"aworkdir{workdir_id}{extension}")
    new_filename = Path(f"a{cluster_id}{extension}")
    
    if old_filename.exists():
        try:
            old_filename.rename(new_filename)
            return f"{old_filename} -> {new_filename}"
        except Exception as e:
            return f"Error renaming {old_filename} to {new_filename}: {e}"
    return f"Skipped: {old_filename} not found, not converting to {new_filename} (likely a decimated cluster)"

def main():
    parser = argparse.ArgumentParser(description="Mass rename files from find_clusters.py into their persistent ID names post process_clusters.py")
    parser.add_argument('--json', type=str, required=False, help="path to all_cluster_information JSON")
    parser.add_argument('--verbose', action='store_true', help="print all renames (recommended unless your logger is extremely slow)")
    args = parser.parse_args()

    id_map = {} # {workdir_cluster_id: cluster_id}
    print(f"Reading {args.json}...", file=sys.stderr)
    with open(args.json, 'r', encoding="utf-8") as f:
        for line in f:
            if line.strip():
                data = json.loads(line)

                # check for error cases
                if data["decimated"] is True and data["workdir_cluster_id"] is not None:
                    raise ValueError(f"Found decimated cluster with a workdir ID: {data}")
                if data["workdir_cluster_id"] is None and data["decimated"] is not True:
                    raise ValueError(f"Found null workdir_cluster_id in non-decimated cluster: {data}")
                
                # handle normal decimation case (wouldn't cause error if we didn't do this but good practice)
                if data["decimated"] is True and data["workdir_cluster_id"] is None:
                    print(f"Will skip decimated cluster with null workdir_cluster_id and persistent cluster_id {data['cluster_id']}")
                else:
                    id_map[data['workdir_cluster_id']] = data['cluster_id']

    print(f"Loaded {len(id_map)} mappings. Starting renaming...", file=sys.stderr)

    for extension in [".nwk", "_dmtrx.tsv", ".pb"]:

        # multi-threaded renaming, just for fun!
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures = [
                executor.submit(rename_file, workdir_id, cluster_id, extension) 
                for workdir_id, cluster_id in id_map.items()
            ]
            
            for i, future in enumerate(futures):
                result = future.result()
                print(result)
                if len(id_map) >= 500 and i % 500 == 0:
                    print(f"Progress: Processed {i} files...", file=sys.stderr)
    
    print("Finished renaming files", file=sys.stderr)

if __name__ == "__main__":
    main()
