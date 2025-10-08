version 1.0

workflow sanity_check_diffs {
  input {
    File preconcatenated_diff
    Array[File] incoming_diffs
  }

  call PreconcatInfo { input: file = preconcatenated_diff }
  call IncomingInfo { input: files = incoming_diffs }

  call CompareTask {
    input:
      existing_diff_sample_level_info = PreconcatInfo.sample_level_info,
      incoming_diff_sample_level_info = IncomingInfo.sample_level_info
  }
}

task PreconcatInfo {
  input {
    File file
  }

  command <<<
    set -euox pipefail
    mkdir split_out

    # Split by sample headers (>sample_name) -- duplicates will be silently overwritten but Python will catch it
    # Yeah, this is kinda redundant with the python code but I'm in a hurry here
    awk -v outdir="split_out" '
      /^>/ {
        if (out) close(out)
        name = substr($0, 2)
        gsub(/[ \t]/, "_", name)
        out = outdir "/" name ".diff"
      }
      { print > out }
    ' "~{file}"

    # Remove terminal newlines (this is easier than doing it in awk)
    for f in split_out/*.diff
    do
      perl -pi -e 'chomp if eof' "$f"
    done

    python3 << CODE
import os, hashlib, sys

seen_samples = set()
tsv_info, all_filenames, all_samplenames, all_md5s = [], [], [], []

with open("~{file}", "r") as f:
    current_lines = []
    current_name = None
    for line in f:
        if line.startswith(">"):
            # process previous sample
            if current_name:
                file_name = f"{current_name}.diff"
                if file_name in seen_samples:
                    raise ValueError(f"Duplicate sample header found: {current_name}")
                seen_samples.add(file_name)
                
                # write this sample to a diff file, but skip the first line
                with open(os.path.join("split_out", file_name), "w") as out:
                    out.writelines(current_lines[1:])
                with open(os.path.join("split_out", file_name), "rb") as f:
                    md5 = hashlib.md5(f.read()).hexdigest()
                tsv_info.append(f"{file_name}\t{current_name}\t{md5}")
                all_filenames.append(file_name), all_samplenames.append(current_name), all_md5s.append(md5)
                if md5 in all_md5s:
                  print(f"WARNING: {current_name} has non-unique MD5sum {md5}")
            # start new sample
            current_name = line[1:].strip().replace(" ", "_")
            current_lines = [line]
        else:
            current_lines.append(line)
    # last sample
    if current_name:
        file_name = f"{current_name}.diff"
        if file_name in seen_samples:
            raise ValueError(f"Duplicate sample header found: {current_name}")
        seen_samples.add(file_name)
        with open(os.path.join("split_out", file_name), "w") as out:
            out.writelines(current_lines[:1])
        md5 = hashlib.md5(open(os.path.join("split_out", file_name), "rb").read()).hexdigest()
        tsv_info.append(f"{file_name}\t{current_name}\t{md5}")
        all_filenames.append(file_name), all_samplenames.append(current_name), all_md5s.append(md5)
with open("file_info.tsv", "w") as out:
    out.write("\n".join(tsv_info))
CODE
  >>>

  output {
    File sample_level_info = "file_info.tsv"
  }

  runtime {
    docker: "python:3.11"
  }
}

task IncomingInfo {
  input {
    Array[File] files
  }

  command <<<
    set -euox pipefail

    python3 << CODE
import os, hashlib, sys, subprocess

def print_dupes(list_with_dupes):
  seen = set()
  duplicates = set()
  for thing in list_with_dupes:
      if thing in seen:
          duplicates.add(thing)
      else:
          seen.add(thing)
  for dup in sorted(duplicates):
      print(f"  {dup}")

incoming_files = ['~{sep="','" files}']

basenames = [os.path.basename(f) for f in incoming_files]
if len(basenames) != len(set(basenames)):
    raise ValueError("Non-unique file basenames found: " + str(basenames))

tsv_info, all_filenames, all_samplenames, all_md5s = [], [], [], []
for path in sorted(incoming_files):
    file_name = os.path.basename(path)
    with open(path, "r") as f:
        sample_name = f.readline().removeprefix(">").removesuffix("\n")
    with open(path, "rb") as f:
        f.readline() # skip the first line, which has the sample name
        data = f.read()
        md5 = hashlib.md5(data).hexdigest()
    tsv_str = f"{file_name}\t{sample_name}\t{md5}"
    tsv_info.append(tsv_str)
    all_filenames.append(file_name), all_samplenames.append(sample_name), all_md5s.append(md5)

if len(set(all_filenames)) == len(all_filenames):
  print("All filenames are unique")
else:
  print("Caught non-unique filenames (basenames, to be specific):")
  print_dupes(all_filenames)
  exit(1)

if len(set(all_samplenames)) == len(all_samplenames):
  print("All sample IDs are unique")
else:
  print("Caught non-unique sample IDs:")
  print_dupes(all_samplenames)
  exit(1)

if len(set(all_md5s)) == len(all_md5s):
  print("All MD5sums are unique")
else:
  print("Caught non-unique sample IDs:")
  print_dupes(all_md5s)
  # do NOT exit early -- this is okay!!

with open("file_info.tsv", "w") as out:
    out.write("\n".join(tsv_info))
CODE
  
  >>>

  output {
    File sample_level_info = "file_info.tsv"
  }

  runtime {
    docker: "python:3.11"
  }
}

task CompareTask {
  input {
    File existing_diff_sample_level_info
    File incoming_diff_sample_level_info
  }

  command <<<
    set -euox pipefail

    python3 << CODE
import sys

def parse_array(arr):
    info = {}
    for line in arr:
        file_name, sample, md5 = line.strip().split("\t")
        info[file_name] = {"sample": sample, "md5": md5}
    return info

# nested dictionaries, like {'17RF.diff': {'sample': '17RF', 'md5': '2390118ad29a0091a'}
existing = parse_array(open("~{existing_diff_sample_level_info}").read().splitlines())
incoming = parse_array(open("~{incoming_diff_sample_level_info}").read().splitlines())

print("Every INCOMING sample is compared against EXISTING samples")

print("----- Comparing by filename -----")
for file_name in incoming:
    if file_name in existing:
        if existing[file_name]["md5"] == incoming[file_name]["md5"]:
            print(f"{file_name}: MATCH")
        else:
            print(f"{file_name}: DIFFERENT:")
            print(f"\tincoming: {incoming[file_name]['md5']}")
            print(f"\texisting: {existing[file_name]['md5']}")
    else:
       print(f"{file_name}: ABSENT")

print("\n----- Comparing by sample name -----")
# We already know sample IDs are unique within incoming, and unique within existing,
# thanks to the previous tasks asserting it
incoming_per_samp = {v["sample"]: v["md5"] for v in incoming.values()}
existing_per_samp = {v["sample"]: v["md5"] for v in existing.values()}

print(incoming_per_samp)
print(existing_per_samp)

for sample_name, md5 in incoming_per_samp.items():
    if sample_name in existing_per_samp.keys():
        if incoming_per_samp[sample_name] == existing_per_samp[sample_name]:
            print(f"{sample_name}: MATCH")
        else:
            print(f"{sample_name}: DIFFERENT:")
            print(f"\tincoming: {incoming_per_samp[sample_name]}")
            print(f"\texisting: {existing_per_samp[sample_name]}")
    else:
       print(f"{sample_name}: ABSENT")
CODE
  >>>

  runtime {
    docker: "python:3.11"
  }
}
