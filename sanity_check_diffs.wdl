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
    set -euo pipefail

    mkdir split_out

    # Split by sample headers (>sample_name)
    awk -v outdir="split_out" '
      /^>/ {
        if (out) close(out)
        name = substr($0, 2)
        gsub(/[ \t]/, "_", name)
        out = outdir "/" name ".diff"
      }
      { print > out }
    ' ~{file}

    python3 << CODE
import os, hashlib, sys
out_lines = []
for fname in sorted(os.listdir("split_out")):
    path = os.path.join("split_out", fname)
    sample_name = os.path.splitext(fname)[0]
    with open(path, "rb") as f:
        md5 = hashlib.md5(f.read()).hexdigest()
    line = f"{fname}\t{sample_name}\t{md5}"
    print(line)
    out_lines.append(line)
with open("file_info.tsv", "w") as out:
    out.write("\n".join(out_lines))
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
    set -euo pipefail

    python3 << CODE
import os, hashlib, sys

lines = []

for path in sorted(['~{sep="','" files}']):
    fname = os.path.basename(path)
    sample_name = os.path.splitext(fname)[0]
    with open(path, "rb") as f:
        md5 = hashlib.md5(f.read()).hexdigest()
    line = f"{fname}\t{sample_name}\t{md5}"
    print(line)
    lines.append(line)
with open("file_info.tsv", "w") as out:
    out.write("\n".join(lines))
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
    set -euo pipefail

    python3 << CODE
import sys

def parse_array(arr):
    info = {}
    for line in arr:
        fname, sample, md5 = line.strip().split("\t")
        info[fname] = {"sample": sample, "md5": md5}
    return info

existing = parse_array(open("~{existing_diff_sample_level_info}").read().splitlines())
incoming = parse_array(open("~{incoming_diff_sample_level_info}").read().splitlines())

print("Every INCOMING sample is compared against EXISTING samples")

print("----- Comparing by filename -----")
for fname in incoming:
    if fname in existing:
        if existing[fname]["md5"] == incoming[fname]["md5"]:
            print(f"{fname}: MATCH")
        else:
            print(f"{fname}: DIFFERENT")
    else:
       print(f"{fname}: ABSENT")

print("\n----- Comparing by sample name -----")
samples2 = {v["sample"]: v["md5"] for v in incoming.values()}
for fname, info1 in existing.items():
    sname = info1["sample"]
    if sname in samples2:
        if info1["md5"] == samples2[sname]:
            print(f"{sname}: MATCH")
        else:
            print(f"{sname}: DIFFERENT")
    else:
      print(f"{sname}: ABSENT")
CODE
  >>>

  runtime {
    docker: "python:3.11"
  }
}
