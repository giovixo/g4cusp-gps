#!/usr/bin/env python3
"""Merge the per-thread CSV files written by cusp-activation into one file per run.

In multi-thread mode every worker thread writes its own file:

    scorefile_run<N>_nt_Events_t<thread>.csv

This script concatenates, for each run N, all the thread files into

    scorefile_run<N>.csv

keeping the Geant4 CSV header (the "#..." lines) of the input files, and sorting the rows by EventID.
The input files are left untouched unless --delete is given.

Usage:
    python3 tools/merge_scorefiles.py [directory] [-o OUTDIR] [--run N] [--delete]
"""

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

THREAD_FILE = re.compile(r"^(?P<prefix>.+)_run(?P<run>\d+)_nt_(?P<ntuple>\w+?)_t(?P<thread>\d+)\.csv$")


def read_csv(path):
    """Return (header lines, data rows) of a Geant4 CSV ntuple file."""
    header, rows = [], []
    with open(path, newline="") as f:
        for line in f:
            if line.startswith("#"):
                header.append(line)
            elif line.strip():
                rows.append(line)
    return header, rows


def event_id(row):
    try:
        return int(next(csv.reader([row]))[0])
    except (ValueError, IndexError):
        return sys.maxsize


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("directory", nargs="?", default=".", help="directory with the thread files (default: .)")
    parser.add_argument("-o", "--outdir", help="output directory (default: same as input)")
    parser.add_argument("--run", type=int, help="merge only this run number")
    parser.add_argument("--delete", action="store_true", help="delete the thread files after merging")
    args = parser.parse_args()

    indir = Path(args.directory)
    if not indir.is_dir():
        print(f"ERROR: {indir} is not a directory", file=sys.stderr)
        return 1
    outdir = Path(args.outdir) if args.outdir else indir
    outdir.mkdir(parents=True, exist_ok=True)

    # Group the thread files by (prefix, run, ntuple)
    groups = defaultdict(list)
    for path in indir.iterdir():
        m = THREAD_FILE.match(path.name)
        if m and (args.run is None or int(m["run"]) == args.run):
            groups[(m["prefix"], int(m["run"]), m["ntuple"])].append((int(m["thread"]), path))

    if not groups:
        print(f"No thread files found in {indir}", file=sys.stderr)
        return 1

    for (prefix, run, ntuple), files in sorted(groups.items()):
        files.sort()
        header, rows = None, []
        for _, path in files:
            h, r = read_csv(path)
            if header is None:
                header = h
            elif h != header:
                print(f"ERROR: {path} has a different header from {files[0][1]}; run {run} not merged",
                      file=sys.stderr)
                return 1
            rows.extend(r)
        rows.sort(key=event_id)

        out = outdir / f"{prefix}_run{run}.csv"
        with open(out, "w", newline="") as f:
            f.writelines(header)
            f.writelines(rows)
        print(f"{out}: {len(rows)} rows from {len(files)} thread file(s)")

        if args.delete:
            for _, path in files:
                path.unlink()

    return 0


if __name__ == "__main__":
    sys.exit(main())
