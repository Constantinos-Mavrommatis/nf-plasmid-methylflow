#!/usr/bin/python3

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        prog="modkit_extract_full",
        description="Run `modkit extract full` for a single BAM",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--in-bam", required=True,
                        help="Input modBAM file (with modification tags)",
    )
    parser.add_argument("--ref", required=True,
                        help="Reference FASTA file",
    )
    parser.add_argument("--out-tsv", required=True,
                        help="Output TSV path for modkit extract full",
    )
    parser.add_argument("--cpg", action="store_true",
                        help="Restrict to CpG context (adds --cpg to modkit extract full)",
    )

    return parser.parse_args()


def run_modkit_extract_full(in_bam, ref, out_tsv, use_cpg=False):
    # Build the command as a list for safety
    cmd = ["modkit", "extract", "full", in_bam, out_tsv, "--reference", ref, "--mapped-only",]
    
    if use_cpg:
        cmd.append("--cpg")

    # This prints the command to stderr so you see it in logs
    print("Running:", " ".join(cmd), file=sys.stderr)

    # Run modkit
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise SystemExit(f"modkit extract full failed with exit code {result.returncode}")


def main():
    args = parse_args()

    # Ensure the output directory exists
    out_path = Path(args.out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    run_modkit_extract_full(
        in_bam=args.in_bam,
        ref=args.ref,
        out_tsv=args.out_tsv,
        use_cpg=args.cpg,
    )


if __name__ == "__main__":
    main()
