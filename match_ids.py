#!/usr/bin/env python3
# Licensed under the Apache License, Version 2.0, see LICENSE file for details


from argparse import ArgumentParser
from collections import defaultdict
import os
import subprocess
import sys
from typing import Dict, List, Tuple
import warnings

from helperlibs.bio.featurematch import find_features
from helperlibs.bio import seqio


REFERENCE_DB_NAME = "reference"

def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("reference", help="GBK file to use as reference")
    parser.add_argument("query", help="GBK file to query with")
    parser.add_argument("-s", "--similarity", type=float, default=90.0,
                        help="Minimum similarity (in %%) to be considered a hit. Default: %(default)s")
    parser.add_argument("-k", "--keep", action="store_true", help="Keep the intermediate files")
    parser.add_argument("-V", "--version", action="version", version="0.1.0")
    args = parser.parse_args()
    warnings.filterwarnings("ignore", ".*Partial codon.*")

    run(args.reference, args.query, args.similarity, args.keep)


def run(reference: str, query: str, similarity: float, keep: bool) -> None:
    referece_fasta = dump_fasta(reference)
    query_fasta = dump_fasta(query)

    subprocess.run(["diamond", "makedb", "--in", referece_fasta, "--db", REFERENCE_DB_NAME],
                   encoding="utf-8", capture_output=True, check=True)
    ret = subprocess.run(["diamond", "blastp", "--db", REFERENCE_DB_NAME, "--query", query_fasta,
                          "--outfmt", "6", "qseqid", "sseqid", "pident", "evalue", "bitscore"],
                         encoding="utf-8", capture_output=True)
    ret.check_returncode()
    matches = find_matches(ret.stdout, similarity)
    for qid, match_list in matches.items():
        for match in match_list:
            print(qid, *match.to_table(), sep="\t")

    if not keep:
        os.remove(referece_fasta)
        os.remove(query_fasta)
        os.remove(f"{REFERENCE_DB_NAME}.dmnd")


class Match:

    __slots__ = (
        'reference_id',
        'pident',
        'evalue',
        'bitscore',
    )

    def __init__(self, reference_id: str, pident: float, evalue: float, bitscore: float) -> None:
        self.reference_id = reference_id
        self.pident = pident
        self.evalue = evalue
        self.bitscore = bitscore

    def to_table(self) -> Tuple[str, str, str, str]:
        return self.reference_id, f"{self.pident:0.1f}", f"{self.evalue}", f"{self.bitscore:0.1f}"


def find_matches(output: str, similarity: float) -> Dict[str, List[Match]]:
    matches: Dict[str, List[Match]] = defaultdict(list)
    for line in output.splitlines():
        try:
            qid, rid, pident_str, evalue_str, bitscore_str = line.split("\t")
            pident = float(pident_str)
            evalue = float(evalue_str)
            bitscore = float(bitscore_str)
            if pident >= similarity:
                matches[qid].append(Match(rid, pident, evalue, bitscore))
        except ValueError:
            print("failed to split", line, file=sys.stderr)
            continue

    return matches


def dump_fasta(filename: str) -> str:
    records = seqio.parse(filename)
    fasta_file = derive_fasta_name(filename)

    with open(fasta_file, 'w', encoding="utf-8") as handle:
        for feature in find_features(records):
            tag = feature.feature.qualifiers['locus_tag'][0]
            sequence = str(feature.aas)
            print(f">{tag}\n{sequence}", file=handle)

    return fasta_file


def derive_fasta_name(filename: str) -> str:
    base = os.path.basename(filename)
    root, _ = os.path.splitext(base)
    return f"{root}.fa"


if __name__ == "__main__":
    main()
