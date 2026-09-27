#!/usr/bin/env python3
"""Pull the read pairs that come from a target locus (the BUB1B region) out of FASTQ files.

This is the laptop-sized alternative to aligning the whole genome. Every k-mer of the target
sequence is loaded in both orientations, so reads are looked up forward only. Low-complexity
k-mers (homopolymers and short tandem repeats) are masked, and an optional file of extra
k-mers can be excluded too (for example k-mers found elsewhere in the genome, once the
reference is available). A pair is kept when either mate has at least --min-hits target
k-mers among positions sampled every --stride bases.

The BUB1 paralog shares no 21-mer with BUB1B coding sequence (locus dossier §5), so baiting
with k >= 21 does not recruit BUB1 reads. Genome-wide repeats in introns remain the main
source of off-target pairs; filter them by alignment afterwards.

Memory holds only the target k-mers (a 125 kb locus is about 250,000 k-mers). FASTQs are
streamed. Outputs are gzip FASTQ written outside the repository.

Run: python -m pipeline.kmer_bait TARGET.fa R1.fq.gz R2.fq.gz OUT_PREFIX [--max-pairs N]
"""
import argparse
import gzip
import json
from itertools import islice
from pathlib import Path

from pipeline.common import assert_outside_repo, open_text

COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def low_complexity(kmer):
    """Homopolymer-rich or built from a repeat unit of 1-3 bases."""
    if max(kmer.count(b) for b in "ACGT") >= 0.8 * len(kmer):
        return True
    return any(kmer == (kmer[:u] * (len(kmer) // u + 1))[:len(kmer)] for u in (1, 2, 3))


def target_kmers(seqs, k, exclude=frozenset()):
    kmers, masked = set(), 0
    for seq in seqs.values():
        for i in range(len(seq) - k + 1):
            km = seq[i:i + k]
            if "N" in km:
                continue
            if low_complexity(km) or km in exclude:
                masked += 1
                continue
            kmers.add(km)
            kmers.add(revcomp(km))
    return kmers, masked


def hits(seq, kmers, k, stride):
    return sum(1 for i in range(0, len(seq) - k + 1, stride) if seq[i:i + k] in kmers)


def fastq_records(fh):
    while True:
        rec = list(islice(fh, 4))
        if len(rec) < 4:
            return
        yield rec


def bait(target_fa, r1, r2, out_prefix, k=31, min_hits=2, stride=4, max_pairs=None, exclude_path=None):
    out_prefix = Path(out_prefix)
    assert_outside_repo(out_prefix.parent / (out_prefix.name + "_R1.fq.gz"))
    exclude = frozenset()
    if exclude_path:
        exclude = frozenset(l.strip().upper() for l in open_text(exclude_path) if l.strip())
    kmers, masked = target_kmers(read_fasta(target_fa), k, exclude)
    stats = {"target_kmers": len(kmers) // 2, "masked_kmers": masked, "pairs_seen": 0, "pairs_kept": 0}
    with open_text(r1) as f1, open_text(r2) as f2, \
            gzip.open(f"{out_prefix}_R1.fq.gz", "wt") as o1, gzip.open(f"{out_prefix}_R2.fq.gz", "wt") as o2:
        for a, b in zip(fastq_records(f1), fastq_records(f2)):
            stats["pairs_seen"] += 1
            if hits(a[1].strip(), kmers, k, stride) >= min_hits or hits(b[1].strip(), kmers, k, stride) >= min_hits:
                o1.writelines(a)
                o2.writelines(b)
                stats["pairs_kept"] += 1
            if max_pairs and stats["pairs_seen"] >= max_pairs:
                break
    return stats


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("target_fa")
    ap.add_argument("r1")
    ap.add_argument("r2")
    ap.add_argument("out_prefix")
    ap.add_argument("-k", type=int, default=31)
    ap.add_argument("--min-hits", type=int, default=2)
    ap.add_argument("--stride", type=int, default=4)
    ap.add_argument("--max-pairs", type=int)
    ap.add_argument("--exclude-kmers")
    a = ap.parse_args()
    print(json.dumps(bait(a.target_fa, a.r1, a.r2, a.out_prefix, a.k, a.min_hits, a.stride, a.max_pairs,
                          a.exclude_kmers), indent=2))


if __name__ == "__main__":
    main()
