#!/usr/bin/env python3
"""Write the inputs under tests/data/inputs/ for the reference-output capture.

Three kinds of input are written:

* Real protein alignments, derived from the MAPGAPS CMA alignments that
  accompany cladecompare and biocma.  Each CMA file is reduced to its
  consensus (match-state) columns by the same rule as
  ``cladecomparelib.core.cma_blocks``: lowercase insert residues are dropped,
  leaving uppercase residues and '-' deletions of equal length.  This is the
  form in which cladecompare passes alignments to biofrills.
* Real unaligned protein sequences and a CMA alignment, copied unchanged,
  and a subset of that alignment, for the scripts that operate on sequence
  sets.
* Small synthetic alignments and sequence sets that exercise edge cases
  (ties, mixed case, all-gap and terminal-gap columns, a single sequence,
  'X' residues, duplicate records, nucleotide input for ORF translation).
  The random ones use a fixed seed, so rerunning this script reproduces
  them byte for byte.

Usage, from the repository root, with cladecompare and biocma checked out
beside this repository:

    python3 tests/data/make_inputs.py [--cladecompare DIR] [--biocma DIR]
"""

import argparse
import random
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
OUT = HERE / "inputs"


def read_fasta(path):
    """Return a list of (header, sequence) pairs; the header excludes '>'."""
    records = []
    header = None
    chunks = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header = line[1:]
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def write_fasta(name, records, width=None):
    """Write (header, sequence) pairs, unwrapped unless a width is given."""
    lines = []
    for header, seq in records:
        lines.append(">" + header)
        if width:
            lines.extend(seq[i : i + width] for i in range(0, len(seq), width))
        else:
            lines.append(seq)
    (OUT / name).write_text("\n".join(lines) + "\n")


def cma_blocks(path):
    """Consensus columns of a CMA alignment, as cladecompare reads it."""
    records = []
    lines = iter(Path(path).read_text().splitlines())
    for line in lines:
        if line.startswith(">"):
            acc = line.split(None, 1)[0][1:]
            seq = next(lines).strip()
            seq = "".join(c for c in seq[3:-4] if not c.islower())
            records.append((acc, seq))
    lengths = {len(seq) for _acc, seq in records}
    assert len(lengths) == 1, f"{path}: unequal row lengths {sorted(lengths)}"
    return records


def cma_subset(path, n_keep):
    """The first n_keep sequences of a single-block CMA file, as text.

    The block header's sequence count is updated, so that the result is a
    valid CMA file sharing identifiers with the original, for a seq-diff
    case that subtracts something.
    """
    lines = Path(path).read_text().splitlines()
    starts = [i for i, line in enumerate(lines) if line.startswith("$")]
    footer = next(i for i, line in enumerate(lines) if line.startswith("_0]"))
    head = lines[: starts[0]]
    name, rest = head[0].split("=", 1)
    total = len(starts)
    assert f"({total})" in rest, f"{path}: header does not count {total}"
    head[0] = name + "=" + rest.replace(f"({total})", f"({n_keep})", 1)
    kept = lines[starts[0] : starts[n_keep]]
    return "\n".join(head + kept + lines[footer:]) + "\n"


# Hand-written edge cases.  Each is one alignment.

EDGE = [
    # All-gap first and last columns, an all-'.' column, N- and C-terminal
    # gaps, a column that is gap in most rows, and a description on one row.
    ("s1 first sequence", "--MKV-LA.YRGH--"),
    ("s2", "--MKVILA.FRGH--"),
    ("s3", "---KVIL-.YKGHQ-"),
    ("s4", "-AMRV-LA.YRGN--"),
    ("s5", "---KV-LS.YRG---"),
]

TIES = [
    # Every column is either a 2-2 split or four different residues, so the
    # consensus is decided by tie-breaking throughout.
    ("t1", "AAKLCDEW"),
    ("t2", "AGKLDEFW"),
    ("t3", "GAMICFGY"),
    ("t4", "GGMIDGHY"),
]

MIXEDCASE = [
    # Counting columns from 1: column 3 is lowercase in every row (excluded
    # from the consensus); column 5 mixes cases; 'w' in column 6 occurs only
    # in lowercase, so its uppercase form is absent from the background
    # frequencies.
    ("m1", "MKaVLwGH"),
    ("m2", "MKaVlwGH"),
    ("m3", "MRa-LwGN"),
    ("m4", "MKa-LYGH"),
]

SINGLE = [
    ("only", "MKV-LAYRGH"),
]

XRES = [
    # 'X' residues, for the gap_chars='-.X' call made by cladecompare.
    ("x1", "MKXVLAYXGH"),
    ("x2", "MKXVLAFRGH"),
    ("x3", "MRXV-AYRGX"),
    ("x4", "XKXVLSYRGH"),
]


def tie_batch(rng, n_alns=200):
    """Many small alignments over a small alphabet, so ties are common.

    Record identifiers are ``t<alignment>_s<row>``; the capture script groups
    records into alignments by the prefix before '_'.
    """
    alphabet = "AACDEG-"
    records = []
    for aln_idx in range(n_alns):
        n_rows = rng.randint(3, 7)
        width = rng.randint(6, 14)
        for row in range(n_rows):
            seq = "".join(rng.choice(alphabet) for _ in range(width))
            records.append((f"t{aln_idx:03d}_s{row}", seq))
    return records


def dna_set(rng):
    """Nucleotide records for 6-frame ORF translation."""
    stops = {"TAA", "TAG", "TGA"}
    codons = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"]
    sense = [c for c in codons if c not in stops]

    def random_dna(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    records = []
    # Random sequence: stop codons every ~21 codons, so short ORFs only.
    records.append(("rand1 random", random_dna(600)))
    records.append(("rand2", random_dna(451)))
    # One long open reading frame on the forward strand, flanked by random
    # sequence.
    orf = "ATG" + "".join(rng.choice(sense) for _ in range(120)) + "TAA"
    records.append(("orf1 long forward ORF", random_dna(31) + orf + random_dna(40)))
    # The same kind of ORF on the reverse strand, with ambiguous bases.
    orf = "ATG" + "".join(rng.choice(sense) for _ in range(90)) + "TGA"
    revcomp = orf[::-1].translate(str.maketrans("ACGT", "TGCA"))
    records.append(("orf2 reverse ORF with Ns", random_dna(20) + "NN" + revcomp + random_dna(25)))
    return records


def dupes_set(real):
    """Records with each kind of duplicate that seq-purge distinguishes.

    ``real`` is a list of (header, sequence) pairs of real proteins.
    """
    (h1, q1), (h2, q2), (h3, q3), (h4, q4) = real[:4]
    id1 = h1.split()[0]
    id2 = h2.split()[0]
    return [
        (h1, q1),
        (h2, q2),
        (h3, q3),
        (h4, q4),
        # Exact duplicate: same identifier, description and sequence.
        (h1, q1),
        # Same identifier and sequence, different description.
        (id2 + " another description", q2),
        # Same identifier, different sequence.
        (id1, q3[:120]),
        # Same sequence, different identifier.
        ("copy_of_3", q3),
        # Contained in record 4.
        ("fragment_of_4", q4[30:150]),
        # Contains 'X'; each X-separated segment is contained in record 1, but
        # the two lie 20 residues apart there, so the current seq-purge -c
        # keeps it, where the earlier segment-wise rule purged it.
        ("x_fragment_of_1", q1[10:60] + "X" + q1[80:140]),
        # Contains 'X' and one segment is novel.
        ("x_novel", q1[10:60] + "X" + "WWWWWWWWWW"),
    ]


def intersect_query(ref):
    """Query identifiers for seq-intersect against ``ref``.

    Includes plain matches, HMMer-style '/start-end' suffixes (matched only
    with -r) and identifiers absent from the reference.
    """
    ids = [h.split()[0] for h, _s in ref]
    return [
        (ids[0], "M"),
        (ids[5], "M"),
        (ids[2] + "/12-140", "M"),
        ("not_in_reference", "M"),
        (ids[9] + "/1-50", "M"),
        (ids[30], "M"),
    ]


def diff_query(ref):
    """Records for seq-diff against ``ref``, one per kind of key match.

    Each of the -i, -h, -b and -a keys matches a different subset of these
    records, so each of those seq-diff cases subtracts different records.
    """
    (h0, s0), (h1, s1), (h2, s2), (_h3, s3) = ref[:4]
    return [
        # Same identifier, description and sequence.
        (h0, s0),
        # Same identifier and description, different sequence.
        (h1, s1[:-5]),
        # Same identifier and sequence, different description.
        (h2.split()[0] + " a different description", s2),
        # Same sequence, different identifier.  Only -s matches this record,
        # and -s cannot complete against cdc2-domain.fasta, which contains a
        # duplicate sequence; the record discriminates only if a port lets
        # -s complete.
        ("renamed_copy", s3),
        # Nothing in common.
        ("novel", "MKVLAAGGHW"),
    ]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cladecompare", default=REPO.parent / "cladecompare", type=Path)
    ap.add_argument("--biocma", default=REPO.parent / "biocma", type=Path)
    args = ap.parse_args()
    # Every file here is generated, so start empty: a renamed input must not
    # leave its predecessor behind for the capture to pick up.
    shutil.rmtree(OUT, ignore_errors=True)
    OUT.mkdir()

    cc_test = args.cladecompare / "test"
    write_fasta("cdc2-domain.aln.fasta", cma_blocks(cc_test / "cdc2-domain.fasta_aln.cma"))
    write_fasta("scttlre-domain.aln.fasta", cma_blocks(cc_test / "scttlre-domain.fasta_aln.cma"))
    write_fasta("fikk.aln.fasta", cma_blocks(args.biocma / "test" / "fikk.cma"))

    shutil.copyfile(cc_test / "cdc2-domain.fasta", OUT / "cdc2-domain.fasta")
    shutil.copyfile(args.biocma / "test" / "fikk.cma", OUT / "fikk.cma")
    (OUT / "fikk-subset.cma").write_text(cma_subset(args.biocma / "test" / "fikk.cma", 15))
    cdc2 = read_fasta(cc_test / "cdc2-domain.fasta")

    write_fasta("synth-edge.aln.fasta", EDGE)
    write_fasta("synth-ties.aln.fasta", TIES)
    write_fasta("synth-mixedcase.aln.fasta", MIXEDCASE)
    write_fasta("synth-single.aln.fasta", SINGLE)
    write_fasta("synth-x.aln.fasta", XRES)

    # Each generator has its own stream, so changing one leaves the others'
    # output, and the fixtures captured from it, unchanged.
    write_fasta("tie-batch.fasta", tie_batch(random.Random(20261023)))
    write_fasta("dna.fasta", dna_set(random.Random(20261024)), width=60)
    write_fasta("dupes.fasta", dupes_set(read_fasta(cc_test / "scttlre-domain.fasta")), width=60)
    write_fasta("intersect-query.fasta", intersect_query(cdc2))
    write_fasta("diff-query.fasta", diff_query(cdc2), width=60)


if __name__ == "__main__":
    main()
