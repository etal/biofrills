"""Pick a subset of sequences from a larger sequence set (FASTA)."""

import argparse
import importlib.util
import itertools
import random
import sys

from Bio import SeqIO

from biofrills import alnutils, pairutils
from biofrills.cli import regular_file


def get_random(index, N):
    return random.sample(list(index), min(N, len(index)))


def get_first(index, N):
    return list(itertools.islice(index, N))


def get_by_length(infname, N, sort_reverse):
    seqlens = []
    for rec in SeqIO.parse(infname, "fasta"):
        seqlens.append((len(rec.seq), rec.id))
    seqlens.sort(reverse=sort_reverse)
    return [acc for _length, acc in seqlens[:N]]


def get_cluster(infname, N):
    """Central sequences by clustering coefficient.

    Input sequences must be aligned. Coefficients are compared after rounding
    to 12 decimal places, and equal ones are ordered by sequence ID: two
    coefficients that are equal in exact arithmetic, as for two identical
    sequences, can differ in the last place depending on the order in which
    networkx iterates over sets, and so on the hash seed.

    Where no two sequences share an aligned residue, every pairwise identity
    is 0 and no sequence is more central than another, so every coefficient
    is taken to be 0 and the sequences are picked in ID order. networkx
    divides each edge weight by the largest, which is then 0, so it cannot
    compute that case. A single sequence has no edge and takes the same
    branch; networkx handles that graph itself, with the same result.
    """
    import networkx

    G = alnutils.to_graph(infname, pairutils.identity_abs)
    if any(weight for _u, _v, weight in G.edges(data="weight")):
        coeffs = networkx.clustering(G, weight="weight")
    else:
        coeffs = dict.fromkeys(G, 0.0)
    hubs = sorted(coeffs, key=lambda node: (-round(coeffs[node], 12), node))
    return hubs[:N]


def get_kmeans(infname, N):
    """K-means clustering. (Actually, uses the k-medoids algorithm.)

    See: https://en.wikipedia.org/wiki/K-medoids

    Input sequences must be aligned. Where no two sequences share an aligned
    residue, every cost is 0, no swap lowers the total, and the medoids stay
    as get_cluster seeded them.
    """
    MAX_ITERATE = 10
    G = alnutils.to_graph(infname, pairutils.identity_abs)
    nodes = G.nodes()

    def associate(medoids):
        """Closes over: nodes, G"""
        # ENH: make assoc a Graph object? also store the distance?
        assoc = {med: set() for med in medoids}
        for otr in nodes:
            if otr in medoids:
                continue
            # The graph is complete. Of equally close medoids, max takes the
            # first, also where otr has identity 0 to every medoid.
            closest_m = max(medoids, key=lambda med: G[otr][med]["weight"])
            assoc[closest_m].add(otr)
        return assoc

    def total_cost(assoc):
        """Sum of all distances to nearest medoid. Closes over: G.

        NB: identity is a measure of similarity, not distance, so costs are
        never positive (pairwise cost = -identity).
        """
        cost = 0.0
        for med, cluster in assoc.items():
            for otr in cluster:
                cost -= G[otr][med]["weight"]
        return cost

    # 1. Select k of the n nodes as medoids
    # NB: For efficiency, take "hubs" or random?
    curr_meds = get_cluster(infname, N)

    for _i in range(MAX_ITERATE):
        # 2. Associate each "other" node to the closest medoid
        best_cost = total_cost(associate(curr_meds))
        best_meds = curr_meds

        # 3,4. For each medoid and each other node, replace the medoid with
        # each node and recalculate the cost of the new association.
        # Choose the lowest-cost association.
        for i in range(len(curr_meds)):
            for otr in nodes:
                if otr in curr_meds:
                    continue
                new_meds = curr_meds[:i] + [otr] + curr_meds[i + 1 :]
                new_cost = total_cost(associate(new_meds))
                if new_cost < best_cost:
                    best_cost = new_cost
                    best_meds = new_meds

        # 5. Repeat steps 2-4 until medoids do not change
        if best_meds == curr_meds:
            break
        curr_meds = best_meds

    return curr_meds


def positive_int(text):
    """Parse a command-line argument as an integer of at least 1."""
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError(f"must be at least 1: {text}")
    return value


def main(argv=None):
    AP = argparse.ArgumentParser(description=__doc__)
    method = AP.add_mutually_exclusive_group(required=True)
    for short, long, help_text in [
        ("-r", "--random", "Pick sequences at random."),
        ("-f", "--first", "Pick the first sequences in the file."),
        ("-l", "--longest", "Pick the longest sequences."),
        ("-s", "--shortest", "Pick the shortest sequences."),
        (
            "-c",
            "--cluster",
            "Pick the sequences with the highest clustering coefficient "
            "(aligned sequences only).",
        ),
        (
            "-k",
            "--k-means",
            "Pick the medoids of k-medoids clustering (aligned sequences only).",
        ),
    ]:
        method.add_argument(
            short,
            long,
            dest="method",
            action="store_const",
            const=long.removeprefix("--"),
            help=help_text,
        )
    AP.add_argument("nseqs", type=positive_int, help="Number of sequences to pick.")
    AP.add_argument(
        "seqfname", type=regular_file, help="Input sequences (FASTA format)."
    )
    args = AP.parse_args(argv)
    N = args.nseqs
    if args.method in ("cluster", "k-means"):
        # networkx comes with the graph extra; other methods need only Biopython
        if importlib.util.find_spec("networkx") is None:
            sys.exit(
                "seq-pick -c and -k require networkx, the graph extra of biofrills"
            )

    index = SeqIO.index(args.seqfname, "fasta")

    if args.method == "random":
        accs = get_random(index, N)
    elif args.method == "first":
        accs = get_first(index, N)
    elif args.method == "longest":
        accs = get_by_length(args.seqfname, N, True)
    elif args.method == "shortest":
        accs = get_by_length(args.seqfname, N, False)
    elif args.method == "cluster":
        accs = get_cluster(args.seqfname, N)
    else:  # k-means
        accs = get_kmeans(args.seqfname, N)

    # SeqIO.index returns each record's text from the file as bytes
    for acc in accs:
        block = index.get_raw(acc)
        sys.stdout.buffer.write(block)
