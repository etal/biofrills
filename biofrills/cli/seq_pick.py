"""Pick a subset of k sequences from a larger sequence set (FASTA).

Usage:
    seq-pick [-r|-f|-l|-s|-c|-k] <nseqs> <seqfname>

Methods for selection:

    - random (-r)
    - first (-f)
    - longest (-l)
    - shortest (-s)

For alignments only:

    - clustering coefficient (-c)
    - k-means (-k)
"""

import random
import sys

from Bio import SeqIO

from biofrills import alnutils, pairutils
from biofrills.cli import help_requested


def get_random(index, N):
    return random.sample(index, N)


def get_first(infname, N):
    accs = []
    with open(infname) as infile:
        for line in infile:
            if line.startswith(">"):
                accs.append(line.split(None, 1)[0][1:])
                if len(accs) == N:
                    return accs


def get_by_length(infname, N, sort_reverse):
    seqlens = []
    for rec in SeqIO.parse(infname, "fasta"):
        seqlens.append((len(rec.seq), rec.id))
    seqlens.sort(reverse=sort_reverse)
    return [acc for _length, acc in seqlens[:N]]


def get_cluster(infname, N):
    """Central sequences by clustering coefficient.

    Input sequences must be aligned.
    """
    import networkx

    G = alnutils.to_graph(infname, pairutils.identity_abs)
    coeffs = networkx.clustering(G, weight="weight")
    hubs = [k for k, _ in sorted(coeffs.items(), key=lambda kv: kv[1], reverse=True)]
    return hubs[:N]


def get_kmeans(infname, N):
    """K-means clustering. (Actually, uses the k-medoids algorithm.)

    See: https://en.wikipedia.org/wiki/K-medoids

    Input sequences must be aligned.
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
            closest_m = None
            closest_m_ident = 0
            for med in medoids:
                this_m_ident = G[otr][med]["weight"]
                if this_m_ident > closest_m_ident:
                    closest_m = med
                    closest_m_ident = this_m_ident
            assoc[closest_m].add(otr)
        return assoc

    def total_cost(assoc):
        """Sum of all distances to nearest medoid. Closes over: G.

        NB: identity is a measure of similarity, not distance, so costs will be
        negative numbers (pairwise cost = -identity)
        """
        cost = 0.0
        for med, cluster in assoc.items():
            for otr in cluster:
                cost -= G[otr][med]["weight"]
        return cost

    def same_clustering(old, new):
        """True if no node-medoid associations changed. Closes over: G"""
        return list(sorted(old.keys())) == list(sorted(new.keys()))

    # 1. Select k of the n nodes as medoids
    # NB: For efficiency, take "hubs" or random?
    # curr_meds = get_random()
    curr_meds = get_cluster(infname, N)

    for _i in range(MAX_ITERATE):
        # 2. Associate each "other" node to the closest medoid
        best_assoc = curr_assoc = associate(curr_meds)
        best_cost = total_cost(curr_assoc)
        best_meds = curr_meds

        # 3,4. For each medoid and each other node, replace the medoid with
        # each node and recalculate the cost of the new association.
        # Choose the lowest-cost association.
        for i, _med in enumerate(curr_meds):
            for otr in nodes:
                if otr in curr_meds:
                    continue
                new_meds = curr_meds[:i] + [otr] + curr_meds[i + 1 :]
                new_assoc = associate(new_meds)
                new_cost = total_cost(new_assoc)
                if new_cost < best_cost:
                    best_cost = new_cost
                    best_assoc = new_assoc
                    best_meds = new_meds

        # 5. Repeat steps 2-4 until medoids do not change
        if same_clustering(curr_assoc, best_assoc):
            break
        else:
            curr_assoc = best_assoc
            curr_assoc = best_assoc
            curr_meds = best_meds

    return curr_meds


def main(argv=None):
    args = sys.argv[1:] if argv is None else list(argv)
    if help_requested(args, __doc__):
        return 0

    try:
        option, N, infname = args[0], int(args[1]), args[2]
    except (IndexError, ValueError):
        sys.exit(__doc__)

    index = SeqIO.index(infname, "fasta")

    if option in ("-r", "--random"):
        accs = get_random(index, N)
    elif option in ("-f", "--first"):
        accs = get_first(infname, N)
    elif option in ("-l", "--longest"):
        accs = get_by_length(infname, N, True)
    elif option in ("-s", "--shortest"):
        accs = get_by_length(infname, N, False)
    elif option in ("-c", "--cluster"):
        accs = get_cluster(infname, N)
    elif option in ("-k", "--k-means"):
        accs = get_kmeans(infname, N)
    else:
        sys.exit(__doc__)

    # SeqIO.index returns each record's text from the file as bytes
    for acc in accs:
        block = index.get_raw(acc)
        sys.stdout.buffer.write(block)
