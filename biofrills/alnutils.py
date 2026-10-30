"""Utilities for manipulating multiple sequence alignments."""

import math
from collections import Counter, defaultdict
from copy import deepcopy

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# A weighted proportion at most this far below its threshold counts as
# reaching it; blocks and consensus.consensus both apply it.  Sequence weights
# are rounded floats, so a column whose exact weighted proportion equals the
# threshold can compute a few units in the last place below it.  In measured
# alignments that error was at most 7e-16, and no proportion not equal to the
# threshold came closer to it than 1.2e-6.
_GAP_TOLERANCE = 1e-12


def aa_counts(aln, weights=None, gap_chars="-."):
    """Calculate the amino acid frequencies in a set of SeqRecords.

    Weights for each sequence in the alignment can be given as a list/tuple,
    usually calculated with the sequence_weights function. For convenience, you
    can also pass "weights=True" and the weights will be calculated with
    sequence_weights here. With weights, the rows may also be given as
    sequence strings; rows of different lengths then raise ValueError.
    """
    if weights is None:
        counts = Counter()
        for rec in aln:
            seq_counts = Counter(str(rec.seq))
            counts.update(seq_counts)
    else:
        if weights is True:
            # For convenience
            weights = sequence_weights(aln)
        else:
            assert len(weights) == len(aln), (
                f"Length mismatch: weights = {len(weights)}, alignment = {len(aln)}"
            )
        counts = defaultdict(float)
        for col in zip(*aln, strict=True):
            for aa, wt in zip(col, weights, strict=True):
                counts[aa] += wt

    # Don't count gaps
    for gap_char in gap_chars:
        if gap_char in counts:
            del counts[gap_char]
    return counts


def aa_frequencies(aln, weights=None, gap_chars="-."):
    """Frequency of each residue type in an alignment.

    Alignment is a MultipleSeqAlignment or iterable of SeqRecords. Raises
    ValueError if the alignment contains no residues, only gap characters.
    """
    counts = aa_counts(aln, weights, gap_chars)
    return _frequencies(counts)


def _frequencies(counts):
    """Scale residue counts to frequencies that total 1.0."""
    if not counts:
        raise ValueError("No residues to count: every character is a gap")
    scale = 1.0 / sum(counts.values())
    return dict((aa, cnt * scale) for aa, cnt in counts.items())


def blocks(aln, threshold=0.5, weights=None):
    """Remove gappy columns from an alignment.

    Keep each column whose proportion of non-gap characters is at least
    threshold; only '-' counts as a gap. Every record is kept; its sequence is
    empty if no column remains. Rows of different lengths, and weights of a
    length other than the number of rows, raise ValueError.

    For given sequence weights, which columns are kept does not depend on the
    order of the alignment's rows: the unweighted proportion is a quotient of
    integer counts, and the weighted proportion a quotient of correctly
    rounded sums (math.fsum), both of which a permutation of the rows leaves
    unchanged.

    Parameters
    ----------

    threshold : float
        Keep a column whose proportion of non-gap characters is greater than
        or equal to this value (after sequence weighting, if applicable). A
        weighted proportion at most 1e-12 below the threshold counts as equal
        to it, so that rounding error in the sequence weights does not decide
        a column whose exact weighted proportion equals the threshold.
    weights : sequence of float, None, True or False
        One weight per row of the alignment, in row order, as
        sequence_weights returns. None or True calculates "avg1" weights
        within this function; False weights the rows equally. A weighted
        proportion is a share of the total weight, so weights that differ by
        a constant factor, as most scalings of sequence_weights do, keep the
        same columns.

    """
    assert len(aln)
    if weights is False:

        def has_enough_residues(col):
            return (len(col) - col.count("-")) / len(col) >= threshold
    else:
        if weights is None or weights is True:
            weights = sequence_weights(aln, "avg1")
        elif len(weights) != len(aln):
            raise ValueError(
                f"Length mismatch: weights = {len(weights)}, alignment = {len(aln)}"
            )
        total_weight = math.fsum(weights)
        min_nongap_fraction = threshold - _GAP_TOLERANCE

        def has_enough_residues(col):
            nongap_weight = math.fsum(
                wt for wt, char in zip(weights, col, strict=True) if char != "-"
            )
            return nongap_weight / total_weight >= min_nongap_fraction

    seqstrs = [str(rec.seq) for rec in aln]
    clean_cols = [col for col in zip(*seqstrs, strict=True) if has_enough_residues(col)]
    clean_recs = []
    for rec, seq in zip(aln, _rows(clean_cols, len(seqstrs)), strict=True):
        newrec = deepcopy(rec)
        newrec.seq = Seq(seq)
        clean_recs.append(newrec)
    return MultipleSeqAlignment(clean_recs)


def _rows(columns, nrows):
    """Join the given columns into nrows sequence strings.

    With no columns, the result is nrows empty strings.
    """
    return ["".join(col[idx] for col in columns) for idx in range(nrows)]


def col_counts(col, weights=None, gap_chars="-."):
    """Absolute counts of each residue type in a single column.

    Without weights, each residue counts 1.0.
    """
    if weights is None:
        weights = [1] * len(col)
    cnt = defaultdict(float)
    for aa, wt in zip(col, weights, strict=True):
        if aa not in gap_chars:
            cnt[aa] += wt
    return cnt


def col_frequencies(col, weights=None, gap_chars="-."):
    """Frequencies of each residue type (totaling 1.0) in a single column.

    Raises ValueError if the column holds only gap characters.
    """
    return _frequencies(col_counts(col, weights, gap_chars))


def remove_empty_cols(records):
    """Remove all-gap columns from aligned SeqRecords.

    Every record is kept; its sequence is empty if no column remains. Rows of
    different lengths raise ValueError.
    """
    # In case it's a generator, turn it into a list
    records = list(records)
    seqstrs = [str(rec.seq) for rec in records]
    clean_cols = [
        col for col in zip(*seqstrs, strict=True) if not all(c == "-" for c in col)
    ]
    for rec, clean_seq in zip(records, _rows(clean_cols, len(seqstrs)), strict=True):
        yield SeqRecord(
            Seq(clean_seq),
            id=rec.id,
            name=rec.name,
            description=rec.description,
            dbxrefs=rec.dbxrefs,
            features=rec.features,
            annotations=rec.annotations,
            letter_annotations=rec.letter_annotations,
        )


def _effective_count(mean_types, n_seqs):
    """Number of independent sequences that show mean_types distinct types.

    Solves 20 * (1 - (19/20)**k) = mean_types for k, where the left side is
    the expected number of distinct types among k residues drawn uniformly
    from the 20 amino acids, and caps the result at n_seqs.  A mean of 20 or
    more has no solution and gives n_seqs.
    """
    if mean_types >= 20:
        return float(n_seqs)
    return min(math.log1p(-mean_types / 20) / math.log1p(-1 / 20), float(n_seqs))


def sequence_weights(aln, scaling="none", gap_chars="-."):
    """Weight aligned sequences to emphasize more divergent members.

    Returns a list of non-negative floating-point numbers, one per sequence
    of the alignment, in the same order.

    Scaling schemes:
        - 'sum1': Weights sum to 1.0.
        - 'max1': Weights are all scaled so the max is 1.0.
        - 'avg1': Average (mean) weight is 1.0.
        - 'andy': The robust weights of CHAIN (Neuwald et al. 2003): the
                  weights are scaled to a mean of one-half, and the
                  occasional weight above 1.0 is then truncated to 1.0, so
                  the mean of the result is at most one-half.
        - 'none': Weights sum to the effective number of independent
                  sequences, estimated as described below.

    If no column holds enough residues to tell the sequences apart, as in a
    single sequence, every sequence receives the same weight: 1.0, or 1/N
    under 'sum1', or 0.5 under 'andy'. An alignment of no sequences, or of
    rows of different lengths, raises ValueError.

    Method: At each column position, award each different residue an equal
    share of the weight, and then divide that weight equally among the
    sequences sharing the same residue.  For each sequence, sum the
    contributions from each position to give a sequence weight.  Each gap
    character counts as a residue type of its own.  Columns in which fewer
    than two residues, or fewer than a fifth of the rows, are not gaps
    contribute nothing.

    See Henikoff & Henikoff (1994): Position-based sequence weights.

    Effective number of sequences: k independent sequences whose residues
    are drawn uniformly from the 20 amino acids show on average
    20 * (1 - (19/20)**k) distinct residue types in a column.  The 'none'
    scaling counts the distinct residue types of each contributing column,
    not counting gaps, and solves that expression for k at the mean count,
    capped at the number of sequences.  Identical sequences thus count as
    one sequence, and a mean of 20 or more types counts as all of them.
    Real residue frequencies are not uniform, so independent sequences show
    fewer distinct types than this expectation, and the estimate tends to
    be low.  PSI-BLAST (Altschul et al. 1997) likewise measures the
    diversity of an alignment by its mean number of distinct residue types
    per column.
    """

    def col_weight(column):
        """Represent the diversity at a position.

        Award each different residue an equal share of the weight, and then
        divide that weight equally among the sequences sharing the same
        residue.

        So, if in a position of a multiple alignment, r different residues
        are represented, a residue represented in only one sequence contributes
        a score of 1/r to that sequence, whereas a residue represented in s
        sequences contributes a score of 1/rs to each of the s sequences.

        Also returns the number of distinct residue types other than gaps,
        which is 0 for a skipped column and at least 1 otherwise.
        """
        # Skip columns of all or mostly gaps (i.e. rare inserts)
        min_nongap = max(2, 0.2 * len(column))
        if len([c for c in column if c not in gap_chars]) < min_nongap:
            return ([0] * len(column), 0)
        # Count the number of occurrences of each residue type
        # (Each gap character is a residue type of its own here, but not in
        # the count of residue types returned for the effective number)
        counts = Counter(column)
        # Get residue weights: 1/rs, where
        # r = nb. residue types, s = count of a particular residue type
        n_residues = len(counts)  # r
        freqs = dict((aa, 1.0 / (n_residues * count)) for aa, count in counts.items())
        weights = [freqs[aa] for aa in column]
        return (weights, sum(1 for aa in counts if aa not in gap_chars))

    if not len(aln):
        raise ValueError("No sequences to weight")
    seq_weights = [0] * len(aln)
    n_cols = 0
    n_types = 0  # Distinct non-gap types, summed over contributing columns
    # Sum the contributions from each position along each sequence
    # -> total weight
    for col in zip(*aln, strict=True):
        wts, col_types = col_weight(col)
        if col_types:
            n_cols += 1
            n_types += col_types
        for idx, wt in enumerate(wts):
            seq_weights[idx] += wt
    if not n_cols:
        # Nothing distinguishes the sequences, so weight them equally: one
        # unit of raw weight each, as N independent sequences, which every
        # scaling below turns into its uniform vector.
        seq_weights = [1.0] * len(aln)
    # Normalize w/ the given scaling criterion
    if scaling == "none":
        if n_cols:
            # Each contributing column distributes one unit of weight in
            # total, so the raw weights sum to n_cols.
            scale = _effective_count(n_types / n_cols, len(aln)) / n_cols
        else:
            scale = 1.0
    elif scaling == "max1":
        scale = 1.0 / max(seq_weights)
    elif scaling == "sum1":
        scale = 1.0 / sum(seq_weights)
    elif scaling == "avg1":
        scale = len(aln) / sum(seq_weights)
    elif scaling == "andy":
        # The robust weights of CHAIN (Neuwald et al. 2003): normalized to a
        # mean of one-half, with weights above one truncated to one.
        scale = 0.5 * len(aln) / sum(seq_weights)
        return [min(scale * wt, 1.0) for wt in seq_weights]
    else:
        raise ValueError(f"Unknown scaling scheme '{scaling}'")
    return [scale * wt for wt in seq_weights]


def to_graph(alnfname, weight_func):
    """Create a NetworkX graph from a sequence alignment.

    Nodes are string sequence IDs; edge weights are the output of weight_func
    between each pair, by default the absolute identity (# identical chars).
    """
    import networkx

    G = networkx.Graph()
    aln = AlignIO.read(alnfname, "fasta")
    for i, arec in enumerate(aln):
        # Every record is a node, including a record that pairs with no other
        G.add_node(arec.id)
        for brec in aln[i + 1 :]:
            ident = weight_func(str(arec.seq), str(brec.seq))
            G.add_edge(arec.id, brec.id, weight=ident)
    return G
