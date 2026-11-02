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

# The gap characters of the alignment convention of HMMER (A2M and Stockholm
# output) and CHAIN (CMA files): an uppercase residue is aligned to a match
# column of a profile, a lowercase residue is an insertion relative to it,
# '-' is a deletion in a match column and '.' a gap in an insert column.
# blocks, remove_empty_cols and consensus.consensus apply the convention.
_GAP_CHARS = "-."


def _case_marks_inserts(rows):
    """Whether letter case in these aligned rows distinguishes inserts.

    Only an alignment with at least one uppercase residue follows the
    convention. One without, such as the lowercase nucleotide alignments that
    MAFFT writes, marks nothing by case, so its lowercase residues are
    aligned residues.
    """
    return any(row != row.lower() for row in rows)


def _has_aligned_residue(col, case_marks_inserts):
    """Whether any residue of an alignment column is an aligned residue.

    A column of gap characters alone has none, and where case marks inserts,
    neither does a column whose residues are all lowercase.
    """
    if case_marks_inserts:
        return any(not c.islower() for c in col if c not in _GAP_CHARS)
    return any(c not in _GAP_CHARS for c in col)


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
    threshold, where a column that holds no aligned residue has a proportion
    of zero. Both '-' and '.' count as gaps. Under the HMMER and CHAIN
    convention, a column whose residues are all lowercase is an insert column
    and holds no aligned residue; the convention applies when the alignment
    contains an uppercase residue, and an alignment without one is read
    without regard to case. In a column that holds an aligned residue,
    lowercase residues count toward the proportion.
    A column that is kept keeps its characters as they are, lowercase
    residues included. Every record is kept; its sequence is empty if no
    column remains. Per-letter annotations, such as FASTQ quality scores,
    keep their values at the kept columns; features are copied unchanged, so
    their locations still refer to the columns of the input. Rows of
    different lengths, and weights of a length other than the number of
    rows, raise ValueError.

    For given sequence weights, which columns are kept does not depend on the
    order of the alignment's rows: the unweighted proportion is a quotient of
    integer counts, and the weighted proportion a quotient of correctly
    rounded sums (math.fsum), both of which a permutation of the rows leaves
    unchanged, and whether a column holds an aligned residue does not depend
    on the order of its characters.

    Parameters
    ----------

    threshold : float
        The least proportion, as defined above, that keeps a column (after
        sequence weighting, if applicable). A weighted proportion at most
        1e-12 below the threshold counts as equal to it, so that rounding
        error in the sequence weights does not decide a column whose exact
        weighted proportion equals the threshold.
    weights : sequence of float, None, True or False
        One weight per row of the alignment, in row order, as
        sequence_weights returns. None or True calculates "avg1" weights
        within this function; False weights the rows equally. A weighted
        proportion is a share of the total weight, so weights that differ by
        a constant factor, as most scalings of sequence_weights do, keep the
        same columns.

    """
    assert len(aln)
    seqstrs = [str(rec.seq) for rec in aln]
    case_marks_inserts = _case_marks_inserts(seqstrs)
    if weights is False:

        def residue_fraction(col):
            nresidues = len(col) - sum(col.count(c) for c in _GAP_CHARS)
            return nresidues / len(col)

        min_fraction = threshold
    else:
        if weights is None or weights is True:
            weights = sequence_weights(aln, "avg1")
        elif len(weights) != len(aln):
            raise ValueError(
                f"Length mismatch: weights = {len(weights)}, alignment = {len(aln)}"
            )
        total_weight = math.fsum(weights)

        def residue_fraction(col):
            residue_weight = math.fsum(
                wt
                for wt, char in zip(weights, col, strict=True)
                if char not in _GAP_CHARS
            )
            return residue_weight / total_weight

        min_fraction = threshold - _GAP_TOLERANCE

    def has_enough_residues(col):
        if not _has_aligned_residue(col, case_marks_inserts):
            return 0.0 >= min_fraction
        return residue_fraction(col) >= min_fraction

    keep = [
        idx
        for idx, col in enumerate(zip(*seqstrs, strict=True))
        if has_enough_residues(col)
    ]
    clean_recs = []
    for rec, seqstr in zip(aln, seqstrs, strict=True):
        newrec = deepcopy(rec)
        # Selecting from the copy rather than from rec keeps every annotation
        # value of the result unshared with the input, also when no column is
        # dropped.
        letter_annotations = _select_letter_annotations(newrec, keep)
        # The seq setter refuses a new length while annotations are present.
        newrec.letter_annotations = {}
        newrec.seq = Seq(_select(seqstr, keep))
        newrec.letter_annotations = letter_annotations
        clean_recs.append(newrec)
    return MultipleSeqAlignment(clean_recs)


def _select(items, keep):
    """The items at the indices in keep, in a sequence of the same type.

    A string stays a string, and a list or tuple keeps its type; any other
    sequence, such as a NumPy array, is indexed by the list of indices.
    """
    if isinstance(items, str):
        return "".join([items[idx] for idx in keep])
    if isinstance(items, (list, tuple)):
        return type(items)(items[idx] for idx in keep)
    return items[keep]


def _select_letter_annotations(rec, keep):
    """The per-letter annotations of rec at the column indices in keep.

    When every column is kept, the annotations are returned unchanged,
    whatever the type of their values.
    """
    if len(keep) == len(rec):
        return dict(rec.letter_annotations)
    return {key: _select(value, keep) for key, value in rec.letter_annotations.items()}


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

    A column goes when every character is a gap, '-' and '.' alike. A column
    whose residues are all lowercase is kept: under the HMMER and CHAIN
    convention those residues are insertions rather than aligned residues,
    but they are residues of their sequences, and removing the column would
    delete them; blocks is the function that drops a column holding residues.

    Every record is kept; its sequence is empty if no column remains.
    Per-letter annotations, such as FASTQ quality scores, keep their values
    at the kept columns; features are passed on unchanged, so their locations
    still refer to the columns of the input. Rows of different lengths raise
    ValueError.
    """
    # In case it's a generator, turn it into a list
    records = list(records)
    seqstrs = [str(rec.seq) for rec in records]
    keep = [
        idx
        for idx, col in enumerate(zip(*seqstrs, strict=True))
        if not all(c in _GAP_CHARS for c in col)
    ]
    for rec, seqstr in zip(records, seqstrs, strict=True):
        yield SeqRecord(
            Seq(_select(seqstr, keep)),
            id=rec.id,
            name=rec.name,
            description=rec.description,
            dbxrefs=rec.dbxrefs,
            features=rec.features,
            annotations=rec.annotations,
            letter_annotations=_select_letter_annotations(rec, keep),
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

    Nodes are string sequence IDs; the weight of the edge between each pair is
    the output of weight_func on their sequence strings, for example
    pairutils.identity_abs, their number of matching sites.
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
