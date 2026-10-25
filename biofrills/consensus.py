"""Calculate the consensus of a sequence alignment."""

import logging
import math
from collections import Counter

from biofrills import alnutils, sequtils

logger = logging.getLogger(__name__)


def consensus(aln, weights=None, gap_threshold=0.5, simple=False, trim_ends=True):
    """Get the consensus of an alignment, as a string.

    Emit gap characters for majority-gap columns; apply various strategies to
    choose the consensus amino acid type for the remaining columns. Ties
    between amino acid types are broken by explicit rules, ending in
    alphabetical order, so the choice among tied types does not depend on the
    order in which residue types were counted; see make_entropy_col_consensus
    and make_simple_col_consensus.

    An alignment whose rows differ in length, or that contains no residue at
    all, raises ValueError.

    Parameters
    ----------

    simple : bool
        If True, use simple plurality to determine the consensus amino acid
        type, without weighting sequences for similarity. Otherwise, weight
        sequences for similarity and use relative entropy to choose the
        consensus amino acid type.
    weights : dict or None
        Sequence weights. If given, used to calculate amino acid frequencies;
        otherwise calculated within this function (i.e. this is a way to speed
        up the function if sequence weights have already been calculated).
        Ignored in 'simple' mode.
    trim_ends : bool
        If False, stretch the consensus sequence to include the N- and C-tails
        of the alignment, even if those flanking columns are mostly gap
        characters. This avoids terminal gaps in the consensus (needed for
        MAPGAPS). A column with no residue to report, either all gaps or all
        lowercase insert residues, still yields a gap character: the
        consensus has exactly one character per alignment column, in column
        order.
    gap_threshold : float
        If the proportion of gap characters in a column is greater than or equal
        to this value (after sequence weighting, if applicable), then the
        consensus character emitted will be a gap instead of an amino acid type.

    """
    # A mixed-case column is uppercased below, so the background frequencies
    # fold case too; otherwise a residue type that occurs elsewhere only in
    # lowercase would have no background frequency.
    upper_rows = [str(rec.seq).upper() for rec in aln]
    # Choose your algorithms!
    if simple:
        # Use the simple, unweighted algorithm
        col_consensus = make_simple_col_consensus(
            sequtils.aa_frequencies("".join(upper_rows))
        )

        def is_majority_gap(col):
            return col.count("-") / len(col) >= gap_threshold

        # ENH (alternatively/additionally): does any aa occur more than once?
        # ENH: choose gap-decisionmaking separately from col_consensus
    else:
        # Use the entropy-based, weighted algorithm
        if weights is None:
            seq_weights = alnutils.sequence_weights(aln, "avg1")
        else:
            seq_weights = weights
        aa_frequencies = alnutils.aa_frequencies(upper_rows, weights=seq_weights)
        col_consensus = make_entropy_col_consensus(aa_frequencies)

        def is_majority_gap(col):
            gap_count = 0.0
            for wt, char in zip(seq_weights, col, strict=True):
                if char == "-":
                    gap_count += wt
            return gap_count / sum(seq_weights) >= gap_threshold

    # Traverse the alignment, handling gaps etc.
    def col_wise_consensus(columns):
        """Calculate the consensus chars for an iterable of columns."""
        # Track if we're in the N-term or C-term end of the sequence; these
        # change only when trim_ends is False
        in_left_end = True
        maybe_right_tail = []
        for col in columns:
            # Lowercase cols mean explicitly, "don't include in consensus";
            # a column of only gap characters has no residue to report either.
            # Its gap goes behind any C-terminal residues already buffered, so
            # that the consensus stays aligned with the columns.
            if all(c.islower() for c in col if c not in ".-"):
                if trim_ends or in_left_end:
                    yield "-"
                else:
                    maybe_right_tail.append("-")
                continue
            if any(c.islower() for c in col):
                logger.warning(
                    "Mixed lowercase and uppercase letters in a column: %s",
                    "".join(col),
                )
                col = tuple(map(str.upper, col))

            # Gap chars
            is_gap = is_majority_gap(col)
            if not trim_ends:
                # Avoid N-terminal gaps in the consensus sequence
                if in_left_end:
                    if not is_gap:
                        # Match -- we're no longer in the left end
                        in_left_end = False
                    is_gap = False

            # When to yield a gap for a column with a residue to report (a
            # column without one is handled above):
            #   -----------     ---------   ------  ----------
            #   in_left_end     trim_ends   is_gap  yield gap?
            #   -----------     ---------   ------  ----------
            #   True            True        (True)  yes
            #   True            False       (False) (no -- def. char)
            #   False           True        T/F     yes, if is_gap
            #   False           False       (T/F)   NO! use maybe_right_tail
            #   -----------     ---------   ------  ----------

            if is_gap and trim_ends:
                yield "-"
                continue

            # Get the consensus character, using the chosen algorithm
            cons_char = col_consensus(col)

            if trim_ends:
                yield cons_char
            else:
                # Avoid C-terminal gaps in the consensus sequence
                if is_gap:
                    maybe_right_tail.append(cons_char)
                else:
                    # Match -> gaps weren't the right tail; emit all gaps
                    yield from "-" * len(maybe_right_tail)
                    maybe_right_tail = []
                    yield cons_char

        # Finally, emit any C-terminal tail kept when trim_ends is False
        yield from maybe_right_tail

    return "".join(col_wise_consensus(zip(*aln, strict=True)))


# --- Algorithm: entropy ---


def entropy_func(f_freqs, b_freqs):
    def entropy_term(aa):
        f = f_freqs[aa]
        b = b_freqs[aa]
        return f * math.log(f / b)

    return entropy_term


def make_entropy_col_consensus(bg_freqs):
    """Consensus according to maximal relative entropy term (MET).

    For a given column i, choose the residue j with the highest relative
    entropy term::

        f_ij ln(f_ij/b_j)

    where f_ij = column aa frequency, b_j = background aa frequency.

    Exact ties go to the alphabetically first residue type, so the choice
    among them does not depend on the order in which residue types were
    counted.

    Source: Przybylski D, Rost B (2008). Powerful fusion: PSI-BLAST and
    consensus sequences. Bioinformatics 24(18):1987-1993.
    https://doi.org/10.1093/bioinformatics/btn384
    """

    def col_consensus(col):
        col_freqs = sequtils.aa_frequencies(col)
        entroper = entropy_func(col_freqs, bg_freqs)
        return min(col_freqs, key=lambda aa: (-entroper(aa), aa))

    return col_consensus


# --- Algorithm: simple ---
# Ironically, this heuristic approach takes more code (and probably time) than
# the entropy-based calculation...


def make_simple_col_consensus(bg_freqs):
    """Consensus by simple plurality, unweighted.

    Resolves ties by three rules, in order:

    1. Prefer the aa that follows the preceding consensus aa type most often
       in the original sequences. The preceding consensus aa type is the one
       chosen for the most recent column passed to the returned function;
       columns that consensus() emits as gaps without choosing a residue do
       not count.
    2. Prefer the less-common aa type.
    3. Finally, prefer the alphabetically first aa type, so that the result
       never depends on the order in which residue types were counted.
    """
    # The column and residue chosen by the previous call, for tie-breaker #1
    prev_col = None
    prev_char = None

    def col_consensus(col):
        nonlocal prev_col, prev_char
        # Count the amino acid types in this column
        aa_counts = sequtils.aa_frequencies(col)
        # Take the most common residue(s)
        best_score = max(aa_counts.values())
        ties = [aa for aa in aa_counts if aa_counts[aa] == best_score]
        # Resolve ties
        if len(ties) > 1:
            # Breaker #1: most common after the prev. consensus char
            # Resolve a tied col by restricting to rows where the preceding
            # char is the consensus type for that (preceding) col
            if prev_char is not None:
                mc_next = Counter(
                    [
                        b
                        for a, b in zip(prev_col, col, strict=True)
                        if a == prev_char and b in ties
                    ]
                ).most_common()
                ties_next = [x[0] for x in mc_next if x[1] == mc_next[0][1]]
                if ties_next:
                    ties = ties_next
            if len(ties) > 1:
                # Breaker #2: lowest overall residue frequency;
                # breaker #3: alphabetical order
                ties.sort(key=lambda aa: (bg_freqs[aa], aa))
        best_char = ties[0]
        # Save values for tie-breaker #1
        prev_col = tuple(col)
        prev_char = best_char
        return best_char

    return col_consensus
