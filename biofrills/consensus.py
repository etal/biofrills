"""Calculate the consensus of a sequence alignment."""

import math
import logging
from collections import Counter, defaultdict

from biofrills import alnutils, sequtils


def consensus(aln, weights=None, gap_threshold=0.5, simple=False, trim_ends=True):
    """Get the consensus of an alignment, as a string.

    Emit gap characters for majority-gap columns; apply various strategies to
    choose the consensus amino acid type for the remaining columns.

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
        MAPGAPS).
    gap_threshold : float
        If the proportion of gap characters in a column is greater than or equal
        to this value (after sequence weighting, if applicable), then the
        consensus character emitted will be a gap instead of an amino acid type.

    """
    # Choose your algorithms!
    if simple:
        # Use the simple, unweighted algorithm
        col_consensus = make_simple_col_consensus(alnutils.aa_frequencies(aln))
        def is_majority_gap(col):
            return (float(col.count('-')) / len(col) >= gap_threshold)
        # ENH (alternatively/additionally): does any aa occur more than once?
        # ENH: choose gap-decisionmaking separately from col_consensus
    else:
        # Use the entropy-based, weighted algorithm
        if weights is None:
            seq_weights = alnutils.sequence_weights(aln, 'avg1')
        else:
            seq_weights = weights
        aa_frequencies = alnutils.aa_frequencies(aln, weights=seq_weights)
        col_consensus = make_entropy_col_consensus(aa_frequencies)
        def is_majority_gap(col):
            gap_count = 0.0
            for wt, char in zip(seq_weights, col):
                if char == '-':
                    gap_count += wt
            return (gap_count / sum(seq_weights) >= gap_threshold)

    # Traverse the alignment, handling gaps etc.
    def col_wise_consensus(columns):
        """Calculate the consensus chars for an iterable of columns."""
        if not trim_ends:
            # Track if we're in the N-term or C-term end of the sequence
            in_left_end = True
            maybe_right_tail = []
        # prev_col = None
        # prev_char = None
        for col in columns:
            # Lowercase cols mean explicitly, "don't include in consensus"
            if all(c.islower() for c in col if c not in '.-'):
                yield '-'
                continue
            if any(c.islower() for c in col):
                logging.warn('Mixed lowercase and uppercase letters in a '
                        'column: ' + ''.join(col))
                col = map(str.upper, col)

            # Gap chars
            is_gap = is_majority_gap(col)
            if not trim_ends:
                # Avoid N-terminal gaps in the consensus sequence
                if in_left_end:
                    if not is_gap:
                        # Match -- we're no longer in the left end
                        in_left_end = False
                    is_gap = False

            # When to yield a gap here:
            #   -----------     ---------   ------  ----------
            #   in_left_end     trim_ends   is_gap  yield gap?
            #   -----------     ---------   ------  ----------
            #   True            True        (True)  yes
            #   True            False       (False) (no -- def. char)
            #   False           True        T/F     yes, if is_gap
            #   False           False       (T/F)   NO! use maybe_right_tail
            #   -----------     ---------   ------  ----------

            if is_gap and trim_ends:
                yield '-'
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
                    for char in maybe_right_tail:
                        yield '-'
                    maybe_right_tail = []
                    yield cons_char

            # prev_col = col
            # prev_char = cons_char

        # Finally, if we were keeping a right (C-term) tail, emit it
        if not trim_ends:
            for char in maybe_right_tail:
                yield char

    return ''.join(col_wise_consensus(zip(*aln)))


# --- Algorithm: entropy ---

def entropy_func(f_freqs, b_freqs):
    def entropy_term(aa):
        f = f_freqs[aa] 
        b = b_freqs[aa]
        return f * math.log(f/b)
    return entropy_term


def make_entropy_col_consensus(bg_freqs):
    """Consensus according to maximal relative entropy term (MET).

    For a given column i, choose the residue j with the highest relative
    entropy term::

        f_ij ln(f_ij/b_j)

    where f_ij = column aa frequency, b_j = background aa frequency.

    Source: http://bioinformatics.oxfordjournals.org/content/24/18/1987.long
    """
    def col_consensus(col):
        col_freqs = sequtils.aa_frequencies(col)
        entroper = entropy_func(col_freqs, bg_freqs)
        return max(col_freqs.keys(), key=entroper)
    return col_consensus


# --- Algorithm: simple ---
# Ironically, this heuristic approach takes more code (and probably time) than
# the entropy-based calculation...

def make_simple_col_consensus(bg_freqs):
    """Consensus by simple plurality, unweighted.

    Resolves ties by two heuristics:

    1. Prefer the aa that follows the preceding consensus aa type most often
       in the original sequences.
    2. Finally, prefer the less-common aa type.
    """
    # Hack: use default kwargs to persist across iterations
    def col_consensus(col, prev_col=[], prev_char=[]):
        # Count the amino acid types in this column
        aa_counts = sequtils.aa_frequencies(col)
        assert aa_counts, "Column is all gaps! That's not allowed."
        # Take the most common residue(s)
        best_char, best_score = max(aa_counts.iteritems(),
                                    key=lambda kv: kv[1])
        # Resolve ties
        ties = [aa for aa in aa_counts if aa_counts[aa] == best_score]
        if len(ties) > 1:
            # Breaker #1: most common after the prev. consensus char
            # Resolve a tied col by restricting to rows where the preceding
            # char is the consensus type for that (preceding) col
            if prev_char and prev_col:
                mc_next = Counter(
                        [b for a, b in zip(prev_col, col)
                            if a == prev_char[0] and b in ties]
                        ).most_common()
                ties_next = [x[0] for x in mc_next
                        if x[1] == mc_next[0][1]]
                if ties_next:
                    ties = ties_next
            if len(ties) > 1:
                # Breaker #2: lowest overall residue frequency
                ties.sort(key=lambda aa: bg_freqs[aa])
            best_char = ties[0]
        else:
            assert best_char == ties[0], \
                    'WTF %s != %s[0]' % (best_char, ties)
        # Save values for tie-breaker #1
        prev_col[:] = col
        prev_char[:] = best_char
        return best_char
    return col_consensus


# --- XXX obsolete stuff ---

# XXX only used for CHAIN -- omit? or, option to return all ties?
def supported(aln):
    """Get only the supported consensus residues in each column.

    Meaning:
    - Omit majority-gap columns
    - Omit columns where no residue type appears more than once
    - In case of a tie, return all the top-scoring residue types
      (no prioritization)

    Returns a *list* -- not a string! -- where elements are strings of the
    consensus character(s), potentially a gap ('-') or multiple chars ('KR').
    """
    def col_consensus(columns):
        """Calculate the consensus chars for an iterable of columns."""
        for col in columns:
            if (# Majority gap chars
                (col.count('-') >= len(col)/2) or
                # Lowercase cols mean "don't include in consensus"
                all(c.islower() for c in col if c not in '.-')
                ):
                yield '-'
                continue
            # Validation - copied from consensus() above
            if any(c.islower() for c in col):
                logging.warn('Mixed lowercase and uppercase letters in a '
                        'column: ' + ''.join(col))
                col = map(str.upper, col)
            # Calculate the consensus character
            most_common = Counter(
                    [c for c in col if c not in '-']
                    ).most_common()
            if not most_common:
                # XXX ever reached?
                logging.warn("Column is all gaps! How did that happen?")
            if most_common[0][1] == 1:
                # No char has frequency > 1; no consensus char
                yield '-'
            elif (len(most_common) > 1 and
                    most_common[0][1] == most_common[1][1]):
                # Tie for most-common residue type
                ties = [x[0] for x in most_common
                        if x[1] == most_common[0][1]]
                yield ''.join(ties)
            else:
                yield most_common[0][0]

    return list(col_consensus(zip(*aln)))

