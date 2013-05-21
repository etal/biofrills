"""Utilities for manipulating multiple sequence alignments.
"""

from collections import Counter, defaultdict
from copy import deepcopy

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def aa_counts(aln, weights=None, gap_chars='-.'):
    """Calculate the amino acid frequencies in a set of SeqRecords.

    Weights for each sequence in the alignment can be given as a list/tuple,
    usually calculated with the sequence_weights function. For convenience, you
    can also pass "weights=True" and the weights will be calculated with
    sequence_weights here.
    """
    if weights is None:
        counts = Counter()
        for rec in aln:
            seq_counts = Counter(str(rec.seq))
            counts.update(seq_counts)
    else:
        if weights == True:
            # For convenience
            weights = sequence_weights(aln)
        else:
            assert len(weights) == len(aln), (
                "Length mismatch: weights = %d, alignment = %d"
                % (len(weights), len(aln)))
        counts = defaultdict(float)
        for col in zip(*aln):
            for aa, wt in zip(col, weights):
                counts[aa] += wt

    # Don't count gaps
    for gap_char in gap_chars:
        if gap_char in counts:
            del counts[gap_char]
    return counts


def aa_frequencies(aln, weights=None, gap_chars='-.'):
    """Frequency of each residue type in an alignment.

    Alignment is a MultipleSeqAlignment or iterable of SeqRecords.
    """
    counts = aa_counts(aln, weights, gap_chars)
    # Reduce to frequencies
    scale = 1.0 / sum(counts.values())
    return dict((aa, cnt * scale) for aa, cnt in counts.iteritems())


def blocks(aln, threshold=0.5, weights=None):
    """Remove gappy columns from an alignment."""
    assert len(aln)
    if weights == False:
        def pct_nongaps(col):
            return 1 - (float(col.count('-')) / len(col))
    else:
        if weights in (None, True):
            weights = sequence_weights(aln, 'avg1')
        def pct_nongaps(col):
            assert len(col) == len(weights)
            ngaps = sum(wt * (c == '-')
                        for wt, c in zip(weights, col))
            return 1 - (ngaps / len(col))

    seqstrs = [str(rec.seq) for rec in aln]
    clean_cols = [col for col in zip(*seqstrs)
                  if pct_nongaps(col) >= threshold]
    alphabet = aln[0].seq.alphabet
    clean_seqs = [Seq(''.join(row), alphabet)
                  for row in zip(*clean_cols)]
    clean_recs = []
    for rec, seq in zip(aln, clean_seqs):
        newrec = deepcopy(rec)
        newrec.seq = seq
        clean_recs.append(newrec)
    return MultipleSeqAlignment(clean_recs, alphabet=alphabet)


def col_counts(col, weights=None, gap_chars='-.'):
    """Absolute counts of each residue type in a single column."""
    cnt = defaultdict(float)
    for aa, wt in zip(col, weights):
        if aa not in gap_chars:
            cnt[aa] += wt
    return cnt


def col_frequencies(col, weights=None, gap_chars='-.'):
    """Frequencies of each residue type (totaling 1.0) in a single column."""
    counts = col_counts(col, weights, gap_chars)
    # Reduce to frequencies
    scale = 1.0 / sum(counts.values())
    return dict((aa, cnt * scale) for aa, cnt in counts.iteritems())


def remove_empty_cols(records):
    """Remove all-gap columns from aligned SeqRecords."""
    # In case it's a generator, turn it into a list
    records = list(records)
    seqstrs = [str(rec.seq) for rec in records]
    clean_cols = [col
                  for col in zip(*seqstrs)
                  if not all(c == '-' for c in col)]
    clean_seqs = [''.join(row)
                  for row in zip(*clean_cols)]
    for rec, clean_seq in zip(records, clean_seqs):
        yield SeqRecord(Seq(clean_seq, rec.seq.alphabet), id=rec.id,
                        name=rec.name, description=rec.description,
                        dbxrefs=rec.dbxrefs, features=rec.features,
                        annotations=rec.annotations,
                        letter_annotations=rec.letter_annotations)


def sequence_weights(aln, scaling='none', gap_chars='-.'):
    """Weight aligned sequences to emphasize more divergent members.

    Returns a list of floating-point numbers between 0 and 1, corresponding to
    the proportional weight of each sequence in the alignment. The first list
    is the weight of the first sequence in the alignment, and so on.

    Scaling schemes:
        - 'sum1': Weights sum to 1.0.
        - 'max1': Weights are all scaled so the max is 1.0.
        - 'avg1': Average (mean) weight is 1.0.
        - 'andy': Average (mean) weight is 0.5, ceiling is 1.0.
        - 'none': Weights are scaled to sum to the effective number of
                  independent sequences.

    Method: At each column position, award each different residue an equal
    share of the weight, and then divide that weight equally among the
    sequences sharing the same residue.  For each sequence, sum the
    contributions from each position to give a sequence weight.

    See Henikoff & Henikoff (1994): Position-based sequence weights.
    """
    # Probability is hard, let's estimate by sampling!
    # Sample k from a population of 20 with replacement; how many unique k were
    # chosen? Average of 10000 runs for k = 0..100
    expectk = [0.0, 1.0, 1.953, 2.861, 3.705, 4.524, 5.304, 6.026, 6.724, 7.397,
               8.04, 8.622, 9.191, 9.739, 10.264, 10.758, 11.194, 11.635,
               12.049, 12.468, 12.806, 13.185, 13.539, 13.863, 14.177, 14.466,
               14.737, 15.005, 15.245, 15.491, 15.681, 15.916, 16.12, 16.301,
               16.485, 16.671, 16.831, 16.979, 17.151, 17.315, 17.427, 17.559,
               17.68, 17.791, 17.914, 18.009, 18.113, 18.203, 18.298, 18.391,
               18.46, 18.547, 18.617, 18.669, 18.77, 18.806, 18.858, 18.934,
               18.978, 19.027, 19.085, 19.119, 19.169, 19.202, 19.256, 19.291,
               19.311, 19.357, 19.399, 19.416, 19.456, 19.469, 19.5, 19.53,
               19.553, 19.562, 19.602, 19.608, 19.629, 19.655, 19.67, 19.681,
               19.7, 19.716, 19.724, 19.748, 19.758, 19.765, 19.782, 19.791,
               19.799, 19.812, 19.82, 19.828, 19.844, 19.846, 19.858, 19.863,
               19.862, 19.871, 19.882]

    def col_weight(column):
        """Represent the diversity at a position.

        Award each different residue an equal share of the weight, and then
        divide that weight equally among the sequences sharing the same
        residue.

        So, if in a position of a multiple alignment, r different residues
        are represented, a residue represented in only one sequence contributes
        a score of 1/r to that sequence, whereas a residue represented in s
        sequences contributes a score of 1/rs to each of the s sequences.
        """
        # Skip columns of all or mostly gaps (i.e. rare inserts)
        min_nongap = max(2, .2*len(column))
        if len([c for c in column if c not in gap_chars]) < min_nongap:
            return ([0] * len(column), 0)
        # Count the number of occurrences of each residue type
        # (Treat gaps as a separate, 21st character)
        counts = Counter(column)
        # Get residue weights: 1/rs, where
        # r = nb. residue types, s = count of a particular residue type
        n_residues = len(counts)    # r
        freqs = dict((aa, 1.0 / (n_residues * count))
                     for aa, count in counts.iteritems())
        weights = [freqs[aa] for aa in column]
        return (weights, n_residues)

    seq_weights = [0] * len(aln)
    tot_nres = 0.0  # Expected no. different types in independent seqs
    # Sum the contributions from each position along each sequence
    # -> total weight
    for col in zip(*aln):
        wts, nres = col_weight(col)
        assert sum(wts) <= 20
        tot_nres += expectk[nres] if nres < len(expectk) else 20
        for idx, wt in enumerate(wts):
            seq_weights[idx] += wt
    # if tot_nres == 0:
    #     raise ValueError("Alignment has no meaningful columns to weight")
    # Normalize w/ the given scaling criterion
    if scaling == 'none':
        avg_seq_len = tot_nres / len(aln)
        return [wt/avg_seq_len for wt in seq_weights]
    if scaling == 'max1':
        scale = 1.0 / max(seq_weights)
    elif scaling == 'sum1':
        scale = 1.0 / sum(seq_weights)
    elif scaling == 'avg1':
        scale = len(aln) / sum(seq_weights)
    elif scaling == 'andy':
        # "Robust" strategy used in CHAIN (Neuwald 2003)
        scale = len(aln) / sum(seq_weights)
        return [min(scale * wt, 1.0) for wt in seq_weights]
    else:
        raise ValueError("Unknown scaling scheme '%s'" % scaling)
    return [scale * wt for wt in seq_weights]


def to_graph(alnfname, weight_func):
    """Create a NetworkX graph from a sequence alignment.

    Nodes are string sequence IDs; edge weights are the output of weight_func
    between each pair, by default the absolute identity (# identical chars).
    """
    import networkx
    G = networkx.Graph()
    aln = AlignIO.read(alnfname, 'fasta')
    for i, arec in enumerate(aln):
        for brec in aln[i+1:]:
            ident = weight_func(str(arec.seq), str(brec.seq))
            G.add_edge(arec.id, brec.id, weight=ident)
    return G

