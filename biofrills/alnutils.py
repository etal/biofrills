"""Utilities for manipulating multiple sequence alignments.
"""
import logging
from collections import Counter, defaultdict

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO


def aa_frequencies(aln, gap_chars='-.', weights=None):
    """Calculate the amino acid frequencies in a set of SeqRecords.

    Weights for each sequence in the alignment can be given as a list/tuple,
    usually calculated with the sequence_weights function. For convenience, you
    can also pass "weights=True" and the weights will be calculated with
    sequence_weights here.
    """
    if weights is None:
        aa_counts = Counter()
        for rec in aln:
            seq_counts = Counter(str(rec.seq))
            aa_counts.update(seq_counts)
    else:
        if weights == True:
            # For convenience
            weights = sequence_weights(aln)
        else:
            assert len(weights) == len(aln), (
                "Length mismatch: weights = %d, alignment = %d"
                % (len(weights), len(aln)))
        aa_counts = defaultdict(float)
        for col in zip(*aln):
            for aa, wt in zip(col, weights):
                aa_counts[aa] += wt

    # Don't count gaps
    for gap_char in gap_chars:
        if gap_char in aa_counts:
            del aa_counts[gap_char]
    # Reduce to frequencies
    scale = 1.0 / sum(aa_counts.values())
    return dict((aa, cnt * scale) for aa, cnt in aa_counts.iteritems())


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


def sequence_weights(aln, scaling='none'):
    """Weight aligned sequences to emphasize more divergent members.

    Returns a list of floating-point numbers between 0 and 1, corresponding to
    the proportional weight of each sequence in the alignment. The first list
    is the weight of the first sequence in the alignment, and so on.

    Scaling schemes:
        - 'sum1': Weights sum to 1.0.
        - 'max1': Weights are all scaled so the max is 1.0.
        - 'avg1': Average (mean) weight is 1.0.
        - 'andy': Average (mean) weight is 0.5, ceiling is 1.0.

    Method: At each column position, award each different residue an equal
    share of the weight, and then divide that weight equally among the
    sequences sharing the same residue.  For each sequence, sum the
    contributions from each position to give a sequence weight.

    See Henikoff & Henikoff (1994): Position-based sequence weights.
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
        """
        # Skip columns with all gaps or unique inserts
        if len([c for c in column if c not in '-.']) < 2:
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
    tot_nres = 0.0
    # Sum the contributions from each position along each sequence -> total weight
    for col in zip(*aln):
        wts, nres = col_weight(col)
        assert sum(wts) <= 20
        tot_nres += min(nres, 20)  # Limited alphabet even for independent seqs
        for idx, wt in enumerate(wts):
            seq_weights[idx] += wt
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

